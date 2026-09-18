"""Paired-end read merging for single-cell TAPE (wraps fastp).

seq3 is paired: R2 is the amplicon read (5' -> 3' from the integration barcode),
R1 is its low-quality mate from the 3'/poly-A end. We overlap-merge mates with
fastp (--merge --correction), which corrects low-quality overlap bases using the
higher-quality mate. Because R2 alone spans the full amplicon at high quality
while most R1 mates are uninformative, we use a merge-with-fallback scheme:

    * where a pair merges  -> use the merged read (R1-corrected, full length)
    * where it does not    -> use the R2 (amplicon) read

fastp emits merged reads in R1 orientation, so merged reads are reverse-
complemented back to the amplicon orientation. The output is a single fastq of
amplicon-oriented reads (headers preserved) ready for the shared SC parser.

Requires the input R1/R2 files to be in matching (paired) order.
"""
from __future__ import annotations
import gzip
import shutil
import subprocess
from pathlib import Path

from config import (FASTP_BIN, FASTP_OVERLAP_LEN_REQUIRE,
                    FASTP_OVERLAP_DIFF_LIMIT, FASTP_OVERLAP_DIFF_PCT_LIMIT)

_COMP = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")


def _revcomp(seq: bytes) -> bytes:
    return seq.translate(_COMP)[::-1]


def run_fastp(r1: Path, r2: Path, workdir: Path) -> dict:
    """Run fastp in merge mode. Returns paths to merged + unmerged outputs."""
    workdir.mkdir(parents=True, exist_ok=True)
    out = {
        "merged": workdir / "merged.fastq.gz",
        "un1": workdir / "unmerged.R1.fastq.gz",
        "un2": workdir / "unmerged.R2.fastq.gz",
        "json": workdir / "fastp.json",
        "html": workdir / "fastp.html",
    }
    cmd = [
        FASTP_BIN, "-i", str(r1), "-I", str(r2),
        "--merge", "--merged_out", str(out["merged"]),
        "--out1", str(out["un1"]), "--out2", str(out["un2"]),
        "--detect_adapter_for_pe", "--correction",
        "--disable_quality_filtering", "--disable_length_filtering",
        "--overlap_len_require", str(FASTP_OVERLAP_LEN_REQUIRE),
        "--overlap_diff_limit", str(FASTP_OVERLAP_DIFF_LIMIT),
        "--overlap_diff_percent_limit", str(FASTP_OVERLAP_DIFF_PCT_LIMIT),
        "-j", str(out["json"]), "-h", str(out["html"]),
        "--thread", "8",
    ]
    subprocess.run(cmd, check=True, capture_output=True)
    return out


def _iter_fastq(path: Path):
    """Iterate (header, seq, qual) from a gzipped fastq, decompressing with pigz
    (multi-core) when available, else Python gzip."""
    proc = None
    if shutil.which("pigz"):
        proc = subprocess.Popen(["pigz", "-dc", str(path)],
                                stdout=subprocess.PIPE, bufsize=1 << 22)
        fh = proc.stdout
    else:
        fh = gzip.open(str(path), "rb")
    try:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().rstrip(b"\n")
            fh.readline()
            q = fh.readline().rstrip(b"\n")
            yield h.rstrip(b"\n"), s, q
    finally:
        fh.close()
        if proc is not None:
            proc.wait()


def combine_oriented(fastp_out: dict, combined: Path, threads: int = 8) -> dict:
    """Write amplicon-oriented reads: reverse-complemented merged reads + the
    unmerged R2 (fallback). Compresses the output with pigz (multi-core) when
    available -- the single-threaded gzip write was the merge bottleneck. Content
    is identical regardless of the compressor. Returns counts."""
    n_merged = n_fallback = 0
    use_pigz = shutil.which("pigz") is not None
    wf = wp = None
    if use_pigz:
        wf = open(str(combined), "wb")
        wp = subprocess.Popen(["pigz", "-p", str(threads), "-c"],
                              stdin=subprocess.PIPE, stdout=wf, bufsize=1 << 22)
        w = wp.stdin
    else:
        w = gzip.open(str(combined), "wb")
    try:
        for h, s, q in _iter_fastq(fastp_out["merged"]):
            w.write(h + b"\n" + _revcomp(s) + b"\n+\n" + q[::-1] + b"\n")
            n_merged += 1
        for h, s, q in _iter_fastq(fastp_out["un2"]):
            w.write(h + b"\n" + s + b"\n+\n" + q + b"\n")
            n_fallback += 1
    finally:
        w.close()
        if wp is not None:
            wp.wait()
            wf.close()
    return {"merged": n_merged, "r2_fallback": n_fallback,
            "total": n_merged + n_fallback}


def merge_sample(r1: Path, r2: Path, workdir: Path, combined: Path) -> dict:
    """Full merge for one paired sample: fastp -> amplicon-oriented combined fastq."""
    fastp_out = run_fastp(r1, r2, workdir)
    stats = combine_oriented(fastp_out, combined)
    stats["merge_rate"] = stats["merged"] / max(1, stats["total"])
    return stats
