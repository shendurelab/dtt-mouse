"""Input/output helpers: fastq iteration, reference loading, pickle save/load."""
from __future__ import annotations
import gzip
import pickle
from pathlib import Path
from typing import Iterator


def iter_reads(path, limit: int | None = None) -> Iterator[bytes]:
    """Yield the sequence line (bytes) of each read in a (gzipped) fastq.

    Reads only the SEQ line of every 4-line record; quality is not needed for
    the array parse. Pass `limit` to sample the first N reads (PASS A).
    """
    n = 0
    with gzip.open(str(path), "rb") as fh:
        while True:
            if not fh.readline():          # header (or EOF)
                break
            seq = fh.readline().rstrip(b"\n")
            fh.readline()                  # '+'
            fh.readline()                  # qual
            yield seq
            n += 1
            if limit and n >= limit:
                break


def load_whitelist(path) -> list[bytes]:
    """Load a `<tag>.tapebc_whitelist.tsv` -> list of barcode byte-strings (in rank order)."""
    wl = []
    with open(path) as f:
        next(f)                            # header
        for line in f:
            wl.append(line.split("\t")[1].encode())
    return wl


def load_vocab(path) -> set[bytes]:
    """Load a `<tag>.nnn_vocab.tsv` -> set of in-vocab insertion byte-strings.

    Only rows flagged in_vocab == '1' are kept (the file also lists the full
    observed insertion spectrum for reference).
    """
    v = set()
    with open(path) as f:
        next(f)                            # header
        for line in f:
            c = line.split("\t")
            if c[2] == "1":
                v.add(c[0].encode())
    return v


def save_pickle(obj, path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump(obj, f)


def load_pickle(path):
    with open(path, "rb") as f:
        return pickle.load(f)
