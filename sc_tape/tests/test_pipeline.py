#!/usr/bin/env python3
"""Single-cell circtape pipeline tests -- no raw data, no network, ~20 seconds.

Three groups:
  1. the per-read parse and the three consensus thresholds, on synthetic reads
     whose ground truth is known;
  2. blastomere routing re-run on the committed example matrix (a 1/40 cell
     subsample of the published callset), asserting it reproduces the published
     defining alleles and every published rate;
  3. the published full-run counts, read back from tables/ so they cannot drift.

    python3 tests/test_pipeline.py
"""
import csv
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import config                                                       # noqa: E402
from tape.barcode import correct                                    # noqa: E402
from tape.denoise import denoise_barcode                            # noqa: E402
from tape.io import load_vocab, load_whitelist                      # noqa: E402
from tape.parser import pattern_of, walk_robust                     # noqa: E402
from tape.sc import (chain_to_str, collapse_umis, consensus_locus,  # noqa: E402
                     extract_bc_sc, parse_header)

FAILURES = []
BC = b"CGGGGAATTGTA"            # a real embryo-#3 integration barcode
SYMS = [b"ACT", b"GCC", b"AAG", b"CCC", b"ACG", b"GAA"]


def check(label, got, want):
    ok = got == want
    print(f"  {'PASS' if ok else 'FAIL'}  {label}")
    if not ok:
        print(f"        got      {got!r}")
        print(f"        expected {want!r}")
        FAILURES.append(label)


def approx(label, got, want, tol):
    ok = abs(got - want) <= tol
    print(f"  {'PASS' if ok else 'FAIL'}  {label}  ({got} vs {want} +/- {tol})")
    if not ok:
        FAILURES.append(label)


# --------------------------------------------------------------------------
# 1. per-read parse + the three consensus thresholds
# --------------------------------------------------------------------------
def sc_read(bc=BC, n_edited=6, syms=SYMS):
    """One single-cell amplicon read: starts AT the barcode, no GGCATG flank."""
    array = b""
    for i in range(config.N_SITES):
        junction = (syms[i] + b"GGA" + b"TGAT") if i < n_edited else b"TGAT"
        array += config.MONO + junction
    return (bc + config.BC_RIGHT_SPACER + config.BC_RIGHT + config.LEADIN
            + array + config.TERM_LONG + b"CCCCCCCCCC")


def genotype_of(read, wl, vocab):
    """The full per-read funnel, exactly as scripts/sc_consensus.py applies it."""
    bc, cstart = extract_bc_sc(read)
    if bc is None:
        return None
    corr = correct(bc, wl)
    if corr is None:
        return None
    state, err, resolved_all = walk_robust(read, cstart)
    if err or not resolved_all:
        return None
    if not all((x == "U") or (x[0] == "E" and x[1] in vocab) for x in state):
        return None
    return corr, pattern_of(state)


print("\n1. per-read parse and barcode extraction")
wl = load_whitelist(config.get_embryo("DTTz_3_S3").whitelist_tsv)
vocab = load_vocab(config.get_embryo("DTTz_3_S3").vocab_tsv)
check("11 whitelist barcodes, 13 vocabulary symbols", (len(wl), len(vocab)), (11, 13))

bc, cstart = extract_bc_sc(sc_read())
check("barcode is the 12 bp 5' of the constant G before CTCTGG", bc, BC)
check("array start is just past CTCTGG",
      sc_read()[cstart:cstart + len(config.LEADIN)], config.LEADIN)

TRUTH = tuple(s.decode() for s in SYMS)
check("fully edited array parses to its ground truth", genotype_of(sc_read(), wl, vocab),
      (BC, TRUTH))
check("partially edited array fills U downstream (ordered editing)",
      genotype_of(sc_read(n_edited=2), wl, vocab)[1],
      ("ACT", "GCC", "U", "U", "U", "U"))
check("an off-vocabulary insertion is rejected",
      genotype_of(sc_read(syms=[b"TTT"] + SYMS[1:]), wl, vocab), None)
check("a barcode absent from the whitelist is rejected",
      genotype_of(sc_read(bc=b"AAAAAAAAAAAA"), wl, vocab), None)

print("\n2. header parsing")
cell, umi, sample = parse_header(b"@exp1_PA-04C.LIG-P8-F01_RT-P5-F08,PCR-PA-04C,CGGAAAGG,S")
check("cell id is the transcriptome identifier", cell, b"exp1_PA-04C.LIG-P8-F01_RT-P5-F08")
check("UMI is 8 bp", (umi, len(umi)), (b"CGGAAAGG", config.UMI_LEN))

print("\n3. UMI collapse -- reads to molecules")
deep, shallow = TRUTH, ("ACT", "GCC", "U", "U", "U", "U")
check("a 1-read UMI is discarded (MIN_READS_PER_UMI=2)",
      collapse_umis({b"AAAAAAAA": {deep: 1}}), Counter())
check("a 2-read UMI becomes one molecule",
      collapse_umis({b"AAAAAAAA": {deep: 2}}), Counter({deep: 1}))
check("UMIs within one mismatch merge into the more abundant one",
      collapse_umis({b"AAAAAAAA": {deep: 5}, b"AAAAAAAT": {deep: 1}}), Counter({deep: 1}))
check("UMIs two mismatches apart stay separate",
      collapse_umis({b"AAAAAAAA": {deep: 5}, b"AAAAAATT": {deep: 2}}), Counter({deep: 2}))
check("a molecule takes the dominant genotype among its reads",
      collapse_umis({b"AAAAAAAA": {deep: 5, shallow: 1}}), Counter({deep: 1}))

print("\n4. de-noising on molecule counts")
check("a single-molecule genotype is dropped (MIN_READS=2 on molecules)",
      denoise_barcode(Counter({deep: 6, ("ACT", "GCC", "AAG", "CCC", "ACG", "GAT"): 1})),
      {deep: 7})   # 1-mismatch neighbour of the last symbol -> collapsed into it
check("a well-supported second genotype survives de-noising",
      set(denoise_barcode(Counter({deep: 6, ("GGC",) + TRUTH[1:]: 4}))),
      {deep, ("GGC",) + TRUTH[1:]})

print("\n5. the three consensus thresholds")
check("under-supported locus is left missing (<3 molecules)",
      consensus_locus(Counter({deep: 2})), None)
res = consensus_locus(Counter({deep: 3}))
check("3 molecules on one lineage are called", (chain_to_str(res[0]), res[1], res[2], res[3]),
      ("-".join(TRUTH), False, 1.0, 3))
res = consensus_locus(Counter({deep: 3, shallow: 5}))
check("prefix (carry-over) molecules count AS support, and the call is the deepest tip",
      (chain_to_str(res[0]), res[2], res[3]), ("-".join(TRUTH), 1.0, 8))
branch = ("ACT", "GCC", "GGC", "U", "U", "U")
res = consensus_locus(Counter({deep: 20, branch: 2}))
check("an off-path branch with >=2 molecules flags a doublet but does not truncate the call",
      (chain_to_str(res[0]), res[1]), ("-".join(TRUTH), True))
approx("...and it lowers dominance to 1 - branching/total", res[2], 20 / 22, 1e-9)
check("a locus whose dominant lineage falls below 0.9 dominance is left missing",
      consensus_locus(Counter({deep: 10, branch: 2})), None)

# --------------------------------------------------------------------------
# 2. blastomere routing on the committed example matrix
# --------------------------------------------------------------------------
print("\n6. blastomere routing, re-run on example_data/ (1/40 cell subsample)")
EXAMPLE = ROOT / "example_data" / "e3v8_subsample.cell_consensus.tsv.gz"
if not EXAMPLE.exists():
    print(f"  SKIP  {EXAMPLE.name} not present")
else:
    with tempfile.TemporaryDirectory() as td:
        cmd = [sys.executable, str(ROOT / "scripts" / "blastomere_route.py"), str(EXAMPLE),
               "--outdir", td, "--out", "sub",
               "--min-n-loci", "4", "--max-doublet-loci", "1", "--min-mean-dom", "0.95"]
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.returncode != 0:
            print("  FAIL  blastomere_route.py exited non-zero")
            print(p.stderr[-2000:])
            FAILURES.append("routing run")
        else:
            report = (Path(td) / "sub.routing_report.txt").read_text()

            def num(label):
                for line in report.splitlines():
                    if line.strip().startswith(label):
                        return line.split(":", 1)[1].strip()
                return None

            check("11 discriminating integrations", num("discriminating tapes").split()[0], "11")
            check("20 defining founder alleles", num("defining alleles"), "20")
            # the alleles themselves must be the published ones
            def load_defs(path):
                rows = [l.split("\t") for l in Path(path).read_text().splitlines()[1:] if l]
                return {(r[0], r[2], r[3]) for r in rows}
            check("the derived alleles are identical to the published set",
                  load_defs(Path(td) / "sub.defining_sites.tsv"),
                  load_defs(ROOT / "tables" / "e3v8.defining_sites.tsv"))
            check("99.0% of cells route", num("total routed").split()[1], "(99.0%)")
            check("1.04% are set aside", num("set aside total").split()[1], "(1.04%)")
            check("84.0% pass QC (>=4 loci, <=1 doublet, dominance >=0.95)",
                  num("total pass_qc=1").split()[-1], "(84.0%)")

            # --- founder-state filter, on the matrices routing just wrote -------
            fs = subprocess.run(
                [sys.executable, str(ROOT / "scripts" / "founder_state_filter.py"), td,
                 "--defs", str(Path(td) / "sub.defining_sites.tsv"), "--out", td],
                capture_output=True, text=True)
            if fs.returncode != 0:
                print("  FAIL  founder_state_filter.py exited non-zero")
                print(fs.stderr[-1500:])
                FAILURES.append("founder-state run")
            else:
                got = {}
                cur = None
                for line in fs.stdout.splitlines():
                    if line.startswith("--- B"):
                        cur = line.split()[1].rstrip(":")
                    elif line.startswith("      ") and cur:
                        b, ch = line.split()
                        got[(cur, b)] = ch
                pub = {}
                for row in csv.DictReader(
                        open(ROOT / "tables" / "e3v8.founder_chains.tsv"), delimiter="\t"):
                    pub[(row["blastomere"], row["integration"])] = row["founder_chain"]
                # blastomere A reproduces exactly from 1/40 of the cells
                check("founder chains for blastomere A match the published set",
                      {k: v for k, v in got.items() if k[0] == "B1"},
                      {k: v for k, v in pub.items() if k[0] == "B1"})
                # B has one integration (TGGGGAACATAT, 3.7% recovery) whose chain
                # over-extends on a subsample; the rest must match
                b2_got = {k: v for k, v in got.items() if k[0] == "B2"
                          and k[1] != "TGGGGAACATAT"}
                check("founder chains for blastomere B match, bar the 3.7%-recovery tape",
                      b2_got, {k: v for k, v in pub.items() if k[0] == "B2"
                               and k[1] != "TGGGGAACATAT"})
                rates = [float(l.split("(")[1].rstrip("%)\n"))
                         for l in fs.stdout.splitlines() if "contradict the founder" in l]
                for got_r, want_r, lbl in zip(rates, [1.01, 1.39, 1.50, 0.91],
                                              ["A >=7", "A 4-6", "B >=7", "B 4-6"]):
                    approx(f"founder-contradiction rate, {lbl} tapes", got_r, want_r, 0.5)

# --------------------------------------------------------------------------
# 3. the published full-run counts
# --------------------------------------------------------------------------
print("\n7. published full-run values (tables/e3v8.routing_report.txt)")
pub = (ROOT / "tables" / "e3v8.routing_report.txt").read_text()
for label, want in [("cells", "1,559,419"),
                    ("defining alleles", "20"),
                    ("routed B1", "889,849"),
                    ("routed B2", "653,335"),
                    ("total routed", "1,543,184 (99.0%)"),
                    ("masked (1 bleed->NA)", "58,711"),
                    ("set aside total", "16,235 (1.04%)")]:
    line = next((l for l in pub.splitlines() if l.strip().startswith(label)), "")
    check(f"{label} = {want}", line.split(":", 1)[1].strip() if ":" in line else None, want)
check("pass_qc total = 1,296,799 of 1,543,184 (84.0%)",
      "1,296,799 / 1,543,184 (84.0%)" in pub, True)

recov = (ROOT / "tables" / "e3v8.sc_recovery_by_integration.csv").read_text()
chains = (ROOT / "tables" / "e3v8.founder_chains.tsv").read_text().splitlines()
check("20 founder chains are shipped, over 11 + 9 integrations", len(chains) - 1, 20)

check("recovery denominator is the 1,584,848 corrected nuclei",
      "1,584,848" in recov, True)
check("1,559,419 cells (98.4%) carry >=1 consensus genotype",
      "1,559,419 (98.4%)" in recov, True)

print()
if FAILURES:
    print(f"{len(FAILURES)} check(s) FAILED: {FAILURES}")
    sys.exit(1)
print("all single-cell pipeline checks passed")
