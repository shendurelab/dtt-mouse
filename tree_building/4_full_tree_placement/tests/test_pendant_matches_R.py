#!/usr/bin/env python3
"""Cross-language oracle: our Python matrix pendant == canonical dtt_distance.R.

For real B2 cells we form ordered pairs (C, T), run them through the PRODUCTION
pendants_from_code (C as query with parent T), and check d(C,T) = pendant + residual
equals dtt_distance(C, T) from 1_build_nj_backbone/dtt_distance.R (the exact metric
that built the NJ backbone). This locks the whole pendant path to the project's
canonical distance. Skips gracefully if Rscript or the data are unavailable.
"""
import os
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
REPO = os.path.dirname(PKG)
sys.path.insert(0, PKG)

CONS = os.path.join(REPO, "processed_data/e3v5v6.B2_tape_consensus.tsv.gz")
NROWS = 120                     # first NROWS data rows -> pairs among them


def _have(cmd):
    from shutil import which
    return which(cmd) is not None


def _r_pairwise(nrows):
    """Run dtt_distance.R over the first `nrows` consensus cells; return {(id1,id2): dist}."""
    out_tsv = tempfile.mktemp(suffix=".tsv")
    rscript = f'''
      suppressMessages({{ }})
      source("{REPO}/1_build_nj_backbone/dtt_distance.R")
      source("{REPO}/1_build_nj_backbone/parse_tape_consensus.R")
      cells <- parse_cells(gzfile("{CONS}"), rows = 1:{nrows})
      ids <- names(cells)
      con <- file("{out_tsv}", "w")
      for (i in seq_along(ids)) for (j in i:length(ids)) {{
        d <- dtt_distance(cells[[i]], cells[[j]])
        if (!is.na(d)) cat(ids[i], ids[j], sprintf("%.10f", d), sep="\\t", file=con);
        if (!is.na(d)) cat("\\n", file=con)
      }}
      close(con)
    '''
    subprocess.run(["Rscript", "-e", rscript], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    d = {}
    for line in open(out_tsv):
        a, b, v = line.rstrip("\n").split("\t")
        d[(a, b)] = float(v)
    os.unlink(out_tsv)
    return d


def _py_dCT(pairs, ids_order):
    """d(C,T) for each (C,T) pair via the PRODUCTION pendants_from_code."""
    from genotypes import load_and_encode
    from dtt_lengths import pendants_from_code
    cells = load_and_encode(CONS, set(), 6, allcells=1)
    pos = {cid: i for i, cid in enumerate(cells.ids)}
    usable = [(a, b) for (a, b) in pairs if a in pos and b in pos]
    # build a code array with two rows per pair (C then T); parent pairs C->T
    rows, parent = [], []
    for a, b in usable:
        rows.append(cells.code[pos[a]]); rows.append(cells.code[pos[b]])
        parent.append(len(rows) - 1); parent.append(-1)     # C's parent = the T row just added
    code = np.array(rows, dtype=np.int16)
    pend, resid = pendants_from_code(code, np.array(parent))
    out = {}
    for k, (a, b) in enumerate(usable):
        out[(a, b)] = float(pend[2 * k] + resid[2 * k])     # d(C,P)+d(T,P) = d(C,T)
    return out


def test_matches_R():
    if not _have("Rscript") or not os.path.exists(CONS):
        print("  SKIP test_matches_R (Rscript or B2 consensus unavailable)")
        return
    r = _r_pairwise(NROWS)
    assert len(r) > 50, f"too few non-NA R pairs ({len(r)})"
    py = _py_dCT(list(r.keys()), None)
    checked = maxdiff = 0
    for k, rv in r.items():
        if k in py:
            maxdiff = max(maxdiff, abs(py[k] - rv)); checked += 1
    assert checked > 50, f"too few comparable pairs ({checked})"
    assert maxdiff < 1e-6, f"max |python - R| = {maxdiff}"
    print(f"  ok  test_matches_R ({checked} pairs, max|Δ|={maxdiff:.2e})")


def run():
    test_matches_R()
    print("test_pendant_matches_R: DONE")


if __name__ == "__main__":
    run()
