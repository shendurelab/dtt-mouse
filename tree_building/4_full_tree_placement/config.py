#!/usr/bin/env python3
"""Central configuration for the DTT-ancestor placement pipeline.

Paths and algorithm parameters shared across the pipeline's steps live here so
those modules stay path-free and importable. Every value can be overridden from
the environment, so the same code runs on a laptop or a chunked parallel job.
A few run-specific values that only one step needs (RATE, DATED_TREE_OUT, QC_OUT,
TOPK) are instead read directly by that step (clockrate.py, finalize_dated.py,
bestmatch.py) -- see run_side.sh for the full list of env vars a run accepts.

Data contracts (what the algorithm assumes about its inputs)
------------------------------------------------------------
CONSENSUS_TSV (one per side, gzipped TSV with a header):
    col 0                      cell_id (unique)
    then INTEGRATION columns   each a 6-site edit chain "a|b|c|d|e|f"
                               (site value "" / "NA" / "U" / "-" = uninformative;
                                a whole cell value "" / "NA" = integration absent)
    then QC columns            n_loci, n_doublet_loci, mean_dominance, pass_qc
                               (names listed in QC_COLS, excluded from genotype)

BACKBONE_NWK (one per side): the DATED per-side tree from 3_date_tree (branch
    lengths in days), whose tip labels are a subset of CONSENSUS_TSV's cell_ids.
    These tips are the anchors the tree is grown onto; they are never moved.
"""
import os

# ---- which lineage side we are processing (B1 / B2 are independent trees) ----
SIDE = os.environ.get("SIDE", "B1")

# ---- input paths (run_side.sh always sets these explicitly; empty default so
# importing config for another value, e.g. QC_COLS, never crashes -- consumers that
# actually need these fail with a clear file-not-found at the point of use) ----
CONSENSUS_TSV = os.environ.get("CONSENSUS_TSV", "")
BACKBONE_NWK = os.environ.get("BACKBONE_NWK", "")

# ---- blastomere-founder edit prefix, per (integration, side) -----------------
# Supp. Table 1: ordered prefix of sites fixed (>=90% consensus) within a blastomere
# at its 2-cell-stage founding. Excluded from the "shared edited prefix" placement-QC
# metric (dtt_lengths.py) -- shared by construction by every cell of a side, so
# uninformative about relatedness within the side.
FOUNDER_TABLE = os.environ.get(
    "FOUNDER_TABLE",
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                "..", "processed_data", "e3v5v6.supp_table1_founder_genotypes.csv"))

# ---- output location ---------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.environ.get("OUTDIR", os.path.join(HERE, f"phase_{SIDE}"))

# ---- algorithm parameters ----------------------------------------------------
N_SITES = int(os.environ.get("N_SITES", 6))       # edit sites per integration chain
QC_COLS = {"n_loci", "n_doublet_loci", "mean_dominance", "pass_qc"}

# Include all cells with >=1 informative locus (ALLCELLS=1) or only pass_qc==1 (=0).
ALLCELLS = int(os.environ.get("ALLCELLS", "0"))

# Rarity index: drop any edit token carried by more than CAP*N cells from the
# inverted index. Very common edits are shared by unrelated lineages (homoplasy /
# early shared edits) and add noise, so they never contribute a match. 0.2 = 20%.
CAP_FRAC = float(os.environ.get("CAP", 0.2))

# Chunking knobs for the (embarrassingly parallel) best-match phase. Each worker
# scores queries [START:END] and writes bestmatch_<OUTTAG>.npz; a shared cell
# ordering (cellids/isbb/nloci) is written once by the OUTTAG=="0" worker.
START = int(os.environ.get("START", 0))
END = int(os.environ.get("END", -1))            # -1 => to the end
OUTTAG = os.environ.get("OUTTAG", "0")

# Token id packing: (integration*N_SITES + site)*STRIDE + allele_id. STRIDE just
# needs to exceed the max alleles-per-site; 8192 is comfortably large.
TOK_STRIDE = 8192
