#!/usr/bin/env python3
"""
Central configuration for the bulk TAPE (DNA Typewriter) analysis pipeline.

RELEASE BUILD -- bulk only. This is the bulk half of the pipeline described in
the Methods ("Primary analysis of bulk tape"). The single-cell front-end
(tape/sc.py, tape/merge.py and their settings) is released separately; the
core library below -- the array parser, the insertion vocabulary, the
de-noising funnel and the edit-chain model -- is shared verbatim between the
two assays. Bulk-derived references (refs/) are what the single-cell pipeline
consumes for the same embryo.

This is the SINGLE SOURCE OF TRUTH for every constant, threshold, path, and
per-embryo setting. Nothing else in the pipeline should hard-code a magic
number -- import it from here. If you re-run on a new library, add an Embryo()
entry below and point the driver at its tag.

Read architecture of one R1 read (318 bp single-end), verified from the data:

    5' ...GGCATG [TAPE-BC 12bp] G CTCTGG ATGAT
       a1 site1 a2 site2 a3 site3 a4 site4 a5 site5 a6 site6  TERM  3'-const...

  - GGCATG / CTCTGG : 6-bp left/right flanks bracketing the barcode. A constant 'G'
                      (BC_RIGHT_SPACER) sits between the 12-bp barcode and CTCTGG and
                      is stripped during extraction.
  - ATGAT           : constant lead-in immediately before monomer 1.
  - GGTGAGCACG      : monomer body anchor (a1..a6).
  - GGTGAGCCAC      : terminal monomer -- marks the end of the 6-site array.
  - junction (per site): unedited = 'TGAT'; edited = <NNN><GGA>'TGAT'
                         (NNN = 3-bp insertion barcode; occasional 4-bp e.g. GATG).

  Editing is strictly ordered/unidirectional: site k is edited only if sites
  1..k-1 are all edited. A read with an edit after an unedited site is an
  ordering violation and is dropped.
"""
from __future__ import annotations
import os
from dataclasses import dataclass, field
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT      = Path(__file__).resolve().parent
REF_DIR   = ROOT / "refs"        # checked-in whitelists + NNN vocabularies
# Generated outputs (gitignored). Outputs land here by default; set
# TAPE_DATA_DIR to keep large intermediates outside the repository.
DATA_DIR  = Path(os.environ.get("TAPE_DATA_DIR", str(ROOT / "data"))).expanduser()
FIG_DIR   = DATA_DIR / "figures"
# Small, version-controlled figure-source tables (checked into git, unlike data/).
TABLES_DIR = ROOT / "tables"

# Where the raw fastqs live. They are large and are NOT distributed with the code.
# Default is data/raw; set the TAPE_RAW_DIR env var to point elsewhere. See the
# README ("Raw data") for how to restage them.
RAW_DIR   = Path(os.environ.get("TAPE_RAW_DIR", str(DATA_DIR / "raw"))).expanduser()

# ---------------------------------------------------------------------------
# Read-structure constants (bytes; we scan the raw fastq as bytes)
# ---------------------------------------------------------------------------
BC_LEFT   = b"GGCATG"          # left flank of the TAPE-BC (6 bp)
BC_RIGHT  = b"CTCTGG"          # right flank of the TAPE-BC (6 bp, symmetric with BC_LEFT)
BC_RIGHT_SPACER = b"G"         # constant base(s) between the 12-bp barcode and BC_RIGHT
                               #   (an invariant 'G' across all barcodes in both embryos);
                               #   stripped during extraction to recover the true 12-bp BC.
LEADIN    = b"ATGAT"           # constant lead-in before monomer 1
MONO      = b"GGTGAGCACG"      # monomer body anchor
TERM      = b"GGTGAGCCAC"      # terminal monomer (end of array)
ANCHOR    = b"GGTGAG"          # short anchor shared by body + terminal monomers
TERM_LONG = b"GGTGAGCCACAAGGCAAT"  # 18-bp terminal landmark for fuzzy 3' location

N_SITES   = 6                  # monomer sites per array (recorder design)
BC_LEN    = 12                 # TAPE-BC length (design spec; 2 internal positions are
                               #   a fixed 'AA' spacer, so ~10 positions are informative)

# Valid spacing (bp) between consecutive anchors: 14 (unedited) up through a
# 3-4 bp edit; the [13, 24] window tolerates small wobble while rejecting a
# missing/garbled monomer (gap >= ~26).
GAP_MIN   = 13
GAP_MAX   = 24

# ---------------------------------------------------------------------------
# Reference-derivation thresholds (PASS A: whitelist + NNN vocab)
# ---------------------------------------------------------------------------
SAMPLE_N       = 2_000_000     # reads sampled to derive whitelist + vocab
BC_MIN_FREQ    = 0.002         # min fraction of BC-extractable reads for a whitelist seed
BC_MIN_SEP     = 3             # min Hamming separation between accepted whitelist members
BC_CORRECT_MAX = 2             # correct a read's BC to whitelist iff unique nearest within this

# Insertion vocabulary rule. An insertion is admitted iff, at some single
# (TAPE-BC x site position) combination, it is observed at least VOCAB_PEAK_MIN
# times. Real edits are clonal and concentrate at a specific integration+site,
# so genuine insertions peak far above sequencing/PCR error, leaving a clean gap
# for the threshold to sit in. Measured on the checked-in refs/: embryo 3, 13
# admitted insertions with a lowest peak of 72,667 against a highest rejected
# peak of 5,446 (13.3x); embryo 2, 11 admitted, 14,589 vs 2,350 (6.2x); embryo
# 6, 1 admitted, 28,813 vs 2,997 (9.6x). tests/test_pipeline.py re-derives these
# from refs/ so this comment cannot drift out of date again. Derived on the
# BULK data only and reused for single cell (see reference.py / DATA.md).
VOCAB_PEAK_MIN = 10_000

# ---------------------------------------------------------------------------
# De-noising thresholds (patterns_full -> clean patterns)
# ---------------------------------------------------------------------------
COLLAPSE_RATIO = 2.0           # error-collapse: fold pattern into a 1-symbol neighbor
                               #   only if the neighbor is >= this x more abundant
MIN_READS      = 2             # drop patterns with < this many reads after collapse
CROSS_BC_DEPTH = 3             # a pattern of edit-depth >= this whose max-read "home"
                               #   barcode differs is a TAPE-BC-swap recombinant -> drop
MIN_BC_READS   = 2000          # only analyze barcodes with >= this many total reads

# ---------------------------------------------------------------------------
# Edit-chain model (shared with the single-cell pipeline)
# ---------------------------------------------------------------------------
# A lineage is an ordered chain of edit symbols (5'->3', up to the first
# unedited site). Because editing is monotonic, a shallower chain is a PREFIX of
# a deeper one; a read/molecule ending at an interior node is either a genuine
# early stop OR 3'-dropout/carry-over of a deeper lineage.
#
# DROPOUT_FOLD_RATIO: fold an interior tip into its dominant deeper child
#   (treat it as dropout/truncation, NOT a real stop) only when the deeper
#   subtree outweighs it by at least this factor. Larger = more conservative
#   (keeps more shallow lineages as genuine stops). Set to 0 to disable folding.
#   The single-cell consensus uses full folding (one genotype per cell x locus);
#   bulk keeps many lineages per barcode, so folding is conservative here.
DROPOUT_FOLD_RATIO = 4.0
# A node is a real BRANCH (two competing edit symbols at the same next site --
# impossible from a single lineage; a doublet/mixture or a true sub-lineage
# split) when a 2nd child clears both of these:
BRANCH_MIN_FRAC  = 0.25        # >= this fraction of the node's through-weight
BRANCH_MIN_READS = 2           # and >= this absolute weight

# ---------------------------------------------------------------------------
# Depth model (reads -> genomic equivalents / cells)
# ---------------------------------------------------------------------------
LAMBDA_MIN_READS = 200         # ignore patterns below this when estimating lambda
LAMBDA_KDE_BW    = 0.3         # gaussian_kde bandwidth (log10 space) for the lambda mode

# ---------------------------------------------------------------------------
# Rarefaction (lineage-diversity accumulation vs cells sampled)
# ---------------------------------------------------------------------------
RAREFACTION_EXTRAP = 5.0       # extrapolate the accumulation curve to this x the observed cells
RAREFACTION_INTERP_PTS = 40    # grid points on the interpolation (<= observed) segment
RAREFACTION_EXTRAP_PTS = 60    # grid points on the Chao1 extrapolation (> observed) segment

# ---------------------------------------------------------------------------
# Tree rendering
# ---------------------------------------------------------------------------
TREE_MIN_READS = 100           # per-integration trees drawn from patterns with >= this many reads
TREE_TOPK      = 28            # max patterns drawn per integration (pattern-trie view)
CAPLEAF        = 45            # cap drawn leaves per clone in the cell-rake view

# Insertion -> color map, assigned in order of genome-wide insertion abundance.
PALETTE = ["#2E6FB7", "#D1495B", "#3DA35D", "#E0A100", "#7E57C2", "#16A0A0",
           "#E0731A", "#C2185B", "#5D8C3A", "#8D6E63", "#0C6B74", "#B0306A"]
UNEDITED_COLOR = "#E2E2E2"     # grey = unedited site
OTHER_COLOR    = "#9AA7AD"     # insertion outside the top-12 palette


# ---------------------------------------------------------------------------
# Per-embryo settings
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Embryo:
    tag: str                   # file/output prefix, e.g. "DTTz_3_S3"
    label: str                 # human name, e.g. "embryo 3"
    fastq: str                 # raw fastq filename (expected under RAW_DIR)
    note: str = ""

    @property
    def fastq_path(self) -> Path:
        return RAW_DIR / self.fastq

    @property
    def whitelist_tsv(self) -> Path:
        return REF_DIR / f"{self.tag}.tapebc_whitelist.tsv"

    @property
    def vocab_tsv(self) -> Path:
        return REF_DIR / f"{self.tag}.nnn_vocab.tsv"

    def out(self, suffix: str) -> Path:
        """Path for a generated output, e.g. embryo.out('patterns_full.pkl')."""
        return DATA_DIR / f"{self.tag}.{suffix}"


EMBRYOS = {
    # Barely edited.
    "DTTz_1_S1": Embryo("DTTz_1_S1", "embryo 1", "DTTz_1_S1_R1_001.fastq.gz",
                        note="scarcely edited"),
    # Cold recorder: site-6 edited in only ~3.2% of resolvable reads; 17 barcodes.
    "DTTz_2_S2": Embryo("DTTz_2_S2", "embryo 2", "DTTz_2_S2_R1_001.fastq.gz",
                        note="cold recorder (site6 ~3.2%), 17 TAPE barcodes"),
    # Hot recorder: site-6 edited in ~53%; 11 barcodes. The lineage-analysis focus.
    "DTTz_3_S3": Embryo("DTTz_3_S3", "embryo 3", "DTTz_3_S3_R1_001.fastq.gz",
                        note="hot recorder (site6 ~53%), 11 TAPE barcodes; analysis focus"),
    # Barely edited.
    "DTTz_4_S4": Embryo("DTTz_4_S4", "embryo 4", "DTTz_4_S4_R1_001.fastq.gz",
                        note="scarcely edited"),
    # Embryo 6: sequenced paired-end (unlike the single-end embryos 1-4). The
    # amplicon (barcode->6 sites->terminal) is carried on R2 in reverse-complement
    # orientation; R1 (70bp) is a partial 3' read with no barcode/terminal and is
    # unused. DTTz_6_S6_R1_001.fastq.gz is the reverse-complement of the original
    # embryo6_S2_R2_001.fastq.gz, so it is a drop-in single-end amplicon matching
    # embryos 1-4. See DATA.md.
    "DTTz_6_S6": Embryo("DTTz_6_S6", "embryo 6", "DTTz_6_S6_R1_001.fastq.gz",
                        note="paired-end source; RC(R2) used as single-end amplicon"),
}


def get_embryo(tag: str) -> Embryo:
    if tag not in EMBRYOS:
        raise KeyError(f"unknown embryo tag {tag!r}; known: {sorted(EMBRYOS)}")
    return EMBRYOS[tag]

