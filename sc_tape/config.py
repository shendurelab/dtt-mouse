#!/usr/bin/env python3
"""
Central configuration for the single-cell circtape (DNA Typewriter) pipeline.

SINGLE SOURCE OF TRUTH for every constant and threshold in this package;
nothing else hard-codes one. This is the single-cell counterpart of the bulk
release's config.py: the read body downstream of the integration barcode is
identical in the two assays, so the read-structure and de-noising constants are
shared verbatim, and only barcode extraction and everything after the UMI
differ. The bulk-only knobs (reference derivation, depth model, rarefaction,
tree rendering) are not reproduced here -- see the bulk package.

Read architecture of one single-cell amplicon read (the R2 / amplicon mate):

    5' [TAPE-BC 12bp] G CTCTGG ATGAT
       a1 site1 a2 site2 a3 site3 a4 site4 a5 site5 a6 site6  TERM  3'-const...

  A single-cell read STARTS at the integration barcode -- there is no 5' GGCATG
  flank as in the bulk amplicon -- so the barcode is the 12 bp immediately 5' of
  the constant G that precedes CTCTGG (tape/sc.py::extract_bc_sc). Everything
  from CTCTGG onward is parsed by the shared tape/parser.py.

  - CTCTGG          : 6-bp right flank of the barcode.
  - BC_RIGHT_SPACER : an invariant 'G' between the 12-bp barcode and CTCTGG,
                      stripped during extraction.
  - ATGAT           : constant lead-in immediately before monomer 1.
  - GGTGAGCACG      : monomer body anchor (a1..a6).
  - GGTGAGCCAC      : terminal monomer -- marks the end of the 6-site array.
  - junction (per site): unedited = 'TGAT'; edited = <NNN><GGA>'TGAT'.

  Editing is strictly ordered/unidirectional: site k is edited only if sites
  1..k-1 are all edited. A read with an edit after an unedited site is an
  ordering violation and is dropped.

The integration-barcode whitelist and insertion vocabulary are NOT derived here.
They are the embryo-#3 reference sets from the bulk analysis, shipped verbatim in
refs/ (see README).
"""
from __future__ import annotations
import os
from dataclasses import dataclass
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT      = Path(__file__).resolve().parent
REF_DIR   = ROOT / "refs"        # embryo-#3 whitelist + vocabulary, from bulk
# Generated outputs (gitignored). Single-cell intermediates are large (a merged
# fastq is ~3 GB); point TAPE_DATA_DIR at a local work dir, not a synced folder.
DATA_DIR  = Path(os.environ.get("TAPE_DATA_DIR", str(ROOT / "data"))).expanduser()
TABLES_DIR = ROOT / "tables"     # small, version-controlled figure-source tables

# Where the raw fastqs live. They are large and NOT checked in.
RAW_DIR   = Path(os.environ.get("TAPE_RAW_DIR", str(DATA_DIR / "raw"))).expanduser()

# ---------------------------------------------------------------------------
# Read-structure constants (bytes; the fastq is scanned as bytes)
# ---------------------------------------------------------------------------
BC_LEFT   = b"GGCATG"          # left flank of the TAPE-BC (BULK reads only)
BC_RIGHT  = b"CTCTGG"          # right flank of the TAPE-BC (6 bp)
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

# Correct an observed barcode to the whitelist iff a single nearest member sits
# within this Hamming distance (the second-nearest strictly farther). Embryo #2's
# 17 barcodes are >=4 mismatches from every embryo-#3 barcode, so embryo-#2
# molecules co-processed in the same sci-RNA-seq3 experiment can never be
# corrected onto an embryo-#3 integration; they are dropped here.
BC_CORRECT_MAX = 2

# ---------------------------------------------------------------------------
# De-noising thresholds -- shared with bulk, applied to MOLECULE counts here
# ---------------------------------------------------------------------------
COLLAPSE_RATIO = 2.0           # error-collapse: fold a genotype into a 1-symbol neighbor
                               #   only if the neighbor is >= this x more abundant
MIN_READS      = 2             # drop genotypes with < this much support after collapse
                               #   (molecules, in the single-cell pipeline)
CROSS_BC_DEPTH = 3             # bulk-only: TAPE-BC-swap recombinant depth. The
                               #   cross-barcode filter is a property of the bulk
                               #   POPULATION and is NOT applied per cell; the constant
                               #   is kept only so tape/denoise.py imports cleanly.

# ---------------------------------------------------------------------------
# Edit-chain model (shared with the bulk pipeline)
# ---------------------------------------------------------------------------
# A lineage is an ordered chain of edit symbols (5'->3', up to the first
# unedited site). Because editing is monotonic, a shallower chain is a PREFIX of
# a deeper one; in a single cell such a molecule is 3'-dropout or carry-over of
# the cell's deeper genotype, not a competing allele.
DROPOUT_FOLD_RATIO = 4.0       # bulk-only tree folding; the single-cell consensus
                               #   folds every prefix (one genotype per cell x locus)
BRANCH_MIN_FRAC  = 0.25        # used by editchain.is_branch/consensus, which the
BRANCH_MIN_READS = 2           #   single-cell caller does NOT use -- it calls
                               #   editchain.dominant_lineage, where an off-path
                               #   branch counts from >=2 molecules with no fraction
                               #   test. Kept for import compatibility.

# ---------------------------------------------------------------------------
# Single-cell settings
# ---------------------------------------------------------------------------
# SC read header, written by the first stage of the sci-RNA-seq3 pipeline:
#   "@<cell_id>,PCR-<well>,<UMI>,<sample>"
#   e.g. exp1_PA-04C.LIG-P8-F01_RT-P5-F08,PCR-PA-04C,CGGAAAGG,<sample>
# <cell_id> is the transcriptome cell identifier, so the character matrix joins
# to the cell x gene metadata with no translation.
UMI_LEN         = 8
UMI_MERGE_DIST  = 1            # merge UMIs within this Hamming distance (PCR/seq error)

# Canonical cell metadata (cell_id -> celltype / trajectory / UMAP), used to
# annotate the character matrix and as the recovery denominator. v8 =
# cell_metadata.v8.txt, the cross-talk-corrected re-demux (1,584,848 nuclei).
CELL_METADATA   = os.environ.get("TAPE_CELL_METADATA", "cell_metadata.v8.txt")

# --- single-cell consensus calling thresholds (the published values) --------
# A (cell x integration) locus is CALLED only when the dominant lineage is both
# well-supported and clearly dominant; otherwise it is left MISSING, on the
# principle that a wrong call harms the tree while a missing one does not.
MIN_READS_PER_UMI = 2          # (#1) a UMI must have >= this many reads to count as a molecule
SC_MIN_MOLECULES  = 3          # (#2) the called (dominant) lineage must have >= this many molecules
SC_DOMINANCE_MIN  = 0.9        # (#3) dominant-lineage molecules / all molecules must be >= this
                               #      (prefixes = 3' dropout/carry-over count AS support; only
                               #       branching lineages -- a different edit at the same site --
                               #       count against)

# fastp read-merging parameters (paired single-cell runs only). fastp emits
# merged reads in R1 orientation, so they are reverse-complemented to the
# amplicon orientation. Merge is best-effort; unmerged pairs fall back to R2.
FASTP_BIN                     = os.environ.get("FASTP_BIN", "fastp")
FASTP_OVERLAP_LEN_REQUIRE     = 15
FASTP_OVERLAP_DIFF_LIMIT      = 8
FASTP_OVERLAP_DIFF_PCT_LIMIT  = 30

# ---------------------------------------------------------------------------
# Reference sets: the embryo-#3 whitelist + vocabulary, derived on the BULK data
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Embryo:
    tag: str
    label: str
    note: str = ""

    @property
    def whitelist_tsv(self) -> Path:
        return REF_DIR / f"{self.tag}.tapebc_whitelist.tsv"

    @property
    def vocab_tsv(self) -> Path:
        return REF_DIR / f"{self.tag}.nnn_vocab.tsv"

    def out(self, suffix: str) -> Path:
        return DATA_DIR / f"{self.tag}.{suffix}"


EMBRYOS = {
    # The only embryo with single-cell data. 11 integrations, site-6 edited in
    # ~53% of resolvable bulk reads; 13 admitted insertion symbols.
    "DTTz_3_S3": Embryo("DTTz_3_S3", "embryo 3",
                        note="11 TAPE barcodes; the single-cell analysis focus"),
}


def get_embryo(tag: str) -> Embryo:
    if tag not in EMBRYOS:
        raise KeyError(f"unknown embryo tag {tag!r}; known: {sorted(EMBRYOS)}")
    return EMBRYOS[tag]


# ---------------------------------------------------------------------------
# Single-cell sequencing runs
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class SCSample:
    name: str                  # run name
    layout: str                # "single" (amplicon read only) | "paired" (R1+R2)
    embryo_tag: str            # which reference set to use
    r2: str                    # amplicon read fastq (single-end read, or R2 of a pair)
    r1: str = ""               # mate fastq (paired only)
    sample_field: str = ""     # optional: keep only reads whose header sample contains this

    @property
    def r1_path(self) -> Path:
        return RAW_DIR / self.r1

    @property
    def r2_path(self) -> Path:
        return RAW_DIR / self.r2

    @property
    def embryo(self) -> Embryo:
        return get_embryo(self.embryo_tag)

    def out(self, suffix: str) -> Path:
        return DATA_DIR / f"{self.name}.{suffix}"


# The two runs behind the published callset, both embryo #3. They are called
# TOGETHER in one pass (tag e3v8), pooling reads by (cell, integration, UMI) so a
# molecule sequenced in both runs is counted once:
#
#   seq8      single-end, 183-bp amplicon read; covers cells from BOTH
#             sci-RNA-seq3 experiments (25 PCR plates: exp1 PA-PW, exp2 PX/PY).
#   seq8_long paired, 2x250; covers cells from the FIRST experiment only
#             (exp1; verified by a full header census), re-sequenced at greater
#             read length. Needs 00_merge_pairs.py first.
#
# These are the CORRECTED re-demultiplexing (received 2026-08-21) that fixes
# sci-RNA-seq3 cell-barcode cross-talk. It REPLACES the earlier seq4/seq5 demux
# rather than pooling with it -- the same underlying reads, so pooling would
# retain the misassignments. Cell ids differ from the earlier demux, so the
# matching cell_metadata.v8.txt is required for any celltype/tree join.
SC_SAMPLES = {
    "seq8":      SCSample("seq8", "single", "DTTz_3_S3",
                          r2="Tape_merged.seq8.R2.fastq.gz"),
    "seq8_long": SCSample("seq8_long", "paired", "DTTz_3_S3",
                          r2="Tape_merged.seq8_long.R2.fastq.gz",
                          r1="Tape_merged.seq8_long.R1.fastq.gz"),
}


def get_sc_sample(name: str) -> SCSample:
    if name not in SC_SAMPLES:
        raise KeyError(f"unknown SC sample {name!r}; known: {sorted(SC_SAMPLES)}")
    return SC_SAMPLES[name]
