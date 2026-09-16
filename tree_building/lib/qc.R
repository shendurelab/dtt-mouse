# lib/qc.R
# QC-pass cell adapter: isolates the v1/v2/v3 schema difference behind one
# version-agnostic function so stage 1 and stage 4 scripts don't need to know
# which QC columns exist in which data version.
#
#   v1: e3_tape_consensus.tsv carries pass_qc (0/1) directly.
#   v2: e3_tape_consensus.tsv.gz dropped pass_qc/doublet_score/lowqual_score;
#       QC now lives in a separate e3.sc_qc.tsv.gz (n_loci, n_doublet_loci,
#       mean_dominance, total_molecules) with no explicit pass/fail flag.
#       QC-pass rule (chosen 2026-07-07): n_loci >= 5 & doublet_frac < 0.3,
#       where doublet_frac = n_doublet_loci / n_loci. On the full v2 delivery
#       this passes 182,574 / 208,439 cells (87.6%).
#   v3: e3.sc_qc.tsv.gz has the SAME columns as v2's, so the same rule is
#       reused as-is. Note this is NOT a free pass on the same result: on the
#       full v3 delivery it passes only 138,784 / 205,072 cells (67.7%) -- a
#       ~20-point drop from v2 despite an unchanged formula, worth treating as
#       a flag to double check against the v3 delivery notes rather than
#       assuming the rule still fits.
#   v4: two pre-vetted subsets (B1/B2 tape consensus matrices from the
#       blastomere-routing stage), each carrying n_loci/n_doublet_loci/
#       mean_dominance/pass_qc as PASSTHROUGH columns from the v3 population
#       QC file (pass_qc is always "1" here -- routing's own inclusion
#       criteria, not a QC verdict). QC-pass rule (chosen 2026-07-17,
#       corrects an earlier oversight where v4 runs used QC_PASS_MODE=all,
#       i.e. no filter at all): n_loci >= 6 & mean_dominance >= 0.9 -- a
#       DIFFERENT rule than v2/v3's (higher n_loci floor, a dominance floor
#       instead of a doublet-fraction ceiling), read directly from
#       tsv_path (there is no single v4 consensus file to default to; see
#       lib/paths.R).
#   v5: e3v5 B1/B2 tape consensus matrices (standalone seq4 call), same
#       passthrough columns as v4 plus a pass_qc flag that already encodes
#       v5's shipped gate (n_loci>=4 & n_doublet_loci<=1 & mean_dominance>=0.95).
#       QC-pass rule: pass_qc == 1 AND n_loci >= 5 -- keeps v5's doublet/
#       dominance gates (via pass_qc) but tightens the loci floor from 4 to 5.
#       Selects 367,047 (B1) / 249,154 (B2) cells. Read directly from tsv_path
#       (no single v5 consensus file; see lib/paths.R).
#   v6: e3v5v6 B1/B2 tape consensus matrices (v6 reprocessing of the seq4 e3v5
#       call), same passthrough columns + pass_qc flag as v5. QC-pass rule:
#       pass_qc == 1 AND n_loci >= V6_MIN_LOCI (env, default 7) -- same shape as
#       v5's rule with a configurable floor (7 is the high-quality backbone
#       threshold; see docs/HANDOVER_ge7_backbone_rebuild.md). IMPORTANT: this
#       QC rule does NOT include the founder-consistency (divergent-cell) filter
#       -- the ge7_founderok.tsv.gz files this repo ships already have that
#       filter applied upstream. Read directly from tsv_path (no single
#       v6 consensus file; see lib/paths.R).
#
# Requires lib/paths.R to already be sourced (uses consensus_tsv_path/data_path).
#
# Usage:
#   source(file.path(REPO_DIR, "lib/paths.R"))
#   source(file.path(REPO_DIR, "lib/qc.R"))
#   pass_ids <- qc_pass_cell_ids(REPO_DIR, DATA_VERSION)                # v1-v3
#   pass_ids <- qc_pass_cell_ids(REPO_DIR, "v4", tsv_path = DTT_TSV)    # v4/v5

qc_pass_cell_ids <- function(repo, version, tsv_path = NULL) {
  if (version == "v1") {
    df <- read.delim(consensus_tsv_path(repo), colClasses = "character", check.names = FALSE)
    df$cell_id[df$pass_qc == "1"]
  } else if (version %in% c("v2", "v3")) {
    qc <- read.delim(data_path(repo, "e3.sc_qc.tsv.gz"), colClasses = "character")
    n_loci         <- as.integer(qc$n_loci)
    n_doublet_loci <- as.integer(qc$n_doublet_loci)
    doublet_frac   <- n_doublet_loci / n_loci
    qc$cell[n_loci >= 5 & doublet_frac < 0.3]
  } else if (version == "v4") {
    if (is.null(tsv_path)) {
      stop("qc_pass_cell_ids: v4 has no single consensus file -- pass tsv_path ",
           "(the specific B1/B2 file being processed) explicitly.")
    }
    df <- read.delim(tsv_path, colClasses = "character", check.names = FALSE)
    n_loci         <- as.integer(df$n_loci)
    mean_dominance <- as.numeric(df$mean_dominance)
    df$cell_id[n_loci >= 6 & mean_dominance >= 0.9]
  } else if (version == "v5") {
    if (is.null(tsv_path)) {
      stop("qc_pass_cell_ids: v5 has no single consensus file -- pass tsv_path ",
           "(the specific B1/B2 file being processed) explicitly.")
    }
    # loci floor is env-configurable (V5_MIN_LOCI, default 5). Raising it shrinks
    # the tree quadratically in memory/time: >=7 (B1 225,490 / B2 135,232 cells)
    # fits the 1.9 TB box via the in-memory shim, avoiding the 3.8 TB box + quota.
    min_loci <- as.integer(Sys.getenv("V5_MIN_LOCI", "5"))
    df      <- read.delim(tsv_path, colClasses = "character", check.names = FALSE)
    n_loci  <- as.integer(df$n_loci)
    pass_qc <- as.integer(df$pass_qc)
    df$cell_id[pass_qc == 1L & n_loci >= min_loci]
  } else if (version == "v6") {
    if (is.null(tsv_path)) {
      stop("qc_pass_cell_ids: v6 has no single consensus file -- pass tsv_path ",
           "(the specific B1/B2 file being processed) explicitly.")
    }
    # Same shape as v5's rule, floor env-configurable via V6_MIN_LOCI (default 7,
    # the high-quality backbone threshold). NOTE: founder-consistency is NOT
    # applied here -- the ge7_founderok.tsv.gz files this repo ships already
    # have that filter applied upstream.
    min_loci <- as.integer(Sys.getenv("V6_MIN_LOCI", "7"))
    df      <- read.delim(tsv_path, colClasses = "character", check.names = FALSE)
    n_loci  <- as.integer(df$n_loci)
    pass_qc <- as.integer(df$pass_qc)
    df$cell_id[pass_qc == 1L & n_loci >= min_loci]
  } else {
    stop(sprintf("qc_pass_cell_ids: no QC rule defined for DATA_VERSION '%s'", version))
  }
}
