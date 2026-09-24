# lib/qc.R
# QC-pass cell adapter for the e3v8 B1/B2 tape-consensus matrices, each
# carrying n_loci/n_doublet_loci/mean_dominance/pass_qc columns. There is no
# single canonical consensus file (two per-blastomere files), so callers must
# always pass tsv_path explicitly.
#
# QC-pass rule: pass_qc == 1 AND n_loci >= V8_MIN_LOCI (env, default 7 -- the
# high-quality backbone threshold). This rule does NOT include the
# founder-consistency (divergent-cell) filter -- that's applied separately by
# 1_build_nj_backbone/00_filter_highqual_consensus.R, which derives the
# ge7_founderok.tsv.gz files the backbone is actually built from (those runs
# use QC_PASS_MODE=all, trusting the file as already filtered, rather than
# this function).
#
# Requires lib/paths.R to already be sourced.
#
# Usage:
#   source(file.path(REPO_DIR, "lib/paths.R"))
#   source(file.path(REPO_DIR, "lib/qc.R"))
#   pass_ids <- qc_pass_cell_ids(REPO_DIR, "v8", tsv_path = DTT_TSV)

qc_pass_cell_ids <- function(repo, version, tsv_path = NULL) {
  if (version != "v8") {
    stop(sprintf("qc_pass_cell_ids: no QC rule defined for DATA_VERSION '%s' (only v8 is supported)", version))
  }
  if (is.null(tsv_path)) {
    stop("qc_pass_cell_ids: v8 has no single consensus file -- pass tsv_path ",
         "(the specific B1/B2 file being processed) explicitly.")
  }
  min_loci <- as.integer(Sys.getenv("V8_MIN_LOCI", "7"))
  df      <- read.delim(tsv_path, colClasses = "character", check.names = FALSE)
  n_loci  <- as.integer(df$n_loci)
  pass_qc <- as.integer(df$pass_qc)
  df$cell_id[pass_qc == 1L & n_loci >= min_loci]
}
