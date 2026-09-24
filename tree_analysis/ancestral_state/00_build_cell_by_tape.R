#!/usr/bin/env Rscript
# 00_build_cell_by_tape.R
# Reshapes the raw wide TAPE-consensus tables into the long, per-cell-per-tape
# table single_fig_3_panels_de_analysis.R needs (cell_id, barcode, Site1..Site6,
# tokens recoded). Input: support_data/e3v8.B{1,2}_tape_consensus.tsv.gz (each
# cell_id + 11 integration columns + n_loci/n_doublet_loci/mean_dominance/pass_qc).
#
# Run from the repo root:
#   Rscript tree_analysis/ancestral_state/00_build_cell_by_tape.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ape)
})

TREE_NWK <- "tree_building/results/4-full-tree/merged_full_placed.nwk"
B1_TSV   <- "support_data/e3v8.B1_tape_consensus.tsv.gz"
B2_TSV   <- "support_data/e3v8.B2_tape_consensus.tsv.gz"
OUT_RDS  <- "tree_analysis/ancestral_state/results/cell_by_tape.rds"
stopifnot(file.exists(TREE_NWK), file.exists(B1_TSV), file.exists(B2_TSV))

SITE_COLS <- paste0("Site", 1:6)

# wide (cell_id + 11 integration columns, "s1|..|s6" or NA) -> long
# (cell_id, barcode, Site1..Site6), recoded: whole-cell NA -> six "None",
# site "U" -> "ETY" (unedited), everything else kept verbatim (an edit).
build_cell_by_tape_long <- function(df, site_cols = SITE_COLS) {
  barcodes <- colnames(df)[2:12]   # column 1 is cell_id
  long <- df |>
    select(cell_id, all_of(barcodes)) |>
    pivot_longer(cols = all_of(barcodes), names_to = "barcode", values_to = "value") |>
    separate(value, into = site_cols, sep = "\\|", fill = "right", remove = TRUE)

  observed  <- rowSums(!is.na(long[site_cols])) > 0
  bad_split <- observed & rowSums(is.na(long[site_cols])) > 0
  if (any(bad_split)) {
    stop(sprintf("%d observed tapes did not split into 6 sites (bad delimiter?).", sum(bad_split)))
  }

  long |>
    mutate(across(all_of(site_cols), ~ case_when(
      is.na(.x)   ~ "None",
      .x == "U"   ~ "ETY",
      TRUE        ~ .x
    )))
}

# ---- 1. read + QC-filter each side separately, then combine ----------------
read_side <- function(tsv) {
  df <- read.delim(tsv, colClasses = "character", check.names = FALSE)
  df[df$pass_qc == "1", , drop = FALSE]
}
b1 <- read_side(B1_TSV)
b2 <- read_side(B2_TSV)
cat(sprintf("pass_qc cells: B1 %d, B2 %d\n", nrow(b1), nrow(b2)))

overlap <- intersect(b1$cell_id, b2$cell_id)
stopifnot("B1/B2 cell_id sets must be disjoint" = length(overlap) == 0)
df <- bind_rows(b1, b2)

# ---- 2. wide -> long + token recode -----------------------------------------
long <- build_cell_by_tape_long(df)
cat(sprintf("cell-by-tape rows: %d (%d cells x %d tapes)\n",
            nrow(long), length(unique(long$cell_id)), length(unique(long$barcode))))

# ---- 3. tip-coverage check, fail loudly here rather than deep inside the ---
#         parsimony reconstruction (which hard-errors on any missing tip) ----
tip_labels <- ape::read.tree(TREE_NWK)$tip.label
missing_tips <- setdiff(tip_labels, unique(long$cell_id))
cat(sprintf("tree tips: %d, missing from cell_by_tape: %d\n",
            length(tip_labels), length(missing_tips)))
stopifnot("tree tips missing from cell_by_tape -- QC filter or data mismatch" =
            length(missing_tips) == 0)

# ---- 4. write ----------------------------------------------------------------
dir.create(dirname(OUT_RDS), showWarnings = FALSE, recursive = TRUE)
saveRDS(long, OUT_RDS)
cat("wrote", OUT_RDS, "\n")
