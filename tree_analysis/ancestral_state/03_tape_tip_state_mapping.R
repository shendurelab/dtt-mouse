# 03_tape_tip_state_mapping.R
# Bridges the real tape data to the "0"/"?" sentinel convention
# 01_parsimony_lib.R / 02_support_lib.R expect (issue #21).
#
# Reuses src/4-time-tree/01_cell_by_tape.R's output as-is: a long table
# (cell_id, barcode, Site1..Site6) already following
# src/1-dist-mat/parse_tape_consensus.R's token conventions --
#   "ETY"  = unedited placeholder (site was recorded but never cut)
#   "None" = whole tape not recovered for that cell (per-site, since an
#            observed tape always has all 6 sites; see 01_cell_by_tape.R
#            step 6)
#   anything else = an edit outcome, kept verbatim
# -- and recodes only the two sentinels: "ETY" -> "0", "None" -> "?".
# No re-parsing of the raw consensus TSV happens here.
#
# Run (example, from repo root):
#   source("src/lib/paths.R")
#   source("src/8-support/02_support_lib.R")   # defines SITE_COLS, sourced first everywhere in this repo
#   source("src/8-support/03_tape_tip_state_mapping.R")
#   tip_states_long <- load_recoded_tip_states(results_path(".", "4-time-tree", "cell_by_tape.rds"))
#   barcodes <- sort(unique(tip_states_long$barcode))
#
# SITE_COLS (the 6 sequential site columns) is defined once in
# 02_support_lib.R and reused here rather than re-literalized.

# Recode one long cell-by-tape table's site columns in place: "ETY" -> "0"
# (unedited/reference sentinel), "None" -> "?" (missing sentinel), anything
# else passed through verbatim as the edit token.
recode_tape_tokens <- function(long_df, site_cols = SITE_COLS) {
  for (col in site_cols) {
    v <- long_df[[col]]
    v[v == "ETY"]  <- "0"
    v[v == "None"] <- "?"
    long_df[[col]] <- v
  }
  long_df
}

# Load 01_cell_by_tape.R's RDS and apply the sentinel recoding.
#   cell_by_tape_rds : path to the cell_by_tape.rds produced by
#                      src/4-time-tree/01_cell_by_tape.R for this DATA_VERSION
load_recoded_tip_states <- function(cell_by_tape_rds, site_cols = SITE_COLS) {
  long_df <- readRDS(cell_by_tape_rds)
  recode_tape_tokens(long_df, site_cols)
}
