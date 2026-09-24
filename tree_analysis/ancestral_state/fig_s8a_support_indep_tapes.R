#!/usr/bin/env Rscript
# fig_s8a_support_indep_tapes.R
# Fig. S8A: node support counted by DISTINCT TAPES edited on the incoming
# branch, not total edited sites (Fig. 2K's metric) -- a branch with 3 edited
# sites in one tape + 1 in another scores 4 under Fig. 2K, 2 here. 
#
# Reuses fig_2jk_analysis.R's ALREADY-CACHED whole-tree reconstruction
# (results/support_df.rds) -- no new tree reconstruction here, just a cheap
# derived per-node column (count of tapes with >0 edits, capped at 1 each)
# rebinned with 02_support_lib.R's support_over_time_table(), same as Fig. 2K.
#
# Run from the repo root (after fig_2jk_analysis.R):
#   Rscript tree_analysis/ancestral_state/fig_s8a_support_indep_tapes.R

suppressMessages({ library(ape) })
source("tree_analysis/ancestral_state/02_support_lib.R")   # support_over_time_table()

TREE_NWK <- "tree_building/results/3-dated-tree/merged_minB2h_lineage_constrained.nwk"
RESULTS_DIR <- "tree_analysis/ancestral_state/results"
RECON <- file.path(RESULTS_DIR, "support_df.rds")
stopifnot("tree not found" = file.exists(TREE_NWK),
          "run fig_2jk_analysis.R first" = file.exists(RECON))

OUT_CSV <- "figures_data/fig_s8_ancestral_panel_A_support_over_time_indep_tapes.csv"

# ---- 1. reload the SAME binarised tree + cache fig_2jk_analysis.R used -----
tree <- ape::multi2di(ape::read.tree(TREE_NWK), random = FALSE)
support_df <- readRDS(RECON)

# ---- 2. per-node independent-tape count -- cap each tape's edit count at 1
#         before summing across tapes (support_df$support sums the raw counts)
tape_cols <- setdiff(colnames(support_df), c("node", "support", "embryonic_day", "blastomere"))
cat(length(tape_cols), "tapes:", paste(tape_cols, collapse = ", "), "\n")
n_indep_tapes <- rowSums(support_df[, tape_cols] > 0)

# ---- 3. same aggregation core as fig_2jk_analysis.R's panel K, thresholding
#         n_indep_tapes instead, with an extra support>=3 threshold ---------
BIN_EDGES <- seq(0, 14, by = 1.0)     # 1-day windows from the zygote (E0), same as panel K
THRESHOLDS <- c("support >= 1" = 1L, "support >= 2" = 2L, "support >= 3" = 3L)
tab_indep <- support_over_time_table(tree, support_df, n_indep_tapes, BIN_EDGES, THRESHOLDS)
cat(sprintf("internal nodes counted (real branches, A+B): %d\n", sum(tab_indep$n_internal_nodes)))

# ---- 4. same table for the ORIGINAL total-edits metric (Fig. 2K's), at the
#         same 3 thresholds -- lets the plot overlay both as a comparison.
#         Recomputed here rather than reusing fig2_ancestral_panel_E_...csv,
#         which only has thresholds 1-2 (Fig. 2K's own scope).
tab_total <- support_over_time_table(tree, support_df, support_df$support, BIN_EDGES, THRESHOLDS)

tab <- rbind(cbind(metric = "indep_tapes", tab_indep), cbind(metric = "total_edits", tab_total))
write.csv(tab, OUT_CSV, row.names = FALSE)
cat("wrote", OUT_CSV, "\n")
print(tab, row.names = FALSE, digits = 3)
