#!/usr/bin/env Rscript
# fig_s7c_node_time_sensitivity.R
# Fig. S7C: sensitivity of node dates to the assumed B1/B2 cell-count-ceiling
# split (0.5/0.5 default vs. asymmetric 0.63/0.37 and 0.58/0.42 re-datings).
#
# Run from the repo root:
#   Rscript tree_analysis/dating_qc/fig_s7c_node_time_sensitivity.R

suppressPackageStartupMessages({ library(ape); library(castor); library(dplyr); library(tidyr) })

TREE_FILES <- c(
  main              = "tree_building/results/3-dated-tree/merged_minB2h_lineage_constrained.nwk",
  "sensitivity_63-37" = "tree_building/results/3-dated-tree-sensitivity/merged_minB2h_lineage_constrained_sidefrac63-37.nwk",
  "sensitivity_58-42" = "tree_building/results/3-dated-tree-sensitivity/merged_minB2h_lineage_constrained_sidefrac58-42.nwk"
)
stopifnot(file.exists(TREE_FILES))
OUT_CSV <- "figures_data/fig_s7_dating_node_time_shifts.csv"

trees <- lapply(TREE_FILES, ape::read.tree)

# The trees have the exact same nodes, but different node dates. For each node
# in the sensitivity trees, quantify the difference in node time (root->node
# depth) vs. the main tree.
root_depths <- lapply(trees, castor::get_all_distances_to_root)

# ---- blastomere labels: larger root-child clade = A, smaller = B ----------
Ntip <- length(trees$main$tip.label)
root_id <- setdiff(trees$main$edge[, 1], trees$main$edge[, 2])
children <- trees$main$edge[trees$main$edge[, 1] == root_id, 2]
stopifnot(length(children) == 2)

sub1 <- castor::get_subtree_at_node(trees$main, children[1] - Ntip)
sub2 <- castor::get_subtree_at_node(trees$main, children[2] - Ntip)
if (length(sub1$new2old_tip) >= length(sub2$new2old_tip)) {
  bigger <- sub1; smaller <- sub2
} else {
  bigger <- sub2; smaller <- sub1
}

blastomere <- rep(NA_character_, Ntip + trees$main$Nnode)
blastomere[bigger$new2old_tip]         <- "A"
blastomere[Ntip + bigger$new2old_node] <- "A"
blastomere[smaller$new2old_tip]        <- "B"
blastomere[Ntip + smaller$new2old_node] <- "B"

# ---- per-node date differences vs. the main (default 0.5/0.5) tree --------
diff_63_37 <- root_depths[["sensitivity_63-37"]] - root_depths$main
diff_58_42 <- root_depths[["sensitivity_58-42"]] - root_depths$main

diffs <- data.frame(original_depth = root_depths$main,
                    blastomere = blastomere,
                    diff_63_37 = diff_63_37,
                    diff_58_42 = diff_58_42)

diffs <- diffs[diffs$original_depth < 13.49, ]   # exclude tip depths
diffs <- diffs[diffs$original_depth > 1.49, ]    # exclude the fixed root depth (blastomere = NA there)

diffs_long <- diffs %>%
  pivot_longer(cols = starts_with("diff_"), names_to = "sensitivity", values_to = "depth_diff") %>%
  mutate(depth_diff_hours = depth_diff * 24)   # days -> hours

write.csv(diffs_long, OUT_CSV, row.names = FALSE)
cat("wrote", OUT_CSV, "\n")
cat(sprintf("rows: %d (blastomere A=%d, B=%d)\n", nrow(diffs_long),
            sum(diffs$blastomere == "A"), sum(diffs$blastomere == "B")))
