# test_mini_tree_lib.R
# Eyeball-checkable test: run the SAME production functions used for the full
# 11-tape build (parse_genotype, genotype_distance_matrix, build_mini_tree) on
# just 5 real genotypes from one real tapebc, and print every intermediate so
# the distance/tree logic can be hand-verified against the raw data.
#
# Run from the repo root: Rscript tree_analysis/mini_trees/test_mini_tree_lib.R

suppressMessages({library(ape); library(phangorn)})
source("tree_analysis/mini_trees/mini_tree_lib.R")

GENOTYPES_CSV <- "bulk_tape/tables/DTTz_3_S3.bulk_lineage_genotypes.csv"
TAPEBC <- "ATATCAAATTGA"
N_TEST <- 5

all_genotypes <- read.csv(GENOTYPES_CSV, stringsAsFactors = FALSE)
five <- all_genotypes[all_genotypes$tapebc == TAPEBC, ]
five <- five[order(-five$n_cells), ][seq_len(N_TEST), ]

# ---- 1. the 5 genotypes' edit status ---------------------------------------
cat("==== 1. Edit status of the", N_TEST, "test genotypes (tapebc", TAPEBC, ") ====\n")
site_mat <- t(vapply(five$genotype, parse_genotype, character(6)))
rownames(site_mat) <- five$genotype_id
colnames(site_mat) <- paste0("site", 1:6)
print(cbind(as.data.frame(site_mat),
            edit_depth = five$edit_depth, n_cells = five$n_cells))

# ---- 2. the 5x5 distance matrix ---------------------------------------------
cat("\n==== 2. Pairwise distance matrix (genotype_distance_matrix) ====\n")
d <- genotype_distance_matrix(five)
print(as.matrix(d))

# ---- 3. the midpoint-rooted NJ tree ------------------------------------------
cat("\n==== 3. Midpoint-rooted NJ tree (build_mini_tree) ====\n")
tree <- build_mini_tree(five)
cat("tips:", paste(tree$tip.label, collapse = ", "), "\n")
cat("rooted:", ape::is.rooted(tree), " | binary:", ape::is.binary(tree), "\n")
cat("Newick:", ape::write.tree(tree), "\n\n")
cat("edges (parent -> child, branch length):\n")
edge_df <- data.frame(
  parent = tree$edge[, 1], child = tree$edge[, 2], length = round(tree$edge.length, 3)
)
edge_df$child_label <- ifelse(edge_df$child <= length(tree$tip.label),
                               tree$tip.label[edge_df$child], paste0("internal_", edge_df$child))
print(edge_df)
