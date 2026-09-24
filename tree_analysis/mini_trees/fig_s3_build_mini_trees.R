# fig_s3_build_mini_trees.R
# Fig. S3: build one midpoint-rooted NJ tree per TAPE integration (11 for
# embryo 3) from the bulk lineage-genotype table, over that integration's
# distinct observed genotypes -- no single-cell data used. Analysis only;
# plotting (tree + tip-aligned tape heatmap) is making_figures/Figure-S3.R.
#
# Run from the repo root:
#   Rscript tree_analysis/mini_trees/fig_s3_build_mini_trees.R

source("tree_analysis/mini_trees/mini_tree_lib.R")   # build_mini_tree

GENOTYPES_CSV <- "bulk_tape/tables/DTTz_3_S3.bulk_lineage_genotypes.csv"
TREES_DIR <- "tree_analysis/mini_trees/results/trees"
dir.create(TREES_DIR, showWarnings = FALSE, recursive = TRUE)

genotypes <- read.csv(GENOTYPES_CSV, stringsAsFactors = FALSE)
by_tapebc <- split(genotypes, genotypes$tapebc)
by_tapebc <- by_tapebc[order(-lengths(lapply(by_tapebc, `[[`, "genotype_id")))]

for (tapebc in names(by_tapebc)) {
  tree <- build_mini_tree(by_tapebc[[tapebc]])
  ape::write.tree(tree, file.path(TREES_DIR, paste0(tapebc, ".nwk")))
  cat(sprintf("[fig_s3] %s: %d genotypes -> %s\n",
              tapebc, nrow(by_tapebc[[tapebc]]), file.path(TREES_DIR, paste0(tapebc, ".nwk"))))
}

cat(sprintf("[fig_s3] built %d mini trees -> %s\n", length(by_tapebc), TREES_DIR))
