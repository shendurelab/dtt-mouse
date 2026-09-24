#!/usr/bin/env Rscript
# fig_s7d_blastomere_support.R
# Fig. S7D: distribution of per-internal-node edit-count support (number of
# tape-sites edited on the branch leading into a node), for Blastomere A (B1)
# and Blastomere B (B2) separately.
#
# Run from the repo root:
#   Rscript tree_analysis/ancestral_state/fig_s7d_blastomere_support.R

suppressMessages({ library(ape) })

source("tree_analysis/ancestral_state/01_parsimony_lib.R")
source("tree_analysis/ancestral_state/02_support_lib.R")   # reconstruct_support()
source("tree_analysis/ancestral_state/03_tape_tip_state_mapping.R")  # load_recoded_tip_states()

# ---- 1. Inputs (constants at top) ----------------------------------------
# The two per-blastomere BAT-corrected divergence trees (dtt-mouse's
# tree_nonneg.nwk equivalent -- see tree_building/3_date_tree/README.md).
TREES <- list(
  B1 = "tree_building/results/2-rooted-nj/divergence_nonneg_B1.nwk",
  B2 = "tree_building/results/2-rooted-nj/divergence_nonneg_B2.nwk"
)
CELL_BY_TAPE_RDS <- "tree_analysis/ancestral_state/results/cell_by_tape.rds"
stopifnot(file.exists(TREES$B1), file.exists(TREES$B2), file.exists(CELL_BY_TAPE_RDS))

OUT_DIR <- "tree_analysis/ancestral_state/results/blastomere_support"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
FIG_S7D_CSV <- "figures_data/fig_s7_ancestral_panel_A_support_histogram.csv"

# ---- 2. Tape tip states (shared by both trees) ---------------------------
tip_states_long <- load_recoded_tip_states(CELL_BY_TAPE_RDS)
barcodes <- sort(unique(tip_states_long$barcode))
cat(length(barcodes), "tapes:", paste(barcodes, collapse = ", "), "\n")

# ---- 3. Per tree: reconstruct support, keep the per-node vector ----------
per_tree <- function(label, tree_path) {
  cat("\n==", label, "==", tree_path, "\n")
  tree <- ape::read.tree(tree_path)
  # multi2di is a no-op on these (already strictly binary) but reconstruct_*()
  # requires binary; resolve deterministically so node ids are reproducible.
  tree_bin <- ape::multi2di(tree, random = FALSE)
  stopifnot(ape::is.binary(tree_bin))
  cat("  tips:", length(tree_bin$tip.label), " internal nodes:", tree_bin$Nnode,
      " max root-dist:", round(max(ape::node.depth.edgelength(tree_bin)), 2), "days\n")

  result <- reconstruct_support(tree_bin, tip_states_long, barcodes)
  support <- result$support_df$support

  Ntip <- length(tree_bin$tip.label)
  internal_idx <- (Ntip + 1):(Ntip + tree_bin$Nnode)
  data.frame(blastomere = label, node = internal_idx, support = support[internal_idx])
}

pernode_df <- do.call(rbind, Map(per_tree, names(TREES), unlist(TREES)))
rownames(pernode_df) <- NULL

# per-node support -> CSV (drives Fig. S7D's histogram; one row per internal
# node). Written both as a local, inspectable artifact and as the staged
# figures_data copy Figure-S7.R's panel D will eventually read.
write.csv(pernode_df, file.path(OUT_DIR, "blastomere_support_pernode.csv"), row.names = FALSE)
write.csv(pernode_df, FIG_S7D_CSV, row.names = FALSE)
cat("wrote", FIG_S7D_CSV, "\n")
