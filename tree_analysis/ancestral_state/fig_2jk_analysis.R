#!/usr/bin/env Rscript
# single_fig_3_panels_de_analysis.R
# Parsimony edit-count reconstruction underlying Figure 2's panels J (editing
# rate over time) and K (internal-node support over time) -- mouse_sprint's
# fig_3 panels D/E. Runs the reconstruction once and writes the two binned,
# plot-ready tables fig2_ancestral_panel_D_editing_rate.csv /
# fig2_ancestral_panel_E_support_over_time.csv that making_figures/Figure-2.R
# eventually reads (staged under those fig2_ancestral_-prefixed names for now,
# alongside the existing panel_D_editing_rate.csv / panel_E_support_over_time.csv,
# so the two can be diffed before replacing the ones Figure-2.R actually reads).
#
# Needs (see this directory's tests/ + README for status):
#   results/cell_by_tape.rds -- long-format (cell_id, barcode, Site1..Site6)
#   tip tape states, reshaped from support_data/e3v8.B{1,2}_tape_consensus.tsv.gz.
#   Not yet built here; see 01_cell_by_tape.R's port status.
#
# Run from tree_analysis/:
#   Rscript ancestral_state/single_fig_3_panels_de_analysis.R

suppressMessages({ library(ape) })
source("tree_analysis/ancestral_state/01_parsimony_lib.R")
source("tree_analysis/ancestral_state/02_support_lib.R")            # SITE_COLS, reconstruct_support(), support_over_time_table()
source("tree_analysis/ancestral_state/03_tape_tip_state_mapping.R") # load_recoded_tip_states()

#---- 1. inputs ---------------------------------------------------------------
TREE_NWK <- "tree_building/results/3-dated-tree/merged_minB2h_lineage_constrained.nwk"
stopifnot("tree not found" = file.exists(TREE_NWK))

RESULTS_DIR <- "tree_analysis/ancestral_state/results"

CELL_BY_TAPE_RDS <- file.path(RESULTS_DIR, "cell_by_tape.rds")
stopifnot("cell_by_tape.rds not found -- build it from support_data/e3v8.B{1,2}_tape_consensus.tsv.gz first (reshape not yet ported, see this directory's README/tests)" =
            file.exists(CELL_BY_TAPE_RDS))

FIG_DATA_DIR <- "figures_data/"
D_CSV <- file.path(FIG_DATA_DIR, "fig2_ancestral_panel_D_editing_rate.csv")
E_CSV <- file.path(FIG_DATA_DIR, "fig2_ancestral_panel_E_support_over_time.csv")

cat("tree  :", TREE_NWK, "\n")

#---- 2. tip tape states (whole embryo) --------------------------------------
tip_states_long <- load_recoded_tip_states(CELL_BY_TAPE_RDS)
barcodes <- sort(unique(tip_states_long$barcode))
cat(length(barcodes), "tapes;", length(unique(tip_states_long$cell_id)), "cells in tape table\n")

#---- 3. tree: binarise (multi2di) for the reconstruction's binary contract --
tree <- ape::multi2di(ape::read.tree(TREE_NWK), random = FALSE)
stopifnot("tree must be binary after multi2di" = ape::is.binary(tree))
Ntip <- length(tree$tip.label); N <- Ntip + tree$Nnode
root <- .get_root(tree)
node_day <- ape::node.depth.edgelength(tree)   # merged root = day 0 exactly
cat(sprintf("tree: %d tips, %d internal; root=E%.2f, tips=E%.2f\n",
            Ntip, tree$Nnode, node_day[root], max(node_day)))

#---- 4. blastomere membership: which root child each node descends from ----
# A = larger root-child clade (Blastomere A = B1), B = smaller.
root_children <- tree$edge[tree$edge[, 1] == root, 2]
stopifnot("expected a bifurcating zygote root" = length(root_children) == 2L)
csize <- vapply(root_children, function(k)
  length(phangorn::Descendants(tree, k, type = "tips")[[1]]), integer(1))
root_children <- root_children[order(-csize)]
A_node <- root_children[1]; B_node <- root_children[2]
memb <- rep(NA_integer_, N); memb[A_node] <- 1L; memb[B_node] <- 2L
edge_pre <- reorder(tree, "cladewise")$edge     # parents before children
for (i in seq_len(nrow(edge_pre))) {
  p <- edge_pre[i, 1]; c <- edge_pre[i, 2]
  if (is.na(memb[c])) memb[c] <- memb[p]
}
blastomere <- c("A", "B")[memb]                 # NA for the root only
cat(sprintf("blastomere split: A(=B1) %d tips  B(=B2) %d tips\n",
            sum(memb[seq_len(Ntip)] == 1L), sum(memb[seq_len(Ntip)] == 2L)))

#---- 5. the (heavy, cached) parsimony reconstruction over the whole tree ----
recon <- get_or_compute_full_tree_reconstruction(
  cache_path  = file.path(RESULTS_DIR, "recon_merged.rds"),
  input_paths = c(TREE_NWK, CELL_BY_TAPE_RDS),
  compute_fn  = function() reconstruct_support(tree, tip_states_long, barcodes))
support_df <- recon$support_df                  # node, <tape cols>, support
stopifnot(nrow(support_df) == N)
support_df$embryonic_day <- node_day
support_df$blastomere    <- blastomere
saveRDS(support_df, file.path(RESULTS_DIR, "support_df.rds"))
write.csv(support_df, file.path(RESULTS_DIR, "support_df.csv"), row.names = FALSE)
cat("wrote", file.path(RESULTS_DIR, "support_df.{rds,csv}"), "\n")

#---- 6. Panel D data: editing rate (edits/day) vs date, per blastomere -----
# Rate per 1-day window = (edits on branches whose midpoint date falls in the
# window, summed over all tapes) / (branch-time in the window, days). Branch
# = one tree edge; edits on it = support of its child node; blastomere =
# child's. The merged root is the zygote at day 0, so the two E0->E1.5
# blastomere stem branches are real edges here, included automatically.
BIN_EDGES <- seq(0, 14, by = 1.0)
child  <- tree$edge[, 2]; parent <- tree$edge[, 1]
edits      <- support_df$support[child]
branch_day <- node_day[child] - node_day[parent]
mid        <- (node_day[child] + node_day[parent]) / 2
blast      <- support_df$blastomere[child]

keep <- branch_day > 1e-9 & !is.na(blast)       # drop zero-length multi2di edges
agg  <- data.frame(blastomere = blast[keep], edits = edits[keep],
                   days = branch_day[keep], mid = mid[keep])
agg$bin <- cut(agg$mid, breaks = BIN_EDGES, include.lowest = TRUE)

lo <- head(BIN_EDGES, -1); hi <- tail(BIN_EDGES, -1)
grid <- expand.grid(blastomere = c("A", "B"), bin = levels(agg$bin), stringsAsFactors = FALSE)
sums <- aggregate(cbind(edits, days) ~ blastomere + bin, data = agg, sum)
rate_df <- merge(grid, sums, all.x = TRUE)
rate_df$edits[is.na(rate_df$edits)] <- 0
rate_df$days[is.na(rate_df$days)]   <- 0
rate_df$day_mid <- ((lo + hi) / 2)[match(rate_df$bin, levels(agg$bin))]
rate_df$rate    <- ifelse(rate_df$days > 0, rate_df$edits / rate_df$days, NA_real_)
rate_df <- rate_df[order(rate_df$blastomere, rate_df$day_mid), ]
rate_df$blastomere <- factor(rate_df$blastomere, levels = c("A", "B"))

write.csv(rate_df, D_CSV, row.names = FALSE)
cat("wrote", D_CSV, "\n")
print(rate_df[, c("blastomere", "day_mid", "edits", "days", "rate")], row.names = FALSE, digits = 3)

#---- 7. Panel E data: % of internal nodes with support >= threshold, per day bin
THRESHOLDS <- c("support >= 1 edit" = 1L, "support >= 2 edits" = 2L)
tab <- support_over_time_table(tree, support_df, support_df$support, BIN_EDGES, THRESHOLDS)
cat(sprintf("internal nodes counted (real branches, A+B): %d\n", sum(tab$n_internal_nodes)))

write.csv(tab, E_CSV, row.names = FALSE)
cat("wrote", E_CSV, "\n")
print(tab, row.names = FALSE, digits = 3)

cat("done.\n")

