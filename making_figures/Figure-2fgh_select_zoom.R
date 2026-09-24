#!/usr/bin/env Rscript
# 02_select_nested_zoom.R  (Figure 3B)
# Pick two NESTED, cell-type-diverse zoom-in clades from the current merged time
# tree: zoom1 ~ TARGET1 tips, and zoom2 ~ TARGET2 tips chosen from WITHIN zoom1
# (so zoom2 is a true zoom into zoom1). Diversity = number of distinct
# major_trajectory cell types among a clade's tips; among size-matched
# candidates the most diverse wins (ties -> closest to target size).
#
# Efficiency: descendant tip counts AND per-node type-diversity are both
# computed in ONE postorder pass over the 360k-tip tree -- diversity via a
# 25-bit type bitmask OR'd up the tree (popcount = distinct types), so no
# per-node descendant-set query is needed.
#
# Writes, under making_figures/Figure-2fgh_output/zoom/ (where Figure-2fgh.py
# expects them):
#   zoom1.nwk, zoom1_day_offset.txt          (clade + its root's embryonic day)
#   zoom2.nwk, zoom2_day_offset.txt
#   zoom2_tips.txt                           (zoom2 tip labels -> highlight box on zoom1)
#
# Run from the repo root:
#   Rscript making_figures/Figure-2fgh_select_zoom.R
# Ported from mouse_sprint's src/3-plot-tree/02_select_nested_zoom.R -- same
# tree Figure-2fgh.py's panel F and single_fig_3_panels_de_analysis.R use.

suppressMessages({ library(ape) })

TARGET1 <- 1000L; RANGE1 <- c(800L, 1300L)
TARGET2 <- 100L;  RANGE2 <- c(70L, 140L)
TYPE_COL <- "major_trajectory"

# ---- resolve tree + metadata ----------------------------------------------
TREE_NWK <- "tree_building/results/3-dated-tree/merged_minB2h_lineage_constrained.nwk"
META     <- "support_data/cell_metadata.v8.txt.gz"
stopifnot(file.exists(TREE_NWK), file.exists(META))
OUT_DIR <- "making_figures/Figure-2fgh_output/zoom"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("tree:", TREE_NWK, "\n"); tree <- read.tree(TREE_NWK)
Ntip <- length(tree$tip.label); Nnode <- tree$Nnode; N <- Ntip + Nnode
cat(sprintf("  %d tips, %d internal nodes\n", Ntip, Nnode))

# ---- tip -> type bit (>=25 categories fit in 32-bit int) -----------------
meta <- read.delim(META, colClasses = "character", check.names = FALSE)
type_of <- setNames(meta[[TYPE_COL]], meta$cell_id)
tip_type <- type_of[tree$tip.label]
cats <- sort(unique(tip_type[!is.na(tip_type)]))
stopifnot("too many cell types for a 31-bit mask" = length(cats) <= 31)
tip_bit <- match(tip_type, cats) - 1L            # 0-based bit index; NA = unknown (no bit)
cat(sprintf("  %d major_trajectory categories; %d/%d tips typed\n",
            length(cats), sum(!is.na(tip_bit)), Ntip))

# ---- one postorder pass: descendant tip counts + type-mask ---------------
count <- integer(N); count[seq_len(Ntip)] <- 1L
mask  <- integer(N)
has_bit <- !is.na(tip_bit)
mask[which(has_bit)] <- bitwShiftL(1L, tip_bit[has_bit])
edge_po <- reorder(tree, "postorder")$edge         # children before parents
for (i in seq_len(nrow(edge_po))) {
  p <- edge_po[i, 1]; c <- edge_po[i, 2]
  count[p] <- count[p] + count[c]
  mask[p]  <- bitwOr(mask[p], mask[c])
}
popcount <- function(m) { k <- 0L; while (any(m > 0)) { k <- k + (m %% 2L); m <- bitwShiftR(m, 1L) }; k }
ndistinct <- popcount(mask)                        # vectorised over all nodes

depth <- node.depth.edgelength(tree)               # embryonic day (merged root = day 0)

# ---- pick the most type-diverse clade in a size window -------------------
pick <- function(candidate_nodes, target) {
  # rank by distinct types (desc), then by closeness to target (asc)
  o <- order(-ndistinct[candidate_nodes], abs(count[candidate_nodes] - target))
  candidate_nodes[o[1]]
}

internal <- (Ntip + 1):N
z1_cand <- internal[count[internal] >= RANGE1[1] & count[internal] <= RANGE1[2]]
stopifnot("no clade near TARGET1" = length(z1_cand) > 0)
zoom1 <- pick(z1_cand, TARGET1)
cat(sprintf("zoom1: node %d  %d tips  %d cell types  root E%.2f\n",
            zoom1, count[zoom1], ndistinct[zoom1], depth[zoom1]))

# zoom2 candidates: internal nodes strictly BELOW zoom1 (so zoom2 nests in zoom1)
below <- unlist(phangorn::Descendants(tree, zoom1, type = "all"))
below <- below[below > Ntip]
z2_cand <- below[count[below] >= RANGE2[1] & count[below] <= RANGE2[2]]
stopifnot("no sub-clade near TARGET2 inside zoom1" = length(z2_cand) > 0)
zoom2 <- pick(z2_cand, TARGET2)
cat(sprintf("zoom2: node %d  %d tips  %d cell types  root E%.2f  (nested in zoom1)\n",
            zoom2, count[zoom2], ndistinct[zoom2], depth[zoom2]))

# ---- write clade artifacts ------------------------------------------------
write_clade <- function(node, stem) {
  sub <- extract.clade(tree, node)
  write.tree(sub, file.path(OUT_DIR, paste0(stem, ".nwk")))
  writeLines(sprintf("%.6f", depth[node]), file.path(OUT_DIR, paste0(stem, "_day_offset.txt")))
  sub
}
sub1 <- write_clade(zoom1, "zoom1")
sub2 <- write_clade(zoom2, "zoom2")
# Tip lists drive the nested-zoom highlight boxes: zoom1's tips -> box on the
# full tree (panel a), zoom2's tips -> box on zoom1 (panel b).
writeLines(sub1$tip.label, file.path(OUT_DIR, "zoom1_tips.txt"))
writeLines(sub2$tip.label, file.path(OUT_DIR, "zoom2_tips.txt"))
cat("wrote zoom1/zoom2 clades to", OUT_DIR, "\n")
