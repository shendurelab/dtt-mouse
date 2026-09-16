# merge_dated_subtrees.R
# Merge two independently time-dated subtrees (B1/B2, each already LSD2-dated
# with the SAME root calibration day) into one tree with a synthetic shared
# MRCA at day 0.
#
suppressPackageStartupMessages(library(ape))

# merge_dated_subtrees(left, right, mrca_age, root_age)
#   left, right : ape::phylo, already rooted, DISJOINT tip label sets (real
#                 cell ids are unique across sides, but this is asserted, not
#                 assumed)
#   mrca_age    : age of the new shared root (e.g. 0 = day 0 / zygote)
#   root_age    : the (shared) age each subtree's own root was calibrated to
#                 (e.g. 1.5) -- must be strictly older (a positive stem)
merge_dated_subtrees <- function(left, right, mrca_age, root_age) {
  stopifnot(inherits(left, "phylo"), inherits(right, "phylo"))
  stopifnot(length(intersect(left$tip.label, right$tip.label)) == 0)
  stem <- root_age - mrca_age
  if (stem <= 0) {
    stop(sprintf("root_age (%.4f) must be strictly older than mrca_age (%.4f) -- got a non-positive stem",
                 root_age, mrca_age))
  }

  # write.tree()'s trailing ";" has to go so we can nest each subtree as a
  # clade inside a new outer set of parens with its own stem branch length.
  left_nwk  <- sub(";\\s*$", "", ape::write.tree(left))
  right_nwk <- sub(";\\s*$", "", ape::write.tree(right))
  merged <- ape::read.tree(text = sprintf("(%s:%.10f,%s:%.10f);", left_nwk, stem, right_nwk, stem))

  stopifnot(ape::is.rooted(merged))
  merged
}

# ---- CLI entry point (skipped when this file is source()'d for its function) ----
# Run from the repo root:
#   LEFT_NWK=results/3-dated-tree/B1/lsd2/nj99478/minB2h/time_tree.nwk \
#   RIGHT_NWK=results/3-dated-tree/B2/lsd2/nj99478/minB2h/time_tree.nwk \
#   OUT_NWK=results/3-dated-tree/merged/merged_time_tree.nwk \
#     Rscript 3_date_tree/merge_dated_subtrees.R
if (sys.nframe() == 0) {
  LEFT_NWK  <- Sys.getenv("LEFT_NWK")
  RIGHT_NWK <- Sys.getenv("RIGHT_NWK")
  OUT_NWK   <- Sys.getenv("OUT_NWK")
  MRCA_AGE  <- as.numeric(Sys.getenv("MRCA_AGE", "0"))
  ROOT_AGE  <- as.numeric(Sys.getenv("ROOT_AGE", "1.5"))
  if (!nzchar(LEFT_NWK) || !nzchar(RIGHT_NWK) || !nzchar(OUT_NWK)) {
    stop("LEFT_NWK, RIGHT_NWK, OUT_NWK must all be set (MRCA_AGE default 0, ROOT_AGE default 1.5)")
  }
  left  <- ape::read.tree(LEFT_NWK)
  right <- ape::read.tree(RIGHT_NWK)

  root_nchild <- function(t) sum(t$edge[, 1] == (ape::Ntip(t) + 1))
  left_depths  <- ape::node.depth.edgelength(left)[1:ape::Ntip(left)]
  right_depths <- ape::node.depth.edgelength(right)[1:ape::Ntip(right)]
  cat(sprintf("[merge] left:  %d tips, root-to-tip depth range %.4f - %.4f, rooted=%s (root has %d children)\n",
              ape::Ntip(left), min(left_depths), max(left_depths), ape::is.rooted(left), root_nchild(left)))
  cat(sprintf("[merge] right: %d tips, root-to-tip depth range %.4f - %.4f, rooted=%s (root has %d children)\n",
              ape::Ntip(right), min(right_depths), max(right_depths), ape::is.rooted(right), root_nchild(right)))

  merged <- merge_dated_subtrees(left, right, MRCA_AGE, ROOT_AGE)

  # verification: every tip's root-to-tip path length should now land at
  # (ROOT_AGE - MRCA_AGE) + that side's own root-to-tip depth.
  merged_depths <- ape::node.depth.edgelength(merged)[1:ape::Ntip(merged)]
  cat(sprintf("[merge] merged: %d tips, rooted=%s binary=%s, root-to-tip depth range %.4f - %.4f\n",
              ape::Ntip(merged), ape::is.rooted(merged), ape::is.binary(merged),
              min(merged_depths), max(merged_depths)))

  dir.create(dirname(OUT_NWK), showWarnings = FALSE, recursive = TRUE)
  ape::write.tree(merged, OUT_NWK)
  cat(sprintf("[merge] wrote %s\n", OUT_NWK))
}
