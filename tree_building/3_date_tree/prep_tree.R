# prep_tree.R
# LSD2 reads branch lengths as substitution counts, so they must be
# non-negative -- but the NJ tree carries a few percent negative artefact
# edges. Corrects those with BAT::tree.zero() (Felsenstein-Kuhner
# redistribution: zeroes a negative edge while shortening its two child
# edges by the same total amount, preserving tip-to-tip distances), which
# gives equal-or-lower cophenetic distortion than naively clamping negative
# edges to 0 on this data.
#
#   input : a rooted, outgroup-dropped newick (env LSD_TREE -- 2_root_tree's
#           output)
#   output: <LSD_OUTDIR>/<tag>/tree_nonneg.nwk  (tag from tip count)
#
# Run from the repo root:
#   LSD_TREE=results/2-rooted-nj/B1/nj_rooted_ingroup.nwk \
#   LSD_OUTDIR=results/3-dated-tree/B1/lsd2 \
#     Rscript 3_date_tree/prep_tree.R

suppressPackageStartupMessages({ library(ape); library(BAT) })

IN_TREE <- Sys.getenv("LSD_TREE")
OUTDIR  <- Sys.getenv("LSD_OUTDIR")
if (!nzchar(IN_TREE)) stop("LSD_TREE must be set (rooted, outgroup-dropped newick)")
if (!nzchar(OUTDIR))  stop("LSD_OUTDIR must be set (output namespace dir)")

# 1. read the rooted tree.
message("[prep] reading ", IN_TREE)
tree <- read.tree(IN_TREE)
ntip <- ape::Ntip(tree)
message(sprintf("[prep] %d tips, %d edges", ntip, nrow(tree$edge)))

# 2. correct negative branch lengths (NJ artefacts; not valid substitution counts).
n_neg <- sum(tree$edge.length < 0)
tree  <- BAT::tree.zero(tree)
message(sprintf("[prep] BAT::tree.zero()-corrected %d negative edges (%.2f%%)", n_neg, 100 * n_neg / nrow(tree$edge)))

# 3. sanity: LSD2 needs a rooted, resolved topology (polytomies from the correction are fine).
stopifnot(ape::is.rooted(tree))
if (!ape::is.binary(tree)) message("[prep] note: tree is not fully binary (expected if any edges collapsed)")
depths <- ape::node.depth.edgelength(tree)[seq_len(ntip)]
message(sprintf("[prep] root->tip substitution depth range: %.3f .. %.3f", min(depths), max(depths)))

# 4. write the cleaned newick into a per-tree-size subfolder (so full and subtree
#    runs never clash).
tag    <- if (ntip > 50000) "nj99478" else sprintf("sub%dk", round(ntip / 1000))
subdir <- file.path(OUTDIR, tag)
dir.create(subdir, showWarnings = FALSE, recursive = TRUE)
out    <- file.path(subdir, "tree_nonneg.nwk")
write.tree(tree, out)
message(sprintf("[prep] wrote %s (%d tips, tag=%s)", out, ntip, tag))
