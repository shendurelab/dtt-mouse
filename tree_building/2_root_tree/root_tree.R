# root_tree.R
# Roots a raw NJ tree (basal trifurcation) on its side's synthetic founder-
# genotype root tip (SYNTHETIC_ROOT_B1 or SYNTHETIC_ROOT_B2, appended to the
# distance matrix by 1_build_nj_backbone/make_synthetic_root_outgroup.R) and
# drops that tip, producing the rooted, in-group-only backbone.
#
# Configuration via environment variables (all required, no defaults -- this
# script is generic over "which side"):
#   IN_NWK   : raw NJ newick (e.g. 1_build_nj_backbone's output for B1/B2)
#   OUTGROUP : tip to root on and drop (SYNTHETIC_ROOT_B1 or SYNTHETIC_ROOT_B2)
#   OUT_NWK  : where to write the rooted, outgroup-dropped newick
# e.g.
#   IN_NWK=results/B1/e3v5v6_nj_ge7.nwk OUTGROUP=SYNTHETIC_ROOT_B1 \
#   OUT_NWK=results/B1/nj_rooted_ingroup.nwk Rscript 2_root_tree/root_tree.R

IN_NWK   <- Sys.getenv("IN_NWK")
OUTGROUP <- Sys.getenv("OUTGROUP")
OUT_NWK  <- Sys.getenv("OUT_NWK")
if (!nzchar(IN_NWK))   stop("IN_NWK must be set (path to the raw NJ newick)")
if (!nzchar(OUTGROUP)) stop("OUTGROUP must be set (SYNTHETIC_ROOT_B1 or SYNTHETIC_ROOT_B2)")
if (!nzchar(OUT_NWK))  stop("OUT_NWK must be set (output path for the rooted, dropped newick)")
if (!file.exists(IN_NWK)) stop(sprintf("IN_NWK not found: %s", IN_NWK))

cat(sprintf("[root_tree] reading %s, rooting on %s ...\n", IN_NWK, OUTGROUP))

# NJ output is UNROOTED (basal trifurcation). resolve.root=TRUE puts a
# dichotomous root at the synthetic founder-genotype tip (the biologically
# correct root, editing being monotonic from that founder state), then that
# tip is dropped since it's synthetic, not a real cell.
tree <- ape::read.tree(IN_NWK)
stopifnot(OUTGROUP %in% tree$tip.label)
tree <- ape::root(tree, outgroup = OUTGROUP, resolve.root = TRUE)
tree <- ape::drop.tip(tree, OUTGROUP)
stopifnot(!(OUTGROUP %in% tree$tip.label))

cat(sprintf("[root_tree] rooted: %d tips, rooted=%s binary=%s (outgroup dropped)\n",
            ape::Ntip(tree), ape::is.rooted(tree), ape::is.binary(tree)))

dir.create(dirname(OUT_NWK), showWarnings = FALSE, recursive = TRUE)
ape::write.tree(tree, OUT_NWK)
cat(sprintf("[root_tree] wrote %s\n", OUT_NWK))
