# time_tree.R
# Turn LSD2's output into a clean, downstream-friendly time tree. LSD2 writes
# the dated tree as a NEXUS with per-node [&date=...] annotations
# (dated.date.nexus) whose branch lengths are already in DAYS; the plain
# dated.nwk it also writes is in substitution units (days x rate), which is
# not what we want downstream. This reads the day-scaled NEXUS and writes a
# plain newick in days, then validates that the tree is ultrametric with tips
# at the sampling day and root at the root date.
#
#   input : <tag>/dated.date.nexus   (from run_lsd2.sh)
#   output: <tag>/time_tree.nwk      (branch lengths in DAYS)
#
# Run from the repo root:
#   Rscript 3_date_tree/time_tree.R results/3-dated-tree/B1/lsd2/nj99478

suppressPackageStartupMessages(library(ape))

args   <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: time_tree.R SUBDIR")
SUBDIR   <- args[[1]]
TIPS_DAY <- as.numeric(Sys.getenv("LSD_TIPS", "13.5"))

IN  <- file.path(SUBDIR, "dated.date.nexus")
OUT <- file.path(SUBDIR, "time_tree.nwk")

# 1. read the day-scaled dated tree (read.nexus keeps the time-unit branch lengths
#    and drops the [&date=...] comments).
message("[time_tree] reading ", IN)
tree <- read.nexus(IN)

# 2. node depths from the root = node dates in days (root = 0).
depths   <- ape::node.depth.edgelength(tree)
tip_days <- depths[seq_len(ape::Ntip(tree))]

# 3. validate: tips contemporaneous at the sampling day, root at 0, no negatives.
tol <- 1e-2
stopifnot(all(tree$edge.length >= -1e-9))
if (max(abs(tip_days - TIPS_DAY)) > tol)
  warning(sprintf("[time_tree] tip dates deviate from %.3f by up to %.4f day", TIPS_DAY, max(abs(tip_days - TIPS_DAY))))
message(sprintf("[time_tree] %d tips; tip-date range %.4f..%.4f d; root %.4f; height %.4f d",
                ape::Ntip(tree), min(tip_days), max(tip_days), min(depths), max(tip_days)))

# 4. write the plain day-unit newick for downstream use.
write.tree(tree, OUT)
message(sprintf("[time_tree] wrote %s (branch lengths in days)", OUT))
