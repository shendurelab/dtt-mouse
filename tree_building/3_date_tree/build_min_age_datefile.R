#!/usr/bin/env Rscript
# =============================================================================
# build_min_age_datefile.R -- build an LSD2 "-d" date file that caps the LTT
# (lineages-through-time) at a literature-derived embryo cell-count ceiling.
#
# THE PROBLEM: plain LSD2 reads "zero edits on a branch" as "~zero elapsed
# time", so it packs hundreds of early divisions into a sliver of calendar
# time just after the root. The resulting time tree can have MORE sampled
# ancestral lineages at ~E5-E7 than the embryo had cells -- physically
# impossible for a subsample's genealogy.
#
# THE FIX ENCODED HERE: a literature-derived LADDER of minimum node ages. A
# subsample's genealogy can have at most N(d) lineages at day d (N = embryo
# cell count). So the split that raises the tree's lineage count to K cannot
# have happened earlier than the first day the embryo had K cells. We turn
# that into per-node lower-age bounds `l(day)` in LSD2's native `-d` date file
# and let LSD2 re-date the whole tree under its own constrained least squares.
#
#
# INPUTS (positional args):
#   1 DATED_TREE   per-side initial LSD2 time tree (time_tree.nwk; branch
#                  lengths in days, root at ROOT_DAY, tips at TIP_DAY). Used
#                  ONLY to read the initial age RANK and to name nodes by two
#                  representative tips -- the bounds themselves come from the
#                  ceiling, not from these ages.
#   2 CEILING_CSV  whole-embryo cell-count ceiling (processed_data's
#                  sample_matched_ceiling_sourced.csv: columns day, cells_min,
#                  cells_max, basis).
#   3 OUT_DATEFILE path to write the LSD2 date file.
#   4 SIDE_FRAC    fraction of the whole embryo this tree represents (default
#                  0.5). B1/B2 are the two first-division daughters but need
#                  not split the embryo evenly -- pass per-side fractions
#                  (e.g. 0.63/0.37) when that's better supported than 0.5/0.5.
#
# ENV (defaults match the per-side pipeline):
#   ROOT_DAY  (default 1.5)  LSD_ROOT used in run_dating.sh
#   CEIL_COL  (default cells_max)  which ceiling column to invert. cells_max is
#             the WEAKEST / most conservative bound: we only forbid lineage
#             counts impossible even under the most generous cell count.
#
# OUTPUTS:
#   OUT_DATEFILE                 the LSD2 -d file (line 1 = #constraints).
#   <OUT_DATEFILE>.audit.csv     one row per constrained node for eyeballing:
#                                node, n_children, cum_lineages, min_day,
#                                cur_day, binding, tipA, tipB.
# =============================================================================

suppressPackageStartupMessages(library(ape))

args     <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3)
  stop("usage: build_min_age_datefile.R DATED_TREE CEILING_CSV OUT_DATEFILE [SIDE_FRAC]")
DATED_TREE   <- args[1]
CEILING_CSV  <- args[2]
OUT_DATEFILE <- args[3]
SIDE_FRAC    <- if (length(args) >= 4) as.numeric(args[4]) else 0.5
ROOT_DAY     <- as.numeric(Sys.getenv("ROOT_DAY", "1.5"))
CEIL_COL     <- Sys.getenv("CEIL_COL", "cells_max")

# ---- 1. per-side cell-count ceiling, made monotone (embryos only grow) ------
ceil <- read.csv(CEILING_CSV, stringsAsFactors = FALSE)
ceil <- ceil[order(ceil$day), ]
days_land  <- ceil$day
# per-side landmark counts; cummax so a later, smaller literature number never
# lowers the ceiling below an earlier one (different literature sources use
# different cell-population definitions, which is not real shrinkage).
cells_land <- cummax(ceil[[CEIL_COL]] * SIDE_FRAC)
CAP <- cells_land[length(cells_land)]   # above this the ceiling is off-table => no bound

# t_min(m): earliest day the (per-side) embryo first had >= m cells. Log-linear
# in log(cells) between landmarks (exponential growth looks linear there).
t_min <- function(m) {
  if (m <= cells_land[1]) return(days_land[1])
  if (m >  CAP)           return(Inf)          # off-table => leave node free
  for (i in 2:length(cells_land)) {
    if (cells_land[i] >= m) {
      c0 <- cells_land[i - 1]; c1 <- cells_land[i]
      if (c1 == c0) return(days_land[i])       # flat segment: reached at its right end
      f <- (log(m) - log(c0)) / (log(c1) - log(c0))
      return(days_land[i - 1] + f * (days_land[i] - days_land[i - 1]))
    }
  }
  Inf
}

# ---- 2. read the initial per-side dated tree -------------------------------
tr <- read.tree(DATED_TREE)
ntip <- length(tr$tip.label)
message(sprintf("[build] %s: %d tips, %d internal nodes",
                basename(DATED_TREE), ntip, tr$Nnode))

# current day of every node = ROOT_DAY + (root-to-node time depth)
depth   <- node.depth.edgelength(tr)          # root = 0
cur_day <- ROOT_DAY + depth                    # indexed 1..(ntip+Nnode)

# children of each internal node, in tree order
kids <- split(tr$edge[, 2], tr$edge[, 1])      # names = parent node ids (chars)

# ---- 3. one representative tip label per node (iterative post-order) --------
# rep_tip[v] = some leaf under v. Post-order guarantees children are filled
# before parents, so rep_tip[parent] <- rep_tip[first child seen].
po      <- reorder(tr, "postorder")$edge       # child rows precede parent rows
rep_tip <- rep(NA_integer_, ntip + tr$Nnode)
rep_tip[seq_len(ntip)] <- seq_len(ntip)        # tips represent themselves
for (e in seq_len(nrow(po))) {
  parent <- po[e, 1]; child <- po[e, 2]
  if (is.na(rep_tip[parent])) rep_tip[parent] <- rep_tip[child]
}

# ---- 4. rank internal nodes oldest->youngest; cumulative lineage count ------
# Just after the j-th split (counting from the root), the tree carries
# 1 + sum_{split s<=j} (n_children(s) - 1) coexisting lineages. That running
# total is the K we feed to t_min().
root_node <- setdiff(tr$edge[, 1], tr$edge[, 2])   # the one node with no parent
stopifnot(length(root_node) == 1)
inodes <- (ntip + 1):(ntip + tr$Nnode)
inodes <- inodes[order(cur_day[inodes])]       # oldest (smallest day) first
n_child  <- vapply(inodes, function(v) length(kids[[as.character(v)]]), integer(1))
cum_lin  <- 1L + cumsum(n_child - 1L)          # lineages present just after node v splits

# ---- 5. emit a lower-bound line per node, skipping vacuous/off-table ones ---
rows <- vector("list", length(inodes))
for (k in seq_along(inodes)) {
  v <- inodes[k]
  K <- cum_lin[k]
  md <- t_min(K)
  # skip if: off-table (Inf), or the bound is at/below the fixed root age
  # (vacuous -- every node is already >= ROOT_DAY), or v is the root itself
  # (the root's date is fixed by LSD2's -a; any bound on it is a hard conflict).
  if (v == root_node || !is.finite(md) || md <= ROOT_DAY + 1e-9) next
  ch <- kids[[as.character(v)]]
  tipA <- tr$tip.label[rep_tip[ch[1]]]
  tipB <- tr$tip.label[rep_tip[ch[2]]]          # different child => mrca(tipA,tipB) == v
  rows[[k]] <- data.frame(
    node = v, n_children = length(ch), cum_lineages = K,
    min_day = round(md, 6), cur_day = round(cur_day[v], 6),
    binding = md > cur_day[v] + 1e-9,           # does the bound actually push this node later?
    tipA = tipA, tipB = tipB, stringsAsFactors = FALSE)
}
audit <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])

if (is.null(audit) || nrow(audit) == 0)
  stop("No constraints generated -- check ceiling / ROOT_DAY.")

# ---- 6. write the LSD2 date file + the audit table -------------------------
dir.create(dirname(OUT_DATEFILE), showWarnings = FALSE, recursive = TRUE)
lines <- c(as.character(nrow(audit)),
           sprintf("mrca(%s,%s) l(%s)", audit$tipA, audit$tipB, format(audit$min_day, trim = TRUE)))
writeLines(lines, OUT_DATEFILE)
write.csv(audit, paste0(OUT_DATEFILE, ".audit.csv"), row.names = FALSE)

# ---- 7. one-screen summary for the human --------------------------------
message(sprintf("[build] ceiling column     : %s  (x side_frac=%.2f, cummax; CAP=%g cells/side)",
                CEIL_COL, SIDE_FRAC, CAP))
message(sprintf("[build] constrained nodes   : %d of %d internal (cum_lineages in [2, %g])",
                nrow(audit), tr$Nnode, CAP))
message(sprintf("[build] BINDING (bound pushes node later than its current age): %d",
                sum(audit$binding)))
message(sprintf("[build] min_day range       : %.3f .. %.3f", min(audit$min_day), max(audit$min_day)))
message(sprintf("[build] wrote %s  (+ .audit.csv)", OUT_DATEFILE))
