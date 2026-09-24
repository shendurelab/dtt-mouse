# 01_parsimony_lib.R
# Postorder "prefix max parsimony" ancestral-state reconstruction for
# sequential lineage recorders (issue #21).
#
# Ported, with error handling simplified to plain stop() (no rlang
# dependency, matching this repo's style), from
# /Users/soseidel/Projects/phyloGramInf/R/parsimony.R -- an already-tested
# package function, not re-derived here. Only pars_ancestral_state() and its
# three tree-helper dependencies are ported; the traversal driver itself
# (02_support_lib.R's reconstruct_states_and_edits()) is new, because it also
# needs to count edits per branch in the same pass -- see that file.
#
# Sequential recorders write edits left-to-right into a fixed-length "tape".
# Once a lineage diverges at position k, positions k+1..L cannot be resolved
# by shared ancestry, so pars_ancestral_state() keeps the common prefix and
# collapses the tail to the unedited baseline "0". Missing data (a whole tape
# not recovered for a cell) is encoded as all-"?" and passed through.

# Root node id: the one node that never appears as an edge child.
.get_root <- function(phy) {
  root <- setdiff(unique(phy$edge[, 1]), unique(phy$edge[, 2]))
  if (length(root) != 1L) {
    stop(sprintf("Tree has %d candidate root node(s) (expected exactly 1).",
                 length(root)), call. = FALSE)
  }
  root
}

# Position-indexed children list: kids[[node]] is the integer vector of
# node's children (NULL for tips). split() does the O(n) grouping once, so
# later lookups are O(1) positional access.
.build_kids <- function(phy, n_total) {
  kids <- vector("list", n_total)
  grouped <- split(phy$edge[, 2], phy$edge[, 1])
  kids[as.integer(names(grouped))] <- grouped
  kids
}

# Internal nodes in post-order (descendants before ancestors). ape::postorder()
# returns edge indices; every internal node appears exactly once as edge[, 2]
# by the time all edges in its subtree have been emitted. The root never
# appears as a child, so it is appended before filtering to internal nodes.
.internal_postorder <- function(phy, root, Ntip) {
  po_edges <- ape::postorder(phy)
  nodes    <- c(phy$edge[po_edges, 2], root)
  nodes[nodes > Ntip]
}

#' Parsimony state at a bifurcation for a sequential recorder
#'
#' Given two child states `cs1` and `cs2` of equal length, returns the most
#' parsimonious parent state under a sequential-editing model: the longest
#' common prefix is retained, and all positions from the first mismatch
#' onward are collapsed to the unedited baseline `"0"`. Missing-data
#' children (all entries equal to `"?"`) are treated as unknown and the
#' non-missing sibling is passed through.
#'
#' @param cs1,cs2 Character vectors of equal length: per-site state of the
#'   two children. Use `"?"` for missing data, not `NA`.
#' @return A character vector of the same length: the reconstructed parent state.
pars_ancestral_state <- function(cs1, cs2) {
  if (!is.character(cs1) || !is.character(cs2)) {
    stop("`cs1` and `cs2` must be character vectors.", call. = FALSE)
  }
  if (length(cs1) != length(cs2)) {
    stop("`cs1` and `cs2` must have equal length.", call. = FALSE)
  }

  tape_length <- length(cs1)
  mismatch    <- which(cs1 != cs2)

  if (length(mismatch) == 0L) return(cs1)

  # All positions differ: treat all-"?" children as missing, otherwise no
  # common prefix exists and the parent must pre-date any edit -- baseline "0".
  if (length(mismatch) == tape_length) {
    if (all(cs1 == "?")) return(cs2)
    if (all(cs2 == "?")) return(cs1)
    return(rep("0", tape_length))
  }

  out <- cs1
  out[mismatch[1]:tape_length] <- "0"
  out
}
