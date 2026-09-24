# 02_support_lib.R
# New, mouse_sprint-specific layer on top of 01_parsimony_lib.R's ported
# pars_ancestral_state(): reconstructs per-tape ancestral tape states AND
# counts, per node, how many sites show an edit-state change on the branch
# leading into that node -- the "edit-count support" issue #21 asks for.
#
# Design: a single postorder pass (reconstruct_states_and_edits()) computes
# the parent state via pars_ancestral_state(cs1, cs2) and immediately counts
# edits for BOTH children right there, since cs1, cs2, and the freshly
# derived parent_state are all already in scope at that point -- no second
# traversal or parent-lookup table is needed. The root (no parent edge) is
# compared against a synthetic, fully-unedited "germline" state once the
# loop finishes.
#
# Run (example, from repo root):
#   source("src/8-support/01_parsimony_lib.R")
#   source("src/8-support/02_support_lib.R")
#   result <- reconstruct_support(tree, tip_states_long, barcodes)

# The 6 sequential sites of every tape. The single source of truth for this
# constant -- 03_tape_tip_state_mapping.R and any Rmd report reference this
# copy rather than re-literalizing paste0("Site", 1:6), so a change to the
# site count only has to happen in one place. (Every entry point in this repo
# sources 01 -> 02 -> 03 in that order, so 03 can rely on this being defined.)
SITE_COLS <- paste0("Site", 1:6)

#' Count edits on the branch from a resolved parent state to a child state
#'
#' A site counts as an edit only if the parent is a RESOLVED unedited state
#' ("0") and the child is a RESOLVED edit token (not "0", not "?"). Any "?"
#' on either side (missing data, however it arose) contributes 0 -- this is
#' the "conservative strategy for handling missing tapes" issue #21 asks
#' for, and it falls out of this one formula with no extra branching.
#'
#' @param parent_state,child_state Character vectors of equal length.
#' @return Integer: number of sites edited on this branch.
count_edits <- function(parent_state, child_state) {
  sum(parent_state == "0" & !(child_state %in% c("0", "?")))
}

#' Reconstruct ancestral tape states and per-node edit counts for one tape
#'
#' Single postorder pass over `phy`: reconstructs the parsimony state at
#' every internal node via [pars_ancestral_state()], and for every node
#' (tip or internal) counts edits on the branch leading into it (the root's
#' branch is measured against a synthetic all-"0" germline state).
#'
#' @param phy        An `ape::phylo` object. Must be rooted and binary.
#' @param tip_states A data frame with an id column (default `"cell_id"`)
#'   matching `phy$tip.label`, and one column per site.
#' @param site_cols  Site column names in `tip_states`. Defaults to
#'   `"Site1".."Site6"`.
#' @param id_col     Id column in `tip_states`. Defaults to `"cell_id"`.
#' @return A list with:
#'   * `states`      -- (Ntip+Nnode) x length(site_cols) character matrix,
#'                      ape-numbered, the reconstructed tape state per node.
#'   * `edit_counts` -- integer vector of length Ntip+Nnode, ape-numbered,
#'                      edits on the branch leading into each node.
reconstruct_states_and_edits <- function(phy, tip_states,
                                          site_cols = SITE_COLS,
                                          id_col    = "cell_id") {
  if (!inherits(phy, "phylo")) stop("`phy` must be an ape::phylo object.", call. = FALSE)
  if (!ape::is.rooted(phy)) stop("`phy` must be a rooted tree.", call. = FALSE)
  if (!id_col %in% colnames(tip_states)) {
    stop(sprintf("`tip_states` has no column named '%s'.", id_col), call. = FALSE)
  }
  missing_cols <- setdiff(site_cols, colnames(tip_states))
  if (length(missing_cols)) {
    stop("`tip_states` is missing site column(s): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(tip_states[[id_col]])) {
    stop(sprintf("`tip_states` has duplicate values in '%s' column.", id_col), call. = FALSE)
  }

  Ntip    <- length(phy$tip.label)
  n_total <- Ntip + phy$Nnode
  L       <- length(site_cols)
  root    <- .get_root(phy)
  kids    <- .build_kids(phy, n_total)

  tip_idx <- match(phy$tip.label, tip_states[[id_col]])
  if (any(is.na(tip_idx))) {
    stop("`tip_states` has no rows for tip(s): ",
         paste(phy$tip.label[is.na(tip_idx)], collapse = ", "), call. = FALSE)
  }

  states <- matrix(NA_character_, nrow = n_total, ncol = L,
                    dimnames = list(NULL, site_cols))
  for (j in seq_along(site_cols)) {
    col_vals <- tip_states[[site_cols[j]]][tip_idx]
    if (anyNA(col_vals)) {
      stop(sprintf("`tip_states` has NA in column '%s' -- recode missing data as \"?\".",
                   site_cols[j]), call. = FALSE)
    }
    states[seq_len(Ntip), j] <- as.character(col_vals)
  }

  edit_counts <- integer(n_total)
  post_order  <- .internal_postorder(phy, root, Ntip)

  for (node in post_order) {
    node_children <- kids[[node]]
    if (length(node_children) != 2L) {
      stop(sprintf("Node %d has %d children; only binary trees are supported.",
                   node, length(node_children)), call. = FALSE)
    }
    c1 <- node_children[1L]; c2 <- node_children[2L]
    parent_state <- pars_ancestral_state(states[c1, ], states[c2, ])
    states[node, ]     <- parent_state
    edit_counts[c1]    <- count_edits(parent_state, states[c1, ])
    edit_counts[c2]    <- count_edits(parent_state, states[c2, ])
  }

  # Root has no real parent edge -- compare to the synthetic, fully-unedited
  # germline state (nothing has been edited before the founder cell).
  edit_counts[root] <- count_edits(rep("0", L), states[root, ])

  list(states = states, edit_counts = edit_counts)
}

#' Run reconstruct_states_and_edits() once per tape (barcode), in parallel,
#' and aggregate into the two issue #21 deliverables.
#'
#' @param phy             An `ape::phylo` object. Must be rooted and binary.
#' @param tip_states_long Long-format data frame: `id_col`, `barcode_col`,
#'   and one column per site (default `Site1..Site6`).
#' @param barcodes        Character vector of tape/barcode names to process.
#' @param site_cols       Site column names. Defaults to `SITE_COLS`.
#' @param id_col          Id column name. Defaults to `"cell_id"`.
#' @param barcode_col     Barcode column name in `tip_states_long`. Defaults
#'   to `"barcode"`, matching `src/4-time-tree/01_cell_by_tape.R`'s output.
#' @param mc.cores        Cores for `parallel::mclapply()` across tapes.
#' @return A list with:
#'   * `tape_states_long` -- data frame (node, barcode, `site_cols`...): the
#'     reconstructed tape state per node per tape (Deliverable A). Named
#'     "tape_states_long" for the deliverable it produces, but its identifier
#'     column is "barcode" throughout -- same name as the input, end to end.
#'   * `support_df` -- data frame (node, one column per tape/barcode,
#'     `support` = row sum): the edit-count support matrix (Deliverable B,
#'     matches the shape of `support.jpeg`'s table).
reconstruct_support <- function(phy, tip_states_long, barcodes,
                                 site_cols   = SITE_COLS,
                                 id_col      = "cell_id",
                                 barcode_col = "barcode",
                                 mc.cores    = max(1, parallel::detectCores() - 1)) {
  n_total <- length(phy$tip.label) + phy$Nnode

  tape_results <- parallel::mclapply(barcodes, function(bc) {
    tip_states_bc <- tip_states_long[tip_states_long[[barcode_col]] == bc, , drop = FALSE]
    reconstruct_states_and_edits(phy, tip_states_bc, site_cols = site_cols, id_col = id_col)
  }, mc.cores = mc.cores)
  names(tape_results) <- barcodes

  failed <- barcodes[vapply(tape_results, inherits, logical(1), what = "try-error")]
  if (length(failed)) {
    stop("reconstruct_states_and_edits() failed for tape(s): ",
         paste(failed, collapse = ", "), call. = FALSE)
  }

  tape_states_long <- do.call(rbind, lapply(barcodes, function(bc) {
    st <- tape_results[[bc]]$states
    data.frame(node = seq_len(n_total), barcode = bc, st,
               row.names = NULL, stringsAsFactors = FALSE)
  }))

  edit_count_matrix <- vapply(barcodes, function(bc) tape_results[[bc]]$edit_counts,
                               numeric(n_total))
  support <- rowSums(edit_count_matrix)
  support_df <- data.frame(node = seq_len(n_total), edit_count_matrix,
                            support = support, row.names = NULL,
                            stringsAsFactors = FALSE)

  list(tape_states_long = tape_states_long, support_df = support_df)
}

#' Map a subtree's own node ids back to the original tree's node ids
#'
#' `ape::extract.clade()` and `ape::keep.tip()` both renumber nodes in the
#' subtree they return, so its node ids don't line up with the original
#' tree's node ids that `reconstruct_support()`'s output is indexed by. Tips
#' match directly by label; each internal node is matched via the MRCA, in
#' the original tree, of its two already-resolved children, computed
#' bottom-up in one postorder pass (reusing the same `.build_kids()`/
#' `.internal_postorder()` helpers `01_parsimony_lib.R` ports) rather than
#' re-deriving each node's full descendant tip set from scratch -- O(n)
#' instead of the O(n^2)-ish cost of calling `extract.clade()` per node.
#'
#' Works for any subtree whose tips are a subset of `original`'s tips, not
#' just a monophyletic clade -- e.g. a root-spanning sample pruned with
#' `ape::keep.tip()` is not a clade, but this still works the same way.
#'
#' @param subtree  A `phylo` object (from `ape::extract.clade()`,
#'   `ape::keep.tip()`, etc.) whose tips are a subset of `original`'s tips.
#' @param original The tree `subtree` was derived from.
#' @return Integer vector, length `Ntip(subtree) + subtree$Nnode`,
#'   ape-numbered by `subtree`'s own node ids: `original`'s node id for each.
map_subtree_nodes_to_original <- function(subtree, original) {
  n_tip_subtree   <- length(subtree$tip.label)
  n_total_subtree <- n_tip_subtree + subtree$Nnode
  orig_id <- integer(n_total_subtree)
  orig_id[seq_len(n_tip_subtree)] <- match(subtree$tip.label, original$tip.label)

  root       <- .get_root(subtree)
  kids       <- .build_kids(subtree, n_total_subtree)
  post_order <- .internal_postorder(subtree, root, n_tip_subtree)
  for (node in post_order) {
    children <- kids[[node]]
    stopifnot("only binary subtrees are supported" = length(children) == 2L)
    # The MRCA of two already-resolved children IS the MRCA of their full
    # tip sets (a basic property of trees), so getMRCA() never needs to see
    # the subtree's actual leaves below this point.
    orig_id[node] <- ape::getMRCA(original, orig_id[children])
  }
  orig_id
}

#' Load a cached full-tree reconstruction if still fresh, else compute + save
#'
#' `reconstruct_support()` on the full tree takes several seconds; every
#' entry point that needs it (both Rmd reports, the Taxonium export script)
#' shares one cache file instead of recomputing per render. Freshness is
#' checked by md5 of the input FILES themselves (not their parsed contents),
#' so `compute_fn` never has to load the tape data just to find out the
#' cache is still valid -- only on an actual cache miss.
#'
#' @param cache_path  Path to the cached `.rds` (created if missing).
#' @param input_paths Character vector of file paths this reconstruction
#'   depends on (the tree Newick, the cell-by-tape RDS). Whenever any of
#'   these files' contents change, the cache is treated as stale and
#'   recomputed -- a changed upstream input can't silently keep serving an
#'   outdated reconstruction forever, unlike a plain `file.exists()` check.
#' @param compute_fn  Zero-argument function that loads whatever
#'   `input_paths` point to and returns a fresh `reconstruct_support()`
#'   result. Only called on a cache miss/staleness, so callers can defer
#'   loading large tape data until they know it's actually needed.
#' @return Whatever `compute_fn()` returns (a `reconstruct_support()`-shaped list).
get_or_compute_full_tree_reconstruction <- function(cache_path, input_paths, compute_fn) {
  input_md5 <- vapply(input_paths, function(p) unname(tools::md5sum(p)), character(1))

  # cat(), not message(): Rmd reports in this repo set message=FALSE at the
  # chunk-option level, which would silently swallow this status line if it
  # went to the message stream -- cat() is the same choice the code this
  # replaces made, for the same reason.
  if (file.exists(cache_path)) {
    # A cache file can be unreadable for reasons that have nothing to do
    # with staleness (an interrupted write, disk truncation, etc.) --
    # readRDS() throws in that case, which should be treated the same as a
    # missing/stale cache (recompute and overwrite), not a hard failure.
    cached <- tryCatch(readRDS(cache_path), error = function(e) {
      cat("[cache]", cache_path, "is unreadable (", conditionMessage(e),
          ") -- recomputing\n")
      NULL
    })
    if (!is.null(cached) && identical(cached$input_md5, input_md5)) {
      cat("[cache] loaded full-tree reconstruction from", cache_path, "\n")
      return(cached$result)
    }
    if (!is.null(cached)) {
      cat("[cache]", cache_path, "is stale (input file(s) changed) -- recomputing\n")
    }
  }

  t0 <- Sys.time()
  result <- compute_fn()
  cat("[cache] reconstruct_support() took",
      round(difftime(Sys.time(), t0, units = "secs"), 1), "s\n")

  dir.create(dirname(cache_path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(list(result = result, input_md5 = input_md5), cache_path)
  result
}

#' % of internal nodes with support >= threshold, aggregated into day bins
#'
#' Shared by both the full-tree report (irregular bins spanning the whole
#' 0-13.5 day range) and the subclade demo report (regular 2-day bins over
#' just that clade's own span) -- same aggregation either way, just
#' different bin edges and, for a subtree, a day offset. Internal nodes
#' only: tips sit at the sampling endpoint, not partway through the
#' branching history, so they wouldn't add time resolution (same convention
#' as the ancestral cell-type reconstruction report).
#'
#' Plotting is deliberately NOT part of this function -- the two reports
#' style their support-over-time charts differently on purpose (one
#' restrained/monochrome, one the project's default blue), so only the data
#' step is shared; each report builds its own ggplot from this data frame.
#'
#' @param tree           A `phylo`. For the full tree this is the tree
#'   itself (root is already day 0); for a subtree from `ape::extract.clade()`
#'   /`ape::keep.tip()`, pass `day_offset` too (see below).
#' @param support        Numeric vector, ape-node-indexed (length
#'   `Ntip(tree)+tree$Nnode`), e.g. `support_df$support` or
#'   `demo_support$support` -- MUST be in `tree`'s own node order.
#' @param bin_edges_days Bin edges (N+1 edges for N bins), in the same
#'   absolute-day reference frame `day_offset` puts `tree`'s own node depths
#'   into.
#' @param day_offset     Added to `tree`'s own node depths before binning.
#'   0 (default) for a full tree, whose root already IS day 0.
#'   `ape::extract.clade()`/`ape::keep.tip()` reset a subtree's root to x=0,
#'   so a subtree needs its root's absolute day (from the ORIGINAL tree)
#'   here instead.
#' @param threshold      Support threshold. Defaults to 2 (issue #21's own
#'   pivot value, also used for the Taxonium color scale).
#' @return data.frame(bin, n_total, n_supported, pct_supported, mid_day,
#'   width) -- one row per bin, in `bin_edges_days` order.
compute_support_over_time <- function(tree, support, bin_edges_days, day_offset = 0, threshold = 2) {
  Ntip <- length(tree$tip.label)
  internal_idx <- (Ntip + 1):(Ntip + tree$Nnode)
  absolute_day <- day_offset + ape::node.depth.edgelength(tree)[internal_idx]
  supported <- support[internal_idx] >= threshold

  bin_lo <- head(bin_edges_days, -1)
  bin_hi <- tail(bin_edges_days, -1)
  bin_labels <- paste0(bin_lo, "-", bin_hi)
  bin <- cut(absolute_day, breaks = bin_edges_days, labels = bin_labels, include.lowest = TRUE)
  bin <- factor(bin, levels = bin_labels)

  df <- data.frame(bin = bin_labels,
                    n_total = as.integer(table(bin)),
                    n_supported = as.integer(tapply(supported, bin, sum)))
  df$n_supported[is.na(df$n_supported)] <- 0L
  df$pct_supported <- ifelse(df$n_total > 0, 100 * df$n_supported / df$n_total, NA_real_)
  df$mid_day <- (bin_lo + bin_hi) / 2
  df$width   <- bin_hi - bin_lo
  df
}

#' Per (blastomere, day-bin): internal-node counts clearing support thresholds
#'
#' Shared by 18 (support = edits summed across tapes) and 19 (support =
#' distinct tapes edited) -- same binning and thresholding, different input
#' vector. Excludes the root (no incoming branch) and multi2di's zero-length
#' branches (not real nodes).
#'
#' @param tree        Binarised tree matching `support_df`'s node ids.
#' @param support_df  08 reconstruction cache (`node`, `embryonic_day`,
#'   `blastomere`); its `node` order fixes the order `support_vec` must be in.
#' @param support_vec Per-node value to threshold, in `support_df$node` order.
#' @param bin_edges   Day-bin edges (N+1 edges for N bins).
#' @param thresholds  Named integer vector: display name -> threshold (>=).
#' @return data.frame(blastomere, bin, day_mid, n_internal_nodes,
#'   n_ge<thr>/pct_ge<thr> per threshold) -- one row per (blastomere, bin).
support_over_time_table <- function(tree, support_df, support_vec, bin_edges, thresholds) {
  Ntip <- length(tree$tip.label); N <- Ntip + tree$Nnode
  stopifnot("support_df rows != tree nodes" = nrow(support_df) == N,
            "support_df not in node order" = all(support_df$node == seq_len(N)),
            "support_vec length != tree nodes" = length(support_vec) == N)
  day <- support_df$embryonic_day

  parent_of <- integer(N); parent_of[tree$edge[, 2]] <- tree$edge[, 1]  # 0 for the root
  is_internal <- seq_len(N) > Ntip
  # only index has-parent positions -- day[0] (the root) would misalign the vector
  has_parent <- parent_of > 0
  branch_len <- rep(NA_real_, N)
  branch_len[has_parent] <- day[has_parent] - day[parent_of[has_parent]]

  nodes <- data.frame(support = support_vec, day = day,
                       blastomere = support_df$blastomere,
                       branch_len = branch_len, internal = is_internal)
  keep <- nodes$internal & !is.na(nodes$blastomere) &
          !is.na(nodes$branch_len) & nodes$branch_len > 1e-9
  nodes <- nodes[keep, ]
  nodes$bin <- cut(nodes$day, breaks = bin_edges, include.lowest = TRUE)

  lo <- head(bin_edges, -1); hi <- tail(bin_edges, -1)
  grid <- expand.grid(blastomere = c("A", "B"), bin = levels(nodes$bin),
                      stringsAsFactors = FALSE)
  n_tot <- aggregate(support ~ blastomere + bin, data = nodes, FUN = length)
  names(n_tot)[3] <- "n_internal_nodes"
  tab <- merge(grid, n_tot, all.x = TRUE)
  for (tn in names(thresholds)) {
    thr <- thresholds[[tn]]
    nge <- aggregate(support ~ blastomere + bin, data = nodes,
                     FUN = function(s) sum(s >= thr))
    names(nge)[3] <- paste0("n_ge", thr)
    tab <- merge(tab, nge, all.x = TRUE)
  }
  tab$n_internal_nodes[is.na(tab$n_internal_nodes)] <- 0L
  out_cols <- c("blastomere", "bin", "day_mid", "n_internal_nodes")
  for (tn in names(thresholds)) {
    thr <- thresholds[[tn]]
    n_col <- paste0("n_ge", thr); pct_col <- paste0("pct_ge", thr)
    tab[[n_col]][is.na(tab[[n_col]])] <- 0L
    tab[[pct_col]] <- ifelse(tab$n_internal_nodes > 0,
                             100 * tab[[n_col]] / tab$n_internal_nodes, NA_real_)
    out_cols <- c(out_cols, n_col, pct_col)
  }
  tab$day_mid <- ((lo + hi) / 2)[match(tab$bin, levels(nodes$bin))]
  tab[order(tab$blastomere, tab$day_mid), out_cols]
}
