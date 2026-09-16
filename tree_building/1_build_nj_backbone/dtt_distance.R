# dtt_distance.R
# Typewriter-aware cell-cell distance for DNA Typewriter (TAPE) recordings.
#
# DTT site conventions (one value per site, Site1..Site6 of a tape):
#   - a 3-char trinucleotide, e.g. "CTT"  -> a lineage EDIT
#   - "ETY"                               -> empty / UNEDITED site
#   - "None" or "?"                       -> site NOT detected (missing)
#
# A cell is a named list of integrations (tapes); each tape is a length-6
# character vector of site states. The list name is the integration id, so the
# same id in two cells refers to the same genomic integration.
#
# Distance (adapted from prior DNA-Typewriter work):
#   For every integration recovered in BOTH cells, the per-integration distance
#   is the number of edit steps to turn one cell's 6-site array into the other:
#   reverse the unshared edits back to the common ancestor, then apply the
#   forward edits, i.e.
#
#       depth1 + depth2 - 2 * shared_edit_prefix          (range 0..12)
#
#   where depth = number of edited sites on the tape and shared_edit_prefix =
#   length of the leading run of sites that are the SAME edit in both cells
#   (the sequential SciPhy editing model: edits accumulate from Site1 onward).
#   The cell-cell distance is the AVERAGE of this over shared integrations.

MISSING_STATES <- c("None", "?")
UNEDITED_STATE <- "ETY"

# A site is an edit iff it is neither missing nor the unedited placeholder.
is_edit <- function(site) {
  !(site %in% MISSING_STATES) & site != UNEDITED_STATE
}

# depth: number of edited sites on a tape.
tape_depth <- function(sites) {
  sum(is_edit(sites))
}

# shared_edit_prefix: length of the leading run where both tapes carry the
# SAME edit. Stops at the first site that is unedited/missing in either tape or
# where the two edits differ (sequential model => no edits beyond the break).
shared_edit_prefix <- function(sites1, sites2) {
  k <- 0L
  for (i in seq_along(sites1)) {
    if (is_edit(sites1[i]) && is_edit(sites2[i]) && sites1[i] == sites2[i]) {
      k <- k + 1L
    } else {
      break
    }
  }
  k
}

# Per-integration distance: depth1 + depth2 - 2 * shared_edit_prefix.
tape_distance <- function(sites1, sites2) {
  tape_depth(sites1) + tape_depth(sites2) - 2L * shared_edit_prefix(sites1, sites2)
}

# A tape is "recovered" in a cell if it is present and not entirely missing
# (an all-"None" tape is the DTT convention for an undetected integration).
tape_recovered <- function(sites) {
  !is.null(sites) && !all(sites %in% MISSING_STATES)
}

# Cell-cell distance: mean per-integration distance over integrations recovered
# in BOTH cells. Returns NA if the two cells share no recovered integration.
dtt_distance <- function(cell1, cell2) {
  tapes1 <- names(cell1)[vapply(cell1, tape_recovered, logical(1))]
  tapes2 <- names(cell2)[vapply(cell2, tape_recovered, logical(1))]
  shared <- intersect(tapes1, tapes2)
  if (length(shared) == 0L) return(NA_real_)
  d <- vapply(shared, function(tp) tape_distance(cell1[[tp]], cell2[[tp]]), numeric(1))
  mean(d)
}

# Convenience: full cell x cell distance matrix from a named list of cells.
dtt_distance_matrix <- function(cells) {
  ids <- names(cells)
  n <- length(ids)
  m <- matrix(0, n, n, dimnames = list(ids, ids))
  for (i in seq_len(n)) for (j in i:n) {
    d <- dtt_distance(cells[[i]], cells[[j]])
    m[i, j] <- d
    m[j, i] <- d
  }
  m
}
