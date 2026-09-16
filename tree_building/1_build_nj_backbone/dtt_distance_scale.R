# dtt_distance_scale.R
# Fast, exactly-equivalent reimplementation of dtt_distance_matrix() from
# dtt_distance.R: the symmetric cell x cell matrix whose [i,j] is the mean over
# integrations recovered in both cells of (depth_i + depth_j - 2*shared_prefix),
# NA when the two cells share no recovered integration. It encodes cells once
# (edit tokens -> integer ids; per-tape depths; a recovered bitmask per cell)
# then runs a RcppParallel kernel over row blocks.
#
# Run:
#   source("1_build_nj_backbone/dtt_distance_scale.R")
#   m <- dtt_distance_matrix_fast(cells)     # drop-in for dtt_distance_matrix
# or reuse the encoding across calls:
#   enc <- encode_cells(cells); m <- dtt_distance_matrix_encoded(enc)

# Robust path to the sidecar .cpp, then compile (Rcpp caches by file mtime/hash).
.dtt_cpp_path <- function() {
  # When source()'d, sys.frame(1)$ofile (or the sourced file) gives our location.
  for (i in seq_len(sys.nframe())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(file.path(dirname(normalizePath(of)),
                                       "dtt_distance_scale.cpp"))
  }
  # Fallback: assume the conventional repo layout.
  "1_build_nj_backbone/dtt_distance_scale.cpp"
}

Rcpp::sourceCpp(.dtt_cpp_path())

# is_edit / tape_recovered semantics, applied as a vectorised encoder.
# Token convention: not-an-edit == "ETY" | "None" | "?"  (id 0); else an edit.
.NOT_EDIT <- c("ETY", "None", "?")

# encode_cells: build the integer representation consumed by the C++ kernel.
#   - global integration order = union of names across cells (first-seen order)
#   - global edit-token table  = distinct edit tokens -> ids 1..K (0 = non-edit)
# Returns a list with code/depth/mask flat integer vectors and n,T,S,W,ids.
encode_cells <- function(cells) {
  ids <- names(cells)
  n <- length(cells)

  # T = union of integration names (a cell missing an integration => bit unset).
  tape_names <- unique(unlist(lapply(cells, names), use.names = FALSE))
  Tn <- length(tape_names)
  tape_index <- setNames(seq_len(Tn), tape_names)

  # S = tape length (sites per tape); derived from the first non-empty tape.
  S <- 6L
  for (c in cells) if (length(c)) { S <- length(c[[1L]]); break }

  W <- as.integer(ceiling(Tn / 32))

  # Global edit-token hash. factor() over all observed site strings gives stable
  # integer ids; we then remap non-edit tokens to id 0 and edits to 1..K.
  all_sites <- unlist(cells, use.names = FALSE)
  fac <- factor(all_sites)
  lev <- levels(fac)
  is_edit_lev <- !(lev %in% .NOT_EDIT)
  # ids: 0 for non-edit levels, else a dense 1..K id.
  edit_id <- integer(length(lev))
  edit_id[is_edit_lev] <- seq_len(sum(is_edit_lev))
  # token id per observed site (vector aligned with all_sites)
  site_id_all <- edit_id[as.integer(fac)]

  code  <- integer(n * Tn * S)          # [cell][tape][site] row-major
  depth <- integer(n * Tn)              # [cell][tape]
  mask  <- integer(n * W)               # packed uint32 words per cell

  pos <- 0L  # running index into site_id_all (cells unlisted in order)
  for (i in seq_len(n)) {
    cell <- cells[[i]]
    cn <- names(cell)
    base_code <- (i - 1L) * Tn * S
    base_depth <- (i - 1L) * Tn
    base_mask <- (i - 1L) * W
    for (k in seq_along(cell)) {
      tp <- cn[k]
      t <- tape_index[[tp]]            # 1-based global tape index
      sid <- site_id_all[(pos + 1L):(pos + S)]
      pos <- pos + S
      off <- base_code + (t - 1L) * S
      code[(off + 1L):(off + S)] <- sid
      d <- sum(sid > 0L)
      depth[base_depth + t] <- d
      # tape_recovered: present AND not all-missing (matches dtt_distance.R).
      # We are iterating tapes the cell actually has, so "present" holds; a tape
      # is recovered iff NOT every site is in MISSING_STATES {"None","?"}.
      recovered <- !all(cell[[k]] %in% c("None", "?"))
      if (recovered) {
        wi <- (t - 1L) %/% 32L
        bi <- (t - 1L) %%  32L
        mask[base_mask + wi + 1L] <- bitwOr(mask[base_mask + wi + 1L],
                                            bitwShiftL(1L, bi))
      }
    }
  }

  list(code = code, depth = depth, mask = mask,
       n = n, T = Tn, S = S, W = W, ids = ids)
}

# Run the kernel on a prebuilt encoding; returns the named numeric matrix.
dtt_distance_matrix_encoded <- function(enc) {
  m <- dtt_kernel(enc$code, enc$depth, enc$mask,
                  enc$n, enc$T, enc$S, enc$W)
  dimnames(m) <- list(enc$ids, enc$ids)
  m
}

# Drop-in for dtt_distance_matrix: encode then run.
dtt_distance_matrix_fast <- function(cells) {
  dtt_distance_matrix_encoded(encode_cells(cells))
}
