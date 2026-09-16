# make_synthetic_root_outgroup.R
# Synthetic outgroup cell constructors for run_full_distance.R's OUTGROUP_MODE
# branching. Two unrelated constructions live here because both are "build one
# fake cell to append before the NJ build," not because they share logic:
#
#   make_synthetic_root_cell()  (OUTGROUP_MODE=synthroot, the one v6 uses)
#   make_alt_edited_cell()      (OUTGROUP_MODE=alt/both -- not used by v6)
#
# ============================================================================
# make_synthetic_root_cell -- blastomere-founder outgroup (SYNTHETIC_ROOT_B1 /
# SYNTHETIC_ROOT_B2), built from a defining-sites table (this repo's v6 data
# ships processed_data/e3v5v6.defining_sites.tsv), in the dtt cell
# representation parse_tape_consensus.R uses (named list of 11 tapes, each a
# length-6 character vector; "ETY" = unedited).
#
# Mirrors an earlier prototype's make_synthetic_root() (not part of this
# extract) -- same "no signal at this tape -> leave it fully unedited, don't
# claim what isn't supported" principle -- but keyed on defining_sites.tsv's
# per-side site-1 allele instead of a >=90%-consensus supporting chain, and
# only ever sets site 1 (deeper sites carry no additional discriminating
# power for the blastomere split).
#
# defining_sites.tsv has 20 rows, not 22 (11 barcodes x 2 sides): two
# barcodes (TGACTAAAGCGG, CAGCTAACGCCT) have no reported B2 allele, so
# make_synthetic_root_cell() leaves those tapes fully unedited for B2 -- this
# is intentional, not a missing-data bug (verified against the source table).
#
# make_synthetic_root_cell(barcodes, defining_sites, side)
#   barcodes       : the integration ids to give the cell (use names() of a
#                    real cell so it shares all integrations with every
#                    other cell -- same convention as make_unedited_cell())
#   defining_sites : data.frame from e3.defining_sites.tsv (columns
#                    integration, site, blastomere, allele, ...)
#   side           : "B1" or "B2"
make_synthetic_root_cell <- function(barcodes, defining_sites, side) {
  stopifnot(side %in% c("B1", "B2"))
  setNames(lapply(barcodes, function(bc) {
    hit <- defining_sites[defining_sites$integration == bc &
                             defining_sites$blastomere == side &
                             defining_sites$site == "site1", ]
    sites <- rep("ETY", 6L)
    if (nrow(hit) == 1) sites[1] <- hit$allele[[1]]
    sites
  }), barcodes)
}

# ============================================================================
# make_alt_edited_cell -- a synthetic "alternatively edited" outgroup cell
#
# Run (example):
#   source("1_build_nj_backbone/parse_tape_consensus.R")   # parse_cells
#   source("1_build_nj_backbone/dtt_distance.R")           # is_edit, tape_depth, tape_recovered
#   source("1_build_nj_backbone/make_synthetic_root_outgroup.R")
#   cells      <- parse_cells(consensus_tsv_path("."))
#   depth_pool <- unlist(lapply(cells, function(cell) {
#     recovered <- vapply(cell, tape_recovered, logical(1))
#     vapply(cell[recovered], tape_depth, integer(1))
#   }))
#   alt_cell <- make_alt_edited_cell(cells, names(cells[[1]]), depth_pool, seed = 1)

# Number of leading tape positions that must get a never-observed-at-this-
# column edit (see the comment block above for why this guarantees full
# divergence).
N_NOVEL_SITES <- 2L

# build_column_vocab: for every barcode (tape/integration), the set of edit
# tokens actually observed at each of its 6 site columns across `cells`, plus
# the dataset-wide set of edit tokens observed anywhere (any tape, any site).
#   cells    : named list of cells (parse_cells() representation)
#   barcodes : integration ids to build vocab for (use names(cells[[1]]))
# Returns list(observed = <barcode -> site(1..6) -> character vector>,
#              global_vocab = <character vector>).
build_column_vocab <- function(cells, barcodes) {
  # 1. dataset-wide edit vocabulary: every site value, any tape, any cell,
  #    that is actually an edit (mirrors dtt_distance_scale.R's encode_cells()
  #    global edit-token construction).
  all_sites <- unlist(cells, use.names = FALSE)
  global_vocab <- unique(all_sites[is_edit(all_sites)])

  # 2. per-(barcode, site) observed edit vocabulary: what real cells actually
  #    carry at that exact tape+position, so we can find the complement.
  observed <- setNames(vector("list", length(barcodes)), barcodes)
  for (b in barcodes) {
    col_vals <- lapply(seq_len(6L), function(s) {
      v <- vapply(cells, function(cell) cell[[b]][s], character(1))
      unique(v[is_edit(v)])
    })
    observed[[b]] <- setNames(col_vals, as.character(seq_len(6L)))
  }

  list(observed = observed, global_vocab = global_vocab)
}

# sample_alt_tape: one synthetic tape for `barcode` at the target `depth`
# (0..6 edited sites).
#   - sites 1..min(depth, N_NOVEL_SITES): drawn from the dataset-wide edit
#     vocabulary MINUS whatever real cells carry at that exact column (never
#     observed there => guaranteed different from every real cell's tape).
#   - sites (N_NOVEL_SITES+1)..depth (if depth > N_NOVEL_SITES): drawn from
#     the REAL observed vocabulary at that column -- novelty is already
#     guaranteed by the leading sites, so these can look like ordinary edits
#     (falls back to the global vocabulary if a column happens to have no
#     observed edits, e.g. a rarely-edited terminal site).
#   - sites depth+1..6: "ETY" (unedited placeholder), matching the sequential
#     editing model and make_unedited_cell()'s convention.
sample_alt_tape <- function(vocab, barcode, depth) {
  sites <- rep(UNEDITED_STATE, 6L)
  if (depth < 1L) return(sites)

  n_novel <- min(depth, N_NOVEL_SITES)
  for (s in seq_len(n_novel)) {
    pool <- setdiff(vocab$global_vocab, vocab$observed[[barcode]][[as.character(s)]])
    if (length(pool) == 0L) {
      stop(sprintf("no never-observed edit tokens available for barcode '%s' site %d",
                    barcode, s))
    }
    sites[s] <- sample(pool, 1L)
  }

  if (depth > N_NOVEL_SITES) {
    for (s in (N_NOVEL_SITES + 1L):depth) {
      pool <- vocab$observed[[barcode]][[as.character(s)]]
      if (length(pool) == 0L) pool <- vocab$global_vocab
      sites[s] <- sample(pool, 1L)
    }
  }

  sites
}

# make_alt_edited_cell: the full 11-tape synthetic "alternatively edited"
# outgroup cell.
#   cells      : named list of real cells (used only to derive the vocab)
#   barcodes   : integration ids to give the cell (use names(cells[[1]]))
#   depth_pool : empirical (cell, tape) edit-depth vector to sample per-tape
#                depths from
#   seed       : RNG seed, set once up front so the whole cell is reproducible
make_alt_edited_cell <- function(cells, barcodes, depth_pool, seed = 1) {
  vocab <- build_column_vocab(cells, barcodes)
  set.seed(seed)
  depths <- sample(depth_pool, length(barcodes), replace = TRUE)
  tapes <- lapply(seq_along(barcodes), function(i) {
    sample_alt_tape(vocab, barcodes[i], depths[i])
  })
  setNames(tapes, barcodes)
}
