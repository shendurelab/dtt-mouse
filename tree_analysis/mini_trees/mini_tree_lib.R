# mini_tree_lib.R
# Build a per-integration ("mini") lineage tree from bulk TAPE genotypes.
#
# Input: rows of tables/DTTz_3_S3.bulk_lineage_genotypes.csv for ONE tapebc.
# Each row is one distinct GENOTYPE observed at that integration -- a length-6
# vector of site states (six '|'-joined tokens, 'U' = unedited) -- with
# EDIT_DEPTH = number of edited (non-'U') sites on that genotype's tape.
#
# Distance: reuses dtt_distance.R's tape_distance() unchanged (the same
# per-integration term used in the whole-cell distance), so the metric here is
# provably identical to the one already tested in src/1-dist-mat. That file
# uses "ETY" for unedited/"None","?" for missing; this table only ever has
# 'U' (never missing), so we translate 'U' -> "ETY" before calling it.

source("tree_building/1_build_nj_backbone/dtt_distance.R")   # tape_depth, shared_edit_prefix, tape_distance

# Genotype string ("A|B|U|U|U|U") -> raw length-6 character vector, 'U' kept
# as-is. For display/plotting; parse_genotype() below builds on this for the
# distance calculation's own recoding.
split_genotype_sites <- function(genotype_str) strsplit(genotype_str, "|", fixed = TRUE)[[1]]

# Genotype string -> length-6 character vector, with 'U' recoded to "ETY" so
# it can be passed straight into dtt_distance.R's
# is_edit()/tape_depth()/shared_edit_prefix().
parse_genotype <- function(genotype_str) {
  sites <- split_genotype_sites(genotype_str)
  ifelse(sites == "U", "ETY", sites)
}

# n x n distance matrix (as a dist object) over the genotypes in one tapebc's
# rows of the long table, named by genotype_id. Pairwise tape_distance() over
# every pair -- fine at this scale (<=163 genotypes per integration).
genotype_distance_matrix <- function(genotypes_df) {
  ids <- genotypes_df$genotype_id
  n <- length(ids)
  parsed <- lapply(genotypes_df$genotype, parse_genotype)
  m <- matrix(0, n, n, dimnames = list(ids, ids))
  for (i in seq_len(n)) for (j in seq_len(n)) {
    m[i, j] <- tape_distance(parsed[[i]], parsed[[j]])
  }
  as.dist(m)
}

# Neighbour-Joining + midpoint rooting over one tapebc's genotypes. Returns a
# rooted `phylo` object with tip labels = genotype_id.
build_mini_tree <- function(genotypes_df) {
  d <- genotype_distance_matrix(genotypes_df)
  tree <- ape::nj(d)
  phangorn::midpoint(tree)
}

# Genotype table -> long df (label, site = factor "1".."6", symbol) for the
# tip-aligned tape/alignment heatmap, joined onto tip_y's y-coordinate (e.g.
# `ggtree(tree)$data`'s isTip rows) so each genotype's tile row lines up with
# its own tip. Raw symbols ('U' for unedited) -- not parse_genotype()'s "ETY"
# recoding, which is for the distance calc only.
build_tape_long_df <- function(genotypes_df, tip_y) {
  site_mat <- t(vapply(genotypes_df$genotype, split_genotype_sites, character(6)))
  rownames(site_mat) <- genotypes_df$genotype_id
  colnames(site_mat) <- paste0("site", 1:6)
  long <- as.data.frame(site_mat)
  long$label <- rownames(long)
  long <- tidyr::pivot_longer(long, cols = tidyr::starts_with("site"),
                               names_to = "site", values_to = "symbol")
  long$site <- factor(long$site, levels = paste0("site", 1:6), labels = 1:6)
  merge(long, tip_y, by = "label")
}
