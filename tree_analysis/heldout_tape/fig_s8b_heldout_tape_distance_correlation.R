#!/usr/bin/env Rscript
# fig_s8b_heldout_tape_distance_correlation.R
# Fig. S8B: for cells sampled from each blastomere (random cells, and a single
# real clade; 5 replicates each), rebuild an NJ tree leaving out each of the
# 11 tape integrations in turn, then ask: does the held-out tape's OWN pairwise
# editing distance still correlate with patristic distance on the tree that
# never saw it ("leave-one-out tree"), more than on the same tree with cell
# identity shuffled ("random tree" null)?
#
#
#
# Run from the repo root:
#   Rscript tree_analysis/heldout_tape/fig_s8b_heldout_tape_distance_correlation.R

suppressPackageStartupMessages({ library(ape) })
source("tree_building/1_build_nj_backbone/parse_tape_consensus.R")  # parse_cells()
source("tree_building/1_build_nj_backbone/dtt_distance.R")          # tape_distance(), tape_recovered(), MISSING_STATES
source("tree_building/1_build_nj_backbone/dtt_distance_scale.R")    # dtt_distance_matrix_fast() -- Rcpp, drop-in for dtt_distance_matrix()

TREE_FILE <- "tree_building/results/3-dated-tree/merged_minB2h_lineage_constrained.nwk"
B1_TSV    <- "tree_building/processed_data/e3v8.B1_tape_consensus.ge7_founderok.tsv.gz"
B2_TSV    <- "tree_building/processed_data/e3v8.B2_tape_consensus.ge7_founderok.tsv.gz"
stopifnot(file.exists(TREE_FILE), file.exists(B1_TSV), file.exists(B2_TSV))

OUT_DIR <- "tree_analysis/heldout_tape/results"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_CSV <- "figures_data/fig_s8_heldout_tape_distance_correlation.csv"

N_REPS  <- 5
N_CELLS <- 1000

# ---- 1. tree + per-blastomere tape data ------------------------------------
tree <- ape::read.tree(TREE_FILE)
cat(sprintf("Tree: %d tips\n", length(tree$tip.label)))

cells_B1 <- parse_cells(B1_TSV)
cells_B2 <- parse_cells(B2_TSV)
cat(sprintf("Parsed %d B1 cells, %d B2 cells\n", length(cells_B1), length(cells_B2)))
stopifnot(all(tree$tip.label %in% c(names(cells_B1), names(cells_B2))))

# Blastomere identity is which tape-consensus file a cell_id appears in.
# Blastomere A = the larger of the two.
tree_tips_B1 <- tree$tip.label[tree$tip.label %in% names(cells_B1)]
tree_tips_B2 <- tree$tip.label[tree$tip.label %in% names(cells_B2)]
stopifnot(length(intersect(tree_tips_B1, tree_tips_B2)) == 0)

if (length(tree_tips_B1) >= length(tree_tips_B2)) {
  blastomereA_tips <- tree_tips_B1; cells_A <- cells_B1
  blastomereB_tips <- tree_tips_B2; cells_B <- cells_B2
} else {
  blastomereA_tips <- tree_tips_B2; cells_A <- cells_B2
  blastomereB_tips <- tree_tips_B1; cells_B <- cells_B1
}
cat(sprintf("Blastomere A: %d tips; Blastomere B: %d tips\n",
            length(blastomereA_tips), length(blastomereB_tips)))

# ---- 2. sampling: random cells and a single clade, 5 reps per blastomere ---
sample_random_cells <- function(tips) {
  lapply(seq_len(N_REPS), function(rep) {
    set.seed(rep)
    sample(tips, N_CELLS)
  })
}

CLADE_MIN <- 900
CLADE_MAX <- 1100

# Adapted from tree_analysis/step3_clade_coincidence.R::maximal_clades(): a
# node qualifies only if its own tip count is in [min,max] AND its parent's
# tip count exceeds max, i.e. it is maximal within the band. Sizes are
# non-decreasing child->parent, so qualifying nodes' tip sets are disjoint.
find_band_clades <- function(subtree, min_size, max_size) {
  n_tips  <- Ntip(subtree)
  n_nodes <- n_tips + Nnode(subtree)
  desc <- Descendants(subtree, 1:n_nodes, type = "tips")
  size <- lengths(desc)
  parent_of <- integer(n_nodes)
  parent_of[subtree$edge[, 2]] <- subtree$edge[, 1]
  candidates <- c()
  for (u in (n_tips + 1):n_nodes) {
    if (size[u] < min_size || size[u] > max_size) next
    if (parent_of[u] != 0 && size[parent_of[u]] <= max_size) next
    candidates <- c(candidates, u)
  }
  list(nodes = candidates, desc = desc)
}

sample_clade_cells <- function(blastomere_tips) {
  suppressPackageStartupMessages(library(phangorn))
  subtree <- ape::keep.tip(tree, blastomere_tips)
  band <- find_band_clades(subtree, CLADE_MIN, CLADE_MAX)
  stopifnot(length(band$nodes) >= N_REPS)
  set.seed(0)
  chosen_nodes <- sample(band$nodes, N_REPS)
  lapply(chosen_nodes, function(node) subtree$tip.label[band$desc[[node]]])
}

random_A <- sample_random_cells(blastomereA_tips)
random_B <- sample_random_cells(blastomereB_tips)
clade_A  <- sample_clade_cells(blastomereA_tips)
clade_B  <- sample_clade_cells(blastomereB_tips)

# ---- 3. leave-one-tape-out NJ trees ----------------------------------------
run_leave_one_tape_out <- function(sampled_tips, cells_source, label, rep) {
  cells_sub <- cells_source[sampled_tips]
  barcodes  <- names(cells_sub[[1]])
  cat(sprintf("%s rep %d: %d cells, %d tape integrations\n",
              label, rep, length(cells_sub), length(barcodes)))

  trees_est <- vector("list", length(barcodes)); names(trees_est) <- barcodes
  for (b in barcodes) {
    cells_loo <- lapply(cells_sub, function(cell) cell[names(cell) != b])
    M <- dtt_distance_matrix_fast(cells_loo)
    n_na <- sum(is.na(M))
    if (n_na > 0) {
      warning(sprintf("%s rep %d, held-out tape %s: %d/%d distance pairs are NA -- imputing with mean observed distance",
                      label, rep, b, n_na, length(M)))
      M[is.na(M)] <- mean(M, na.rm = TRUE)
    }
    tree_est <- ape::nj(M)
    trees_est[[b]] <- tree_est
    ape::write.tree(tree_est, file.path(OUT_DIR, sprintf("est_tree_%s_rep%d_heldout-%s.nwk", label, rep, b)))
  }
  trees_est
}

build_trees_for_group <- function(group, cells_source, label) {
  lapply(seq_len(N_REPS), function(rep) run_leave_one_tape_out(group[[rep]], cells_source, label, rep))
}

trees_est_random_A <- build_trees_for_group(random_A, cells_A, "random_A")
trees_est_random_B <- build_trees_for_group(random_B, cells_B, "random_B")
trees_est_clade_A  <- build_trees_for_group(clade_A,  cells_A, "clade_A")
trees_est_clade_B  <- build_trees_for_group(clade_B,  cells_B, "clade_B")

# ---- 4. held-out tape distance vs. tree patristic distance -----------------
# rho_est: Spearman rho between the held-out tape's own pairwise distance and
# patristic distance on the tree built WITHOUT it. rho_null: the same tree
# (identical shape/branch lengths) with cell<->tip identity shuffled, so the
# tape-state<->tree-position association is randomized while the tree itself
# is held fixed -- benchmarks rho_est against chance structure in the tree.
N_PERM    <- 100
MAX_PAIRS <- 20000

compute_replicate_rho <- function(sampled_tips, trees_est, cells_source, scenario_label, blastomere_label, rep) {
  cells_sub <- cells_source[sampled_tips]
  barcodes  <- names(trees_est)

  rows <- lapply(barcodes, function(b) {
    tree_est <- trees_est[[b]]
    recovered <- names(cells_sub)[vapply(cells_sub, function(cell) tape_recovered(cell[[b]]), logical(1))]
    rho_est <- NA_real_; rho_null_mean <- NA_real_; rho_null_sd <- NA_real_

    if (length(recovered) >= 10) {
      pairs <- t(combn(recovered, 2))
      n_pairs_total <- nrow(pairs)
      if (n_pairs_total > MAX_PAIRS) {
        set.seed(1)
        pairs <- pairs[sample(n_pairs_total, MAX_PAIRS), , drop = FALSE]
        cat(sprintf("  %s_%s rep %d, tape %s: capping distance-correlation pairs to %d (of %d)\n",
                    scenario_label, blastomere_label, rep, b, MAX_PAIRS, n_pairs_total))
      }
      tape_d <- apply(pairs, 1, function(p) tape_distance(cells_sub[[p[1]]][[b]], cells_sub[[p[2]]][[b]]))

      if (sd(tape_d) == 0) {
        cat(sprintf("  %s_%s rep %d, tape %s: tape_distance is constant (0 variance) across %d recovered cells -- rho left NA\n",
                    scenario_label, blastomere_label, rep, b, length(recovered)))
      } else {
        cophenetic_est <- ape::cophenetic.phylo(tree_est)
        patr_est <- cophenetic_est[cbind(pairs[, 1], pairs[, 2])]
        rho_est  <- unname(suppressWarnings(cor.test(tape_d, patr_est, method = "spearman")$estimate))

        rho_null <- vapply(seq_len(N_PERM), function(s) {
          set.seed(s)
          perm      <- setNames(sample(recovered), recovered)
          patr_null <- cophenetic_est[cbind(perm[pairs[, 1]], perm[pairs[, 2]])]
          unname(suppressWarnings(cor.test(tape_d, patr_null, method = "spearman")$estimate))
        }, numeric(1))
        rho_null_mean <- mean(rho_null)
        rho_null_sd   <- sd(rho_null)
      }
    }

    data.frame(
      scenario = scenario_label, blastomere = blastomere_label, rep = rep, heldout_tape = b,
      n_recovered = length(recovered),
      rho_leave_one_out = rho_est,
      rho_random_mean = rho_null_mean, rho_random_sd = rho_null_sd,
      rho_z = (rho_est - rho_null_mean) / rho_null_sd
    )
  })
  do.call(rbind, rows)
}

compute_group_rho <- function(group, trees_est_group, cells_source, scenario_label, blastomere_label) {
  do.call(rbind, lapply(seq_len(N_REPS), function(rep)
    compute_replicate_rho(group[[rep]], trees_est_group[[rep]], cells_source, scenario_label, blastomere_label, rep)))
}

summary_df <- rbind(
  compute_group_rho(random_A, trees_est_random_A, cells_A, "random", "A"),
  compute_group_rho(random_B, trees_est_random_B, cells_B, "random", "B"),
  compute_group_rho(clade_A,  trees_est_clade_A,  cells_A, "clade",  "A"),
  compute_group_rho(clade_B,  trees_est_clade_B,  cells_B, "clade",  "B")
)

write.csv(summary_df, OUT_CSV, row.names = FALSE)
cat(sprintf("wrote %s (%d rows)\n", OUT_CSV, nrow(summary_df)))
