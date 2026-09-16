# Tree sensitivity
# Leave-one-tape-out sensitivity analysis: for cells sampled from each
# blastomere, rebuild trees with ape::nj() using only 10 of the cells' 11 DNA
# Typewriter tape integrations, then compare topology and the held-out
# tape's editing pattern against the ground-truth tree.
#
# Run from the repo root (paths below are relative to it).

library(ape)
library(phangorn)
library(TreeDist)

source("tree_building/1_build_nj_backbone/parse_tape_consensus.R")  # parse_cells()
source("tree_building/1_build_nj_backbone/dtt_distance.R")          # dtt_distance_matrix(), tape_distance()

out_dir <- "tree_building/results/tape_sensitivity"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

N_REPS  <- 5
N_CELLS <- 1000

#-------------------------------------------------------------------------#
# 1. Load the tree and per-blastomere tape data (parsed once, reused across
#    all replicates and held-out-tape iterations below).

tree_file <- "tree_building/results/3-dated-tree/merged_minB2h_lineage_constrained.nwk"
tree <- ape::read.tree(tree_file)
cat(sprintf("Tree: %d tips\n", length(tree$tip.label)))

cells_B1 <- parse_cells("tree_building/processed_data/e3v8.B1_tape_consensus.ge7_founderok.tsv.gz")
cells_B2 <- parse_cells("tree_building/processed_data/e3v8.B2_tape_consensus.ge7_founderok.tsv.gz")
cat(sprintf("Parsed %d B1 cells, %d B2 cells\n", length(cells_B1), length(cells_B2)))

stopifnot(all(tree$tip.label %in% c(names(cells_B1), names(cells_B2))))

# Blastomere identity is which tape-consensus file a cell_id appears in, NOT
# a structural tree split (same convention as tree_analysis/step2_sibling_cells.R
# and step3_clade_coincidence.R). Blastomere A = the larger of the two.
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

#-------------------------------------------------------------------------#
# 2. Sampling: random cells (5x per blastomere) and a single clade (5x per
# blastomere). Both strategies are run once per blastomere (A and B) below,
# so `scenario` (random/clade) and `blastomere` (A/B) are independent columns
# in the final summary rather than one label conflating the two.

# Scenario "random": n_cells random cells. Output per rep: sampled cell ids
# and the original tree topology restricted to those cells (keep.tip()).
sample_random_cells <- function(tips, blastomere_label) {
  cell_names   <- vector("list", N_REPS)
  true_subtree <- vector("list", N_REPS)
  tag <- sprintf("random_%s", blastomere_label)

  for (rep in seq_len(N_REPS)) {
    set.seed(rep)
    sampled_tips <- sample(tips, N_CELLS)

    cell_names[[rep]]   <- sampled_tips
    true_subtree[[rep]] <- ape::keep.tip(tree, sampled_tips)

    write.csv(data.frame(cell_id = sampled_tips),
              file.path(out_dir, sprintf("cell_names_%s_rep%d.csv", tag, rep)),
              row.names = FALSE)
    ape::write.tree(true_subtree[[rep]],
                    file.path(out_dir, sprintf("true_subtree_%s_rep%d.nwk", tag, rep)))
  }
  cat(sprintf("Random sampling (blastomere %s): wrote %d replicate cell-name lists + ground-truth subtrees to %s\n",
              blastomere_label, N_REPS, out_dir))
  list(cell_names = cell_names, true_subtree = true_subtree)
}

# Scenario "clade": a single real clade of ~N_CELLS tips (5x per blastomere).
# Adapted from tree_analysis/step3_clade_coincidence.R::maximal_clades() --
# same Descendants()/parent-lookup approach, but with a SIZE BAND instead of
# a ceiling: a node qualifies only if its own tip count is in [min,max] AND
# its parent's tip count exceeds max, i.e. it is maximal within the band.
# Sizes are non-decreasing from child to parent, so no qualifying node can be
# nested inside another -- the resulting clades' tip sets are automatically
# disjoint, no extra overlap-avoidance needed.

CLADE_MIN <- 900
CLADE_MAX <- 1100

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
    if (parent_of[u] != 0 && size[parent_of[u]] <= max_size) next  # not maximal
    candidates <- c(candidates, u)
  }
  list(nodes = candidates, sizes = size[candidates], desc = desc)
}

sample_clade_cells <- function(blastomere_tips, blastomere_label) {
  subtree <- ape::keep.tip(tree, blastomere_tips)
  band    <- find_band_clades(subtree, CLADE_MIN, CLADE_MAX)
  cat(sprintf("Clade sampling (blastomere %s): %d candidate clades with %d-%d tips\n",
              blastomere_label, length(band$nodes), CLADE_MIN, CLADE_MAX))
  stopifnot(length(band$nodes) >= N_REPS)

  set.seed(0)
  chosen_nodes <- sample(band$nodes, N_REPS)

  cell_names   <- vector("list", N_REPS)
  true_subtree <- vector("list", N_REPS)
  tag <- sprintf("clade_%s", blastomere_label)

  for (rep in seq_len(N_REPS)) {
    node         <- chosen_nodes[rep]
    sampled_tips <- subtree$tip.label[band$desc[[node]]]

    cell_names[[rep]]   <- sampled_tips
    true_subtree[[rep]] <- ape::keep.tip(tree, sampled_tips)

    write.csv(data.frame(cell_id = sampled_tips),
              file.path(out_dir, sprintf("cell_names_%s_rep%d.csv", tag, rep)),
              row.names = FALSE)
    ape::write.tree(true_subtree[[rep]],
                    file.path(out_dir, sprintf("true_subtree_%s_rep%d.nwk", tag, rep)))
  }
  cat(sprintf("Clade sampling (blastomere %s): wrote %d replicate cell-name lists + ground-truth subtrees to %s\n",
              blastomere_label, N_REPS, out_dir))
  list(cell_names = cell_names, true_subtree = true_subtree)
}

random_A <- sample_random_cells(blastomereA_tips, "A")
random_B <- sample_random_cells(blastomereB_tips, "B")
clade_A  <- sample_clade_cells(blastomereA_tips, "A")
clade_B  <- sample_clade_cells(blastomereB_tips, "B")

#-------------------------------------------------------------------------#
# 3. Leave-one-tape-out NJ tree: for each replicate (all 4 scenario x
# blastomere groups) and each of the cells' 11 tape integrations, drop that
# one integration and rebuild an NJ tree from the remaining 10 with
# ape::nj() -- the R-native NJ step, as opposed to decentTree's C++ RapidNJ
# used in production. The distance matrix itself still comes from
# dtt_distance_matrix() (dtt_distance.R) completely unmodified: "leave tape b
# out" is just removing element b from every cell's tape list before that
# (untouched) function is called.
#
# NOTE on runtime: this is pure-R dtt_distance_matrix() on ~1000 cells x 11
# held-out-tape variants x 20 replicates (4 groups x 5 reps) -- ~11M
# dtt_distance() calls total. If this proves too slow,
# tree_building/1_build_nj_backbone/dtt_distance_scale.R provides an
# Rcpp-accelerated drop-in replacement (dtt_distance_matrix_fast()) using the
# identical algorithm.

run_leave_one_tape_out <- function(sampled_tips, cells_source, label, rep) {
  cells_sub <- cells_source[sampled_tips]
  barcodes  <- names(cells_sub[[1]])
  cat(sprintf("%s rep %d: %d cells, %d tape integrations\n",
              label, rep, length(cells_sub), length(barcodes)))

  trees_est <- vector("list", length(barcodes))
  names(trees_est) <- barcodes

  for (b in barcodes) {
    cells_loo <- lapply(cells_sub, function(cell) cell[names(cell) != b])
    M <- dtt_distance_matrix(cells_loo)

    # QC (n_loci >= 7 of 11) guarantees every pair still shares several of
    # the remaining 10 integrations, so NA should not occur -- but impute
    # defensively (mean observed distance) rather than let ape::nj() error.
    n_na <- sum(is.na(M))
    if (n_na > 0) {
      warning(sprintf("%s rep %d, held-out tape %s: %d/%d distance pairs are NA -- imputing with mean observed distance",
                      label, rep, b, n_na, length(M)))
      M[is.na(M)] <- mean(M, na.rm = TRUE)
    }

    tree_est       <- ape::nj(M)
    trees_est[[b]] <- tree_est

    ape::write.tree(tree_est,
                    file.path(out_dir, sprintf("est_tree_%s_rep%d_heldout-%s.nwk", label, rep, b)))
  }
  trees_est
}

build_trees_for_group <- function(group, cells_source, label) {
  lapply(seq_len(N_REPS), function(rep) {
    run_leave_one_tape_out(group$cell_names[[rep]], cells_source, label, rep)
  })
}

trees_est_random_A <- build_trees_for_group(random_A, cells_A, "random_A")
trees_est_random_B <- build_trees_for_group(random_B, cells_B, "random_B")
trees_est_clade_A  <- build_trees_for_group(clade_A,  cells_A, "clade_A")
trees_est_clade_B  <- build_trees_for_group(clade_B,  cells_B, "clade_B")

#-------------------------------------------------------------------------#
# 4. Topology comparison + held-out-tape concordance.
#
# For each (scenario, blastomere, replicate, held-out tape b):
#   a) topology: TreeDist::PhylogeneticInfoDistance(), normalized -- an
#      information-theoretic generalization of Robinson-Foulds distance
#      (Smith 2020) between the ground-truth subtree (keep.tip() on the full
#      tree, step 2) and the NJ tree built without tape b (step 3). Unlike
#      phangorn::RF.dist(), it is well-defined and correctly normalized even
#      when one tree has polytomies -- the ground-truth subtree, inherited
#      from the dated/lineage-constrained tree, can carry real
#      multifurcations that ape::nj()'s always-fully-resolved output does
#      not; plain RF.dist's normalize=TRUE assumes both trees are fully
#      bifurcating and warns "Some trees are not binary" in that case.
#   b) does tape b's OWN editing pattern (never used to build the tree in
#      step 3) still agree with that tree's topology? Two complementary
#      checks, EACH benchmarked against a permutation null rather than the
#      true tree, so both report the same kind of comparison:
#      - parsimony: phangorn::phyDat() + phangorn::parsimony() score tape b's
#        6 sites on the held-out-tape tree, vs. a tip-label permutation null
#        (same tree, shuffled tip<->data mapping) -- reported as a z-score.
#      - distance correlation: dtt_distance.R's own tape_distance() gives
#        tape b's pairwise distance directly; Spearman cor.test() (the same
#        idiom already used in tree_analysis/step2_sibling_cells.R and
#        step3_clade_coincidence.R) correlates it against the held-out-tape
#        tree's own patristic distance ("leave-one-out tree"), benchmarked
#        against the SAME tree with cell<->tip identity shuffled ("random
#        tree" null).

N_PERM    <- 100   # permutation-null reps (parsimony z-score, random-tree null)
MAX_PAIRS <- 20000 # cap on cell pairs for the distance-correlation check

compute_replicate_metrics <- function(sampled_tips, true_tree, trees_est, cells_source,
                                       scenario_label, blastomere_label, rep) {
  cells_sub <- cells_source[sampled_tips]
  barcodes  <- names(trees_est)

  rows <- lapply(barcodes, function(b) {
    tree_est <- trees_est[[b]]

    # ---- a) topology: estimated tree vs. ground truth ----
    phylo_info_dist <- TreeDist::PhylogeneticInfoDistance(true_tree, tree_est, normalize = TRUE)

    # ---- b-i) parsimony of the held-out tape's 6 sites on the tree ----
    site_mat <- t(vapply(cells_sub, function(cell) cell[[b]], character(6)))
    rownames(site_mat) <- names(cells_sub)
    # phyDat() treats R-level NA as an "unknown character" and DELETES the
    # whole site (column) if any taxon has one there, rather than treating
    # just that taxon as ambiguous at that site -- so missing states must be
    # recoded as an explicit token matched via `ambiguity`, not NA (this was
    # silently zeroing out parsimony scores for almost every tape/site).
    site_mat[site_mat %in% MISSING_STATES] <- "?"
    levels_b <- sort(setdiff(unique(as.vector(site_mat)), "?"))
    # Some integrations are fixed at a single state across an entire
    # blastomere (e.g. an edit that occurred once, very early, before the
    # blastomere split) -- invariant within this sample, so they carry no
    # lineage information here: parsimony is trivially 0 on ANY tree, and
    # tape_distance() below is constant 0 for every pair (rho undefined, not
    # a bug). Flag it explicitly rather than let those rows read as NA/0
    # with no explanation.
    tape_invariant <- length(levels_b) <= 1
    pd <- phangorn::phyDat(site_mat, type = "USER", levels = levels_b, ambiguity = "?")

    p_est  <- phangorn::parsimony(tree_est, pd)
    p_true <- phangorn::parsimony(true_tree, pd)

    p_null <- vapply(seq_len(N_PERM), function(s) {
      set.seed(s)
      pd_perm <- pd
      names(pd_perm) <- sample(names(pd))   # shuffle tip<->data mapping, same tree
      phangorn::parsimony(tree_est, pd_perm)
    }, numeric(1))

    # ---- b-ii) single-tape distance vs. leave-one-out tree patristic
    # distance, benchmarked against a "random tree" null: the SAME tree_est
    # (identical shape/branch lengths), but with which cell sits at which
    # tip shuffled, so the tape-state<->tree-position association is
    # randomized while the tree itself is held fixed.
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
          perm       <- setNames(sample(recovered), recovered)
          patr_null  <- cophenetic_est[cbind(perm[pairs[, 1]], perm[pairs[, 2]])]
          unname(suppressWarnings(cor.test(tape_d, patr_null, method = "spearman")$estimate))
        }, numeric(1))
        rho_null_mean <- mean(rho_null)
        rho_null_sd   <- sd(rho_null)
      }
    }

    data.frame(
      scenario = scenario_label, blastomere = blastomere_label, rep = rep, heldout_tape = b,
      n_cells = length(sampled_tips), n_recovered = length(recovered),
      tape_invariant = tape_invariant,
      phylo_info_dist = phylo_info_dist,
      parsimony_est = p_est, parsimony_true = p_true,
      parsimony_null_mean = mean(p_null), parsimony_null_sd = sd(p_null),
      parsimony_z = (p_est - mean(p_null)) / sd(p_null),
      rho_leave_one_out = rho_est,
      rho_random_mean = rho_null_mean, rho_random_sd = rho_null_sd,
      rho_z = (rho_est - rho_null_mean) / rho_null_sd
    )
  })

  do.call(rbind, rows)
}

compute_group_metrics <- function(group, trees_est_group, cells_source, scenario_label, blastomere_label) {
  do.call(rbind, lapply(seq_len(N_REPS), function(rep) {
    compute_replicate_metrics(group$cell_names[[rep]], group$true_subtree[[rep]], trees_est_group[[rep]],
                               cells_source, scenario_label, blastomere_label, rep)
  }))
}

metrics_random_A <- compute_group_metrics(random_A, trees_est_random_A, cells_A, "random", "A")
metrics_random_B <- compute_group_metrics(random_B, trees_est_random_B, cells_B, "random", "B")
metrics_clade_A  <- compute_group_metrics(clade_A,  trees_est_clade_A,  cells_A, "clade",  "A")
metrics_clade_B  <- compute_group_metrics(clade_B,  trees_est_clade_B,  cells_B, "clade",  "B")

summary_df <- rbind(metrics_random_A, metrics_random_B, metrics_clade_A, metrics_clade_B)
write.csv(summary_df, file.path(out_dir, "summary.csv"), row.names = FALSE)
cat(sprintf("Wrote summary with %d rows to %s\n", nrow(summary_df), file.path(out_dir, "summary.csv")))

# Tapes invariant within the sampled cells (see tape_invariant above) always
# read as parsimony 0 / rho NA -- not a bug, but worth seeing at a glance
# which barcode(s) those are and how often.
invariant_by_tape <- table(summary_df$heldout_tape[summary_df$tape_invariant])
if (length(invariant_by_tape) > 0) {
  cat("Invariant (uninformative) held-out tape, by barcode:\n")
  print(invariant_by_tape)
}

# Is a negative leave-one-out rho concentrated in a particular held-out tape,
# or spread evenly? (heldout_tape is already a summary_df column.)
neg_rho_by_tape <- table(summary_df$heldout_tape[summary_df$rho_leave_one_out < 0])
if (length(neg_rho_by_tape) > 0) {
  cat("Replicates with negative rho_leave_one_out, by held-out tape:\n")
  print(neg_rho_by_tape)
}

#-------------------------------------------------------------------------#
# 5. Summary plots. One point per (replicate, held-out tape); conventions
# (fig_theme, out_dir, ggsave sizing) follow making_figures/single_blastomere_replicates.R.

library(ggplot2)

fig_out_dir <- "figures/tape_sensitivity"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)

base_size <- 10
fig_theme <- theme_classic(base_size = base_size, base_family = "sans") +
  theme(
    axis.title  = element_text(size = base_size),
    axis.text   = element_text(size = base_size - 1, color = "black"),
    axis.line   = element_line(linewidth = 0.3, color = "black"),
    axis.ticks  = element_line(linewidth = 0.3, color = "black"),
    legend.text  = element_text(size = base_size - 1),
    legend.title = element_text(size = base_size - 1),
    plot.margin  = margin(3, 4, 3, 3)
  )

scenario_labels <- c(random = "Random cells", clade = "Clade")
summary_df$scenario_label <- factor(scenario_labels[summary_df$scenario], levels = scenario_labels)
summary_df$blastomere     <- factor(summary_df$blastomere, levels = c("A", "B"))
blastomere_labeller <- labeller(blastomere = c(A = "Blastomere A", B = "Blastomere B"))

# Panel 1: topology recovery -- normalized phylogenetic information distance
p_topology <- ggplot(summary_df, aes(x = scenario_label, y = phylo_info_dist)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, fill = "grey90") +
  geom_jitter(aes(color = factor(rep)), width = 0.15, size = 1.5, alpha = 0.8) +
  facet_wrap(~blastomere, labeller = blastomere_labeller) +
  labs(x = NULL, y = "Normalized phylogenetic info. distance\n(leave-one-tape-out NJ tree vs. ground truth)", color = "replicate") +
  fig_theme

p_topology
write.csv(summary_df[, c("scenario", "blastomere", "rep", "heldout_tape", "phylo_info_dist")],
          file.path(fig_out_dir, "panel_topology_distance.csv"), row.names = FALSE)
ggsave(file.path(fig_out_dir, "panel_topology_distance.pdf"), p_topology, width = 5.5, height = 3.5, dpi = 300)

# Panel 2: does the held-out tape's own editing pattern still cluster on the
# tree that never saw it?
# p_parsimony <- ggplot(summary_df, aes(x = scenario_label, y = parsimony_z)) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
#   geom_boxplot(outlier.shape = NA, width = 0.5, fill = "grey90") +
#   geom_jitter(aes(color = factor(rep)), width = 0.15, size = 1.5, alpha = 0.8) +
#   facet_wrap(~blastomere, labeller = blastomere_labeller) +
#   labs(x = NULL, y = "Held-out tape parsimony z-score\n(vs. tip-label permutation null)", color = "replicate") +
#   fig_theme
#
# p_parsimony
# write.csv(summary_df[, c("scenario", "blastomere", "rep", "heldout_tape", "parsimony_est", "parsimony_true",
#                          "parsimony_null_mean", "parsimony_null_sd", "parsimony_z")],
#           file.path(fig_out_dir, "panel_heldout_tape_parsimony.csv"), row.names = FALSE)
# ggsave(file.path(fig_out_dir, "panel_heldout_tape_parsimony.pdf"), p_parsimony, width = 5.5, height = 3.5, dpi = 300)

# Panel 3: does the held-out tape's own pairwise distance correlate with
# patristic distance on the leave-one-out tree more than on a random tree?
# "Random tree" = the SAME tree (identical shape/branch lengths) with
# cell<->tip identity shuffled -- one point per (scenario, blastomere, rep,
# held-out tape) x category, so the two categories are directly comparable.
rho_long <- rbind(
  data.frame(scenario_label = summary_df$scenario_label, blastomere = summary_df$blastomere, rep = summary_df$rep,
             heldout_tape = summary_df$heldout_tape,
             category = "Leave-one-out tree", rho = summary_df$rho_leave_one_out),
  data.frame(scenario_label = summary_df$scenario_label, blastomere = summary_df$blastomere, rep = summary_df$rep,
             heldout_tape = summary_df$heldout_tape,
             category = "Random tree", rho = summary_df$rho_random_mean)
)
rho_long$category <- factor(rho_long$category, levels = c("Leave-one-out tree", "Random tree"))

cols = RColorBrewer::brewer.pal(n = 11, name = "Set3")
p_rho <- ggplot(rho_long, aes(x = category, y = rho)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, fill = "grey90") +
  geom_jitter(aes(color = factor(heldout_tape)), width = 0.15, size = 1.5, alpha = 0.8) +
  scale_color_manual(values = cols)+
  facet_grid(blastomere ~ scenario_label, labeller = blastomere_labeller) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  labs(x = NULL, y = "Spearman rho between \nheld-out tape distance and tree patristic distance", color = "Left-out-tape") +
  fig_theme

p_rho
write.csv(summary_df[, c("scenario", "blastomere", "rep", "heldout_tape", "n_recovered",
                         "rho_leave_one_out", "rho_random_mean", "rho_random_sd", "rho_z")],
          file.path(fig_out_dir, "panel_heldout_tape_distance_correlation.csv"), row.names = FALSE)
ggsave(file.path(fig_out_dir, "panel_heldout_tape_distance_correlation.pdf"), p_rho, width = 7, height = 6, dpi = 300)

cat(sprintf("Wrote 2 summary plots (+ underlying CSVs) to %s\n", fig_out_dir))
