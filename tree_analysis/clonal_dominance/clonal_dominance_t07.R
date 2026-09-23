library(castor)
library(ggplot2)
library(patchwork)
library(scales)
library(RColorBrewer) 


#-------------------------------------------------------------------------#
# Load tree
n_lineages_to_find <- 200
T0 = 7.0

#backbone tree
tree_file = "tree_building/results/3-dated-tree/merged_minB2h_lineage_constrained.nwk"
bt = ape::read.tree(tree_file)

# full tree
FULL_TREE_PATH <- "tree_building/results/4-full-tree/merged_full_placed.nwk"
ft_tree <- ape::read.tree(FULL_TREE_PATH)

out_dir = "figures_data/"
out_dir_backbone = "figures_data/"
#-------------------------------------------------------------------------#
# Analyse static clonal dominance on full tree, panels a-c
ft_split      <- castor::split_tree_at_height(ft_tree, height = T0)
ft_clone_size <- vapply(ft_split$subtrees, function(s) length(s$tree$tip.label), integer(1))
ft_N0         <- length(ft_clone_size)

ft_total_tips <- length(ft_tree$tip.label)
ft_tip_depths <- ape::node.depth.edgelength(ft_tree)[seq_along(ft_tree$tip.label)]

stopifnot(sum(ft_clone_size) == ft_total_tips)
ft_sorted_sizes <- sort(ft_clone_size, decreasing = TRUE)
ft_rank_df      <- data.frame(rank = seq_along(ft_sorted_sizes), size = ft_sorted_sizes)

write.csv(ft_rank_df,
          file.path(out_dir, "fig4_panel_A.csv"), row.names = FALSE)

ft_dT     <- max(ft_tip_depths) - T0
ft_lambda <- log(ft_total_tips / ft_N0) / ft_dT
ft_p_null <- exp(-ft_lambda * ft_dT)
cat(sprintf("full tree: %d founders at T0=%.2f, %s total tips, fitted p_null=%.5f\n",
            ft_N0, T0, format(ft_total_tips, big.mark = ","), ft_p_null))


# How many founders are required to populate 50% of the embryo?
k_for_frac <- function(sizes, frac) {
  s <- sort(sizes, decreasing = TRUE)
  which(cumsum(s) >= frac * sum(s))[1]
}
FRAC_B  <- 0.5
k_obs_b <- k_for_frac(ft_sorted_sizes, FRAC_B)


#-------------------------------------------------------------------------#
# Analyse static clonal dominance on backbone tree, panels supporting a-c
bt_split      <- castor::split_tree_at_height(bt, height = T0)
bt_clone_size <- vapply(bt_split$subtrees, function(s) length(s$tree$tip.label), integer(1))
bt_N0         <- length(bt_clone_size)

bt_total_tips <- length(bt$tip.label)
bt_tip_depths <- ape::node.depth.edgelength(bt)[seq_along(bt$tip.label)]

stopifnot(sum(bt_clone_size) == bt_total_tips)
bt_sorted_sizes <- sort(bt_clone_size, decreasing = TRUE)
bt_rank_df      <- data.frame(rank = seq_along(bt_sorted_sizes), size = bt_sorted_sizes)

write.csv(bt_rank_df,
          file.path(out_dir_backbone, "fig_s13_panel_A_bt.csv"), row.names = FALSE)

bt_dT     <- max(bt_tip_depths) - T0
bt_lambda <- log(bt_total_tips / bt_N0) / bt_dT
bt_p_null <- exp(-bt_lambda * bt_dT)
cat(sprintf("backbone tree: %d founders at T0=%.2f, %s total tips, fitted p_null=%.5f\n",
            bt_N0, T0, format(bt_total_tips, big.mark = ","), bt_p_null))


# How many founders are required to populate 50% of the embryo?
k_obs_b_bt <- k_for_frac(bt_sorted_sizes, FRAC_B)

#-------------------------------------------------------------------------#
# Panel b: how many founders (observed) are needed to reach 50% of all E13.5
# cells, against its Yule-null Monte Carlo distribution. 

set.seed(1)   
R_b <- 10000
k_null_b <- vapply(seq_len(R_b), function(i) k_for_frac(1 + rgeom(ft_N0, ft_p_null), FRAC_B), integer(1))
p_val_b  <- mean(k_null_b <= k_obs_b)   # one-sided: P(null needs as few or fewer founders than observed)

write.csv(data.frame(k_null = k_null_b, k_obs = k_obs_b, p_value = p_val_b, N0 = ft_N0, frac = FRAC_B),
          file.path(out_dir, "fig4_panel_B_founders_for_half.csv"), row.names = FALSE)



# Panel b Supp on backbone
set.seed(1)   
R_b <- 10000
k_null_b_bt <- vapply(seq_len(R_b), function(i) k_for_frac(1 + rgeom(bt_N0, bt_p_null), FRAC_B), integer(1))
p_val_b_bt  <- mean(k_null_b_bt <= k_obs_b_bt)   # one-sided: P(null needs as few or fewer founders than observed)

write.csv(data.frame(k_null = k_null_b, k_obs = k_obs_b, p_value = p_val_b, N0 = ft_N0, frac = FRAC_B),
          file.path(out_dir_backbone, "fig_s13_panel_B_founders_for_half.csv"), row.names = FALSE)


#-------------------------------------------------------------------------#
# Panel c: Gini of founder clone sizes vs its Monte-Carlo Yule null.

gini_obs <- ineq::Gini(ft_clone_size)

set.seed(1)
R <- 10000
null_ginis <- vapply(seq_len(R), function(i) {
  ineq::Gini(1 + rgeom(ft_N0, ft_p_null))
}, numeric(1))

p_gini_val <- mean(null_ginis >= gini_obs)
gini_df <- data.frame(gini = null_ginis)

write.csv(data.frame(gini_df, gini_obs = gini_obs), file.path(out_dir, "fig4_panel_C_gini_vs_null.csv"), row.names = FALSE)


## and on backbone tree
gini_obs_bt <- ineq::Gini(bt_clone_size)

set.seed(1)
R <- 10000
null_ginis <- vapply(seq_len(R), function(i) {
  ineq::Gini(1 + rgeom(bt_N0, bt_p_null))
}, numeric(1))

p_gini_val <- mean(null_ginis >= gini_obs_bt)
gini_df_bt <- data.frame(gini = null_ginis)

write.csv(data.frame(gini_df_bt, gini_obs = gini_obs_bt), file.path(out_dir_backbone, "fig_s13_gini_vs_null.csv"), row.names = FALSE)



#-------------------------------------------------------------------------#
# Panel d: Dynamic analysis of when rise in Gini occurs over time
# uses the backbone tree!
total_tips <- length(bt$tip.label)

tip_founder <- integer(length(bt$tip.label))
names(tip_founder) <- bt$tip.label
split_bt <- castor::split_tree_at_height(bt, height = T0)

for (i in seq_along(split_bt$subtrees)) {
  tip_founder[split_bt$subtrees[[i]]$tree$tip.label] <- i
}

# get cell type info for tips already here because we will compute the gini
# over time directly globally and for all cell types. 
meta <- read.delim("support_data/cell_metadata.v8.txt.gz", stringsAsFactors = FALSE)
traj_of <- setNames(meta$major_trajectory, meta$cell_id)

type_of_tip <- traj_of[bt$tip.label]
all_types   <- sort(unique(na.omit(type_of_tip)))

## Get gini null over time

lambda <- log(total_tips / bt_N0) / bt_dT
T_max  <- max(ft_tip_depths) - 0.01
n_grid <- 50
R_gini_time <- 10000

T_seq_fixed <- seq(T0, T_max, length.out = n_grid)

gini_null_time_df <- do.call(rbind, lapply(T_seq_fixed, function(T) {
  p_T <- exp(-lambda * (T - T0))
  null_draws <- vapply(seq_len(R_gini_time), function(i) ineq::Gini(1 + rgeom(bt_N0, p_T)), numeric(1))
  ci <- quantile(null_draws, c(0.025, 0.975))
  data.frame(T = T, gini_null_median = median(null_draws), gini_null_lo = ci[[1]], gini_null_hi = ci[[2]])
}))


## Get observed gini over time
N0 = bt_N0
raw_obs_end <- list()

gini_over_time <- lapply(seq_along(T_seq_fixed), function(t_idx) {
  print(t_idx)
  T <- T_seq_fixed[t_idx]
  split_T <- castor::split_tree_at_height(bt, height = T)
  sizes   <- vapply(split_T$subtrees, function(s) length(s$tree$tip.label), integer(1))
  # One lookup per T on the WHOLE vector of first-tips at once (single hash
  # pass over tip_founder's names), not one tip_founder[[name]] call per
  # lineage.
  first_tip <- vapply(split_T$subtrees, function(s) s$tree$tip.label[1], character(1))
  lineage_founder <- tip_founder[first_tip]
  counts  <- tabulate(lineage_founder, nbins = N0)   # clone_size_i(T) per founder, in founder order
  stopifnot(all(counts >= 1))        # every founder should still have >=1 living lineage at T > T0
  whole   <- data.frame(T = T, gini = ineq::Gini(counts))
  
  # Per-type: cross-tab this SAME cut's lineages x cell types in one
  # vectorized pass, plus
  # lineage_founder already computed above
  all_tips   <- unlist(lapply(split_T$subtrees, function(s) s$tree$tip.label))
  lineage_id <- rep(seq_along(split_T$subtrees), times = sizes)
  idx        <- match(all_tips, bt$tip.label)
  type_tab   <- table(lineage = factor(lineage_id, levels = seq_along(split_T$subtrees)),
                      type    = factor(type_of_tip[idx], levels = all_types))
  has_type   <- type_tab > 0
  by_type <- do.call(rbind, lapply(all_types, function(ty) {
    print(ty)
    fc_full <- tabulate(lineage_founder[has_type[, ty]], nbins = N0)   # zero-inclusive, length N0
    if (t_idx == n_grid) raw_obs_end[[ty]] <<- fc_full
    fc_nz <- fc_full[fc_full > 0]
    data.frame(T = T, major_trajectory = ty, n_founders = length(fc_nz),
               gini = if (length(fc_nz) >= 2) ineq::Gini(fc_nz) else NA_real_
               )
  }))
  
  list(whole = whole, by_type = by_type)
})
gini_fixed_df      <- do.call(rbind, lapply(gini_over_time, `[[`, "whole"))
gini_fixed_type_df <- do.call(rbind, lapply(gini_over_time, `[[`, "by_type"))
gini_fixed_type_df <- gini_fixed_type_df[!is.na(gini_fixed_type_df$gini), ]

write.csv(gini_null_time_df, file.path(out_dir, "fig4_panel_D_gini_over_time_null.csv"), row.names = FALSE)
write.csv(gini_fixed_df, file.path(out_dir, "fig4_panel_D_gini_over_time.csv"), row.names = FALSE)


#-------------------------------------------------------------------------#
# Panel E: Observed versus expected contribution to different cell types.
#
bt_root      <- length(bt$tip.label) + 1L
bt_root_kids <- bt$edge[bt$edge[, 1] == bt_root, 2]
stopifnot(length(bt_root_kids) == 2)   # bifurcating root by construction

tips_under_root_kid <- function(node) {
  if (node <= length(bt$tip.label)) return(bt$tip.label[node])
  castor::get_subtree_at_node(bt, node - length(bt$tip.label))$subtree$tip.label
}
root_kid_tips        <- lapply(bt_root_kids, tips_under_root_kid)
root_kid_sizes       <- vapply(root_kid_tips, length, integer(1))
larger_root_kid_tips <- root_kid_tips[[which.max(root_kid_sizes)]]

blastomere_larger_label  <- "B1"
blastomere_smaller_label <- "B2"
blast_labels  <- c(blastomere_larger_label, blastomere_smaller_label)
blast_display <- setNames(c("Blastomere A", "Blastomere B"), blast_labels)
disp_levels   <- unname(blast_display[blast_labels])
to_disp       <- function(x) factor(blast_display[as.character(x)], levels = disp_levels)

# One blastomere label per founder, aligned with tip_founder's index space
# (split_bt$subtrees order, i.e. 1:bt_N0); carried down to one label per tip
# through that same founder assignment.
founder_blast <- vapply(split_bt$subtrees, function(s)
  if (s$tree$tip.label[1] %in% larger_root_kid_tips) blastomere_larger_label else blastomere_smaller_label,
  character(1))
tip_blast <- founder_blast[tip_founder]

# N0/total_tips per blastomere, for the Yule-null fit below. 
blast_gini_df <- data.frame(
  blastomere = blast_labels,
  N0         = as.integer(table(founder_blast)[blast_labels]),
  total_tips = as.integer(table(tip_blast)[blast_labels])
)

# Cell types with >=100 cells in BOTH blastomeres -- a type clearing the bar
# in only one can't be compared on equal footing between the two.
type_n_cells_by_blast <- table(blastomere = tip_blast, type = as.character(type_of_tip))
n_by_blast <- function(ty) setNames(type_n_cells_by_blast[, ty], rownames(type_n_cells_by_blast))
cons_types <- Filter(function(ty) all(n_by_blast(ty) >= 100), all_types)
cons_type_order <- names(sort(colSums(type_n_cells_by_blast[, cons_types, drop = FALSE]), decreasing = TRUE))
cat(sprintf("panel E: %d / %d trajectories have >=100 cells in both blastomeres\n",
            length(cons_types), length(all_types)))

# Gini(T) per (blastomere, cell type), on the same time grid as panel d
# (T_seq_fixed). Same split-and-cross-tab pattern as gini_over_time() above,
# re-run here stratified by blastomere and restricted to cons_types.
gini_by_blast_list <- lapply(seq_along(T_seq_fixed), function(t_idx) {
  T <- T_seq_fixed[t_idx]
  split_T   <- castor::split_tree_at_height(bt, height = T)
  sizes     <- vapply(split_T$subtrees, function(s) length(s$tree$tip.label), integer(1))
  first_tip <- vapply(split_T$subtrees, function(s) s$tree$tip.label[1], character(1))
  lineage_founder  <- tip_founder[first_tip]
  blast_of_lineage <- founder_blast[lineage_founder]

  all_tips   <- unlist(lapply(split_T$subtrees, function(s) s$tree$tip.label))
  lineage_id <- rep(seq_along(split_T$subtrees), times = sizes)
  idx        <- match(all_tips, bt$tip.label)
  type_tab   <- table(lineage = factor(lineage_id, levels = seq_along(split_T$subtrees)),
                      type    = factor(type_of_tip[idx], levels = cons_types))
  has_type   <- type_tab > 0

  do.call(rbind, lapply(blast_labels, function(b) {
    keep_b <- blast_of_lineage == b
    do.call(rbind, lapply(cons_types, function(ty) {
      fc_nz <- tabulate(lineage_founder[keep_b & has_type[, ty]], nbins = bt_N0)
      fc_nz <- fc_nz[fc_nz > 0]
      data.frame(T = T, blastomere = b, major_trajectory = ty,
                 gini = if (length(fc_nz) >= 2) ineq::Gini(fc_nz) else NA_real_)
    }))
  }))
})
gini_fixed_type_blast_df <- do.call(rbind, gini_by_blast_list)
gini_fixed_type_blast_df <- gini_fixed_type_blast_df[!is.na(gini_fixed_type_blast_df$gini), ]

# Median + 95% interval of a numeric vector, NA-safe.
summarize_row <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(c(median = NA_real_, lo = NA_real_, hi = NA_real_))
  ci <- quantile(x, c(0.025, 0.975))
  c(median = median(x), lo = unname(ci[1]), hi = unname(ci[2]))
}

# Observed Gini(T) per type per blastomere as above against a null
# that does not depend on the trajectory's own cell count -- only on the
# blastomere's founder count and total cell count.
set.seed(1)
R_blast_null <- 10000
gini_null_time_blast_only_df <- do.call(rbind, lapply(blast_labels, function(b) {
  N0_b     <- blast_gini_df$N0[blast_gini_df$blastomere == b]
  total_b  <- blast_gini_df$total_tips[blast_gini_df$blastomere == b]
  lambda_b <- log(total_b / N0_b) / bt_dT
  do.call(rbind, lapply(T_seq_fixed, function(T) {
    p_T <- exp(-lambda_b * (T - T0))
    null_draws <- vapply(seq_len(R_blast_null), function(i) ineq::Gini(1 + rgeom(N0_b, p_T)), numeric(1))
    s <- summarize_row(null_draws)
    data.frame(blastomere = b, T = T,
               gini_null_median = s[["median"]], gini_null_lo = s[["lo"]], gini_null_hi = s[["hi"]])
  }))
}))

# same null value applies to every trajectory at a given (blastomere, T) --
# merge on blastomere+T only, NOT major_trajectory
cons_df_blastnull <- merge(gini_fixed_type_blast_df, gini_null_time_blast_only_df, by = c("blastomere", "T"))
cons_df_blastnull$diff    <- cons_df_blastnull$gini - cons_df_blastnull$gini_null_median
# excess over the null's 97.5% upper bound rather than its median -- a more
# conservative "beats the null envelope" comparison
cons_df_blastnull$diff_hi <- cons_df_blastnull$gini - cons_df_blastnull$gini_null_hi
cons_df_blastnull$blast <- to_disp(cons_df_blastnull$blastomere)
cons_df_blastnull$major_trajectory <- factor(cons_df_blastnull$major_trajectory, levels = cons_type_order)
write.csv(cons_df_blastnull, file.path(out_dir, "fig4_panel_F_gini_excess.csv"), row.names = FALSE)

#-------------------------------------------------------------------------#
# Panel G: founder x cell-type heatmap (full tree), fill = log2(obs/exp)
# relative to each founder's overall contribution to the embryo. Cells with
# < MIN_CELLS_G observed cells are masked to NA (grey) -- too few cells for
# a log2 ratio to be meaningful.
MIN_CELLS_G <- 10
MIN_VISIBLE_TYPES_G <- 1   # founders with fewer visible cell types go in a terminal grey block

# Founder's blastomere = which root child its first tip falls under (T0 cut
# sits entirely inside one root child), same logic as panel E/F but applied
# to the full tree.
ft_root      <- length(ft_tree$tip.label) + 1L
ft_root_kids <- ft_tree$edge[ft_tree$edge[, 1] == ft_root, 2]
stopifnot(length(ft_root_kids) == 2)   # bifurcating root by construction

ft_tips_under_root_kid <- function(node) {
  if (node <= length(ft_tree$tip.label)) return(ft_tree$tip.label[node])
  castor::get_subtree_at_node(ft_tree, node - length(ft_tree$tip.label))$subtree$tip.label
}
ft_root_kid_tips        <- lapply(ft_root_kids, ft_tips_under_root_kid)
ft_root_kid_sizes       <- vapply(ft_root_kid_tips, length, integer(1))
ft_larger_root_kid_tips <- ft_root_kid_tips[[which.max(ft_root_kid_sizes)]]

ft_founder_blast <- vapply(ft_split$subtrees, function(s)
  if (s$tree$tip.label[1] %in% ft_larger_root_kid_tips) blastomere_larger_label else blastomere_smaller_label,
  character(1))
blast_colors <- setNames(c("#3573b9", "#e8863b"), disp_levels)   # pairs with to_disp() above

# per-founder x cell-type count matrix
g_all_tips   <- unlist(lapply(ft_split$subtrees, function(s) s$tree$tip.label))
g_all_types  <- sort(unique(na.omit(traj_of[g_all_tips])))
g_ncells_all <- table(factor(traj_of[g_all_tips], levels = g_all_types))
g_type_order <- names(sort(g_ncells_all, decreasing = TRUE))   # largest cell types first

M_all <- t(vapply(ft_split$subtrees, function(s)
  as.numeric(table(factor(traj_of[s$tree$tip.label], levels = g_all_types))), numeric(length(g_all_types))))
colnames(M_all) <- g_all_types

founder_total <- rowSums(M_all)
type_total    <- colSums(M_all)
N_grand       <- sum(M_all)
expected_all  <- outer(founder_total, type_total) / N_grand
log2fc_raw    <- ifelse(M_all > 0, log2(M_all / expected_all), NA_real_)
log2fc_all    <- ifelse(M_all >= MIN_CELLS_G, log2fc_raw, NA_real_)

# order founders by their strongest VISIBLE cell type; must use the masked
# matrix so a founder never gets ordered by a tile that's actually drawn
# grey. NAs (too sparse to place meaningfully) sort last by default.
has_enough_visible <- rowSums(!is.na(log2fc_all)) >= MIN_VISIBLE_TYPES_G
if (!all(has_enough_visible)) {
  cat(sprintf("panel G: %d founder(s) placed in a terminal block (< %d visible cell types)\n",
              sum(!has_enough_visible), MIN_VISIBLE_TYPES_G))
}
best_type_idx <- rep(NA_integer_, nrow(M_all))
best_type_idx[has_enough_visible] <- apply(log2fc_all[has_enough_visible, , drop = FALSE], 1, which.max)

argmax_type <- rep(NA_character_, nrow(M_all))
argmax_type[has_enough_visible] <- colnames(M_all)[best_type_idx[has_enough_visible]]

maxval <- rep(NA_real_, nrow(M_all))
maxval[has_enough_visible] <- log2fc_all[cbind(which(has_enough_visible), best_type_idx[has_enough_visible])]

founder_ord <- order(match(argmax_type, g_type_order), -maxval)

g_col_ids <- sprintf("f%03d", seq_along(founder_ord))
g_long <- data.frame(
  founder   = factor(rep(g_col_ids, times = ncol(M_all)), levels = g_col_ids),
  celltype  = factor(rep(colnames(M_all), each = length(founder_ord)), levels = rev(g_type_order)),
  log2fc    = as.vector(log2fc_all[founder_ord, , drop = FALSE]),
  n_cells   = as.vector(M_all[founder_ord, , drop = FALSE]),
  type_n    = rep(type_total[colnames(M_all)], each = length(founder_ord)),
  founder_n = rep(founder_total[founder_ord], times = ncol(M_all)),
  blast     = rep(ft_founder_blast[founder_ord], times = ncol(M_all))
)
write.csv(g_long, "figures_data/fig4_panel_E_prevalent_founders_heatmap.csv", row.names = FALSE)

#-------------------------------------------------------------------------#
# Sensitivity for Blastomeres A, B, later to be found in Fig S13 D-F.

#-------------------------------------------------------------------------#
# Panel D: rank-abundance of founder clone sizes, one blastomere per facet
# (full tree). Also derives the top-k cutoff that panel E's null is compared
# against (founders for FRAC_B of THIS blastomere's own cells).
ft_blast_stats <- lapply(blast_labels, function(b) {
  sizes_b  <- ft_clone_size[ft_founder_blast == b]
  N0_b     <- length(sizes_b)
  total_b  <- sum(sizes_b)
  sorted_b <- sort(sizes_b, decreasing = TRUE)
  k_obs_b  <- k_for_frac(sorted_b, FRAC_B)

  lambda_b <- log(total_b / N0_b) / ft_dT
  p_null_b <- exp(-lambda_b * ft_dT)
  set.seed(1)
  R_b <- 10000
  k_null_b <- vapply(seq_len(R_b), function(i) k_for_frac(1 + rgeom(N0_b, p_null_b), FRAC_B), integer(1))
  p_val_b  <- mean(k_null_b <= k_obs_b)   # one-sided: P(null needs as few or fewer founders than observed)

  list(blastomere = b, N0 = N0_b, total_tips = total_b, sorted_sizes = sorted_b,
       k_obs = k_obs_b, k_null = k_null_b, p_val = p_val_b)
})
names(ft_blast_stats) <- blast_labels

ft_blast_rank_df <- do.call(rbind, lapply(ft_blast_stats, function(x)
  data.frame(blastomere = x$blastomere, rank = seq_along(x$sorted_sizes), size = x$sorted_sizes)))
ft_blast_rank_df$blast <- to_disp(ft_blast_rank_df$blastomere)

ft_blast_cutoff_df <- do.call(rbind, lapply(ft_blast_stats, function(x)
  data.frame(blastomere = x$blastomere, k = x$k_obs, N0 = x$N0, total_tips = x$total_tips)))
ft_blast_cutoff_df$blast <- to_disp(ft_blast_cutoff_df$blastomere)

write.csv(ft_blast_rank_df, file.path(out_dir, "fig_s13_panel_D_rank_abundance.csv"), row.names = FALSE)
write.csv(ft_blast_cutoff_df, file.path(out_dir, "fig_s13_panel_D_rank_abundance_cutoff.csv"), row.names = FALSE)

#-------------------------------------------------------------------------#
# Panel E: founders needed for FRAC_B of each blastomere's own E13.5 cells,
# against its own Yule-null Monte Carlo distribution (full tree). Reuses
# ft_blast_stats computed for panel D above.
ft_blast_null_long <- do.call(rbind, lapply(ft_blast_stats, function(x)
  data.frame(blastomere = x$blastomere, k_null = x$k_null)))
ft_blast_null_long$blast <- to_disp(ft_blast_null_long$blastomere)

ft_blast_obs_df <- do.call(rbind, lapply(ft_blast_stats, function(x)
  data.frame(blastomere = x$blastomere, k_obs = x$k_obs, p_val = x$p_val, N0 = x$N0, total_tips = x$total_tips)))
ft_blast_obs_df$blast <- to_disp(ft_blast_obs_df$blastomere)

write.csv(ft_blast_null_long, file.path(out_dir, "fig_s13_panel_E_founders_for_half_null.csv"), row.names = FALSE)
write.csv(ft_blast_obs_df, file.path(out_dir, "fig_s13_panel_E_founders_for_half_obs.csv"), row.names = FALSE)

#-------------------------------------------------------------------------#
# Panel F: Gini of founder clone sizes vs its Monte-Carlo Yule null, one
# blastomere per facet (backbone tree).

# tabulate(), not bt_clone_size -- that comes from a separate
# split_tree_at_height() call (bt_split) and isn't guaranteed to align by
# index with founder_blast, which comes from split_bt/tip_founder.
bt_founder_clone_size <- tabulate(tip_founder, nbins = bt_N0)

bt_blast_stats <- lapply(blast_labels, function(b) {
  sizes_b    <- bt_founder_clone_size[founder_blast == b]
  N0_b       <- length(sizes_b)
  total_b    <- sum(sizes_b)
  lambda_b   <- log(total_b / N0_b) / bt_dT
  p_null_b   <- exp(-lambda_b * bt_dT)

  gini_obs_b <- ineq::Gini(sizes_b)
  set.seed(1)
  R <- 10000
  null_ginis_b <- vapply(seq_len(R), function(i) ineq::Gini(1 + rgeom(N0_b, p_null_b)), numeric(1))

  list(blastomere = b, N0 = N0_b, total_tips = total_b, gini_obs = gini_obs_b, null_ginis = null_ginis_b)
})
names(bt_blast_stats) <- blast_labels

bt_blast_null_long <- do.call(rbind, lapply(bt_blast_stats, function(x)
  data.frame(blastomere = x$blastomere, gini = x$null_ginis)))
bt_blast_null_long$blast <- to_disp(bt_blast_null_long$blastomere)

# z-score is comparable across blastomeres despite their different N; raw
# Gini or the Monte-Carlo p (which saturates at 1/R) is not.
bt_blast_obs_df <- do.call(rbind, lapply(bt_blast_stats, function(x)
  data.frame(blastomere = x$blastomere, N0 = x$N0, total_tips = x$total_tips,
             gini_obs = x$gini_obs, gini_z = (x$gini_obs - mean(x$null_ginis)) / sd(x$null_ginis),
             p_gini = mean(x$null_ginis >= x$gini_obs))))
bt_blast_obs_df$blast <- to_disp(bt_blast_obs_df$blastomere)

write.csv(bt_blast_null_long, file.path(out_dir_backbone, "fig_s13_panel_F_gini_vs_null_null.csv"), row.names = FALSE)
write.csv(bt_blast_obs_df, file.path(out_dir_backbone, "fig_s13_panel_F_gini_vs_null_obs.csv"), row.names = FALSE)

#-------------------------------------------------------------------------#
# Fig S14 panels A-F.
#
# A/B: blastomere-split Gini(T) across three differently-dated trees (main +
# two root-split-ratio sensitivity trees). C/D: sweep T0 itself on the
# (undivided) backbone tree. E: founder x cell-type heatmap, backbone tree.
# F: per-cell-type Gini consistency between the full placement tree and the
# backbone tree. Shared time grid (T_max_sens/n_grid) for A-D so all four
# panels share one x-axis.
SENS_TREE_63_37 <- "tree_building/results/3-dated-tree-sensitivity/merged_minB2h_lineage_constrained_sidefrac63-37.nwk"
SENS_TREE_58_42 <- "tree_building/results/3-dated-tree-sensitivity/merged_minB2h_lineage_constrained_sidefrac58-42.nwk"

T_max_sens <- max(bt_tip_depths) - 0.01
R_gini_time_sens <- 1000   # coarser than R_gini_time (10000) above -- matches source

# root-split blastomere assignment for an arbitrary tree/T0 -- founder_blast/
# split_bt above are bt-specific, this generalizes it to the sensitivity trees.
assign_blast_founders <- function(tree, T0) {
  root <- length(tree$tip.label) + 1L
  kids <- tree$edge[tree$edge[, 1] == root, 2]
  stopifnot(length(kids) == 2)
  tips_under <- function(node) {
    if (node <= length(tree$tip.label)) return(tree$tip.label[node])
    castor::get_subtree_at_node(tree, node - length(tree$tip.label))$subtree$tip.label
  }
  kid_tips    <- lapply(kids, tips_under)
  larger_tips <- kid_tips[[which.max(vapply(kid_tips, length, integer(1)))]]

  split <- castor::split_tree_at_height(tree, height = T0)
  blast <- vapply(split$subtrees, function(s)
    if (s$tree$tip.label[1] %in% larger_tips) blastomere_larger_label else blastomere_smaller_label,
    character(1))
  list(split = split, blast = blast)
}

# Gini(T) per blastomere for one tree, over a shared time grid.
gini_over_time_by_blast <- function(tree, blast_founders, T_seq) {
  N0 <- length(blast_founders$blast)
  tip_founder_i <- integer(length(tree$tip.label))
  names(tip_founder_i) <- tree$tip.label
  for (i in seq_along(blast_founders$split$subtrees)) {
    tip_founder_i[blast_founders$split$subtrees[[i]]$tree$tip.label] <- i
  }

  do.call(rbind, lapply(T_seq, function(T) {
    split_T   <- castor::split_tree_at_height(tree, height = T)
    first_tip <- vapply(split_T$subtrees, function(s) s$tree$tip.label[1], character(1))
    counts    <- tabulate(tip_founder_i[first_tip], nbins = N0)
    stopifnot(all(counts >= 1))
    do.call(rbind, lapply(blast_labels, function(b) {
      idx <- which(blast_founders$blast == b)
      data.frame(T = T, blastomere = b, gini = ineq::Gini(counts[idx]))
    }))
  }))
}

T_seq_sens <- seq(T0, T_max_sens, length.out = n_grid)

# "main" tree series -- computed here (needed by both A and B below) but its
# own CSV isn't written until the panel B section, to keep write.csv() calls
# in A/B/C/D/E/F order.
gini_blast_main <- gini_over_time_by_blast(bt, list(split = split_bt, blast = founder_blast), T_seq_sens)
gini_blast_main$tree  <- "main"
gini_blast_main$blast <- to_disp(gini_blast_main$blastomere)

#-------------------------------------------------------------------------#
# Panel A: main + two sensitivity trees, combined.
sens_tree_6337   <- ape::read.tree(SENS_TREE_63_37)
assign_6337      <- assign_blast_founders(sens_tree_6337, T0)
gini_blast_6337  <- gini_over_time_by_blast(sens_tree_6337, assign_6337, T_seq_sens)
gini_blast_6337$tree <- "sensitivity_63-37"

sens_tree_5842   <- ape::read.tree(SENS_TREE_58_42)
assign_5842      <- assign_blast_founders(sens_tree_5842, T0)
gini_blast_5842  <- gini_over_time_by_blast(sens_tree_5842, assign_5842, T_seq_sens)
gini_blast_5842$tree <- "sensitivity_58-42"

gini_blast_sens_df <- rbind(gini_blast_main[, 1:4], gini_blast_6337, gini_blast_5842)
gini_blast_sens_df$blast <- to_disp(gini_blast_sens_df$blastomere)
write.csv(gini_blast_sens_df, file.path(out_dir, "fig_s14_panel_A_gini_per_blastomere_sensitivity.csv"), row.names = FALSE)

#-------------------------------------------------------------------------#
# Panel B: "main" tree only -- reuses gini_blast_main computed above instead
# of re-deriving it.
write.csv(gini_blast_main, file.path(out_dir, "fig_s14_panel_B_gini_emergence.csv"), row.names = FALSE)

#-------------------------------------------------------------------------#
# Panels C/D: T0-choice sensitivity, whole backbone tree (no blastomere
# split). Gini(T) obs + Yule null per candidate T0.
T0_sweep_values <- c(5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 10)

compute_gini_for_T0 <- function(T0_i, tree, total_tips_val, T_max, n_grid, R) {
  split_T0 <- castor::split_tree_at_height(tree, height = T0_i)
  N0 <- length(split_T0$subtrees)
  tip_founder_i <- integer(length(tree$tip.label))
  names(tip_founder_i) <- tree$tip.label
  for (i in seq_along(split_T0$subtrees)) tip_founder_i[split_T0$subtrees[[i]]$tree$tip.label] <- i

  dT_i     <- T_max - T0_i
  lambda_i <- log(total_tips_val / N0) / dT_i
  T_seq_i  <- seq(T0_i, T_max, length.out = n_grid)

  obs_df <- do.call(rbind, lapply(T_seq_i, function(T) {
    split_T   <- castor::split_tree_at_height(tree, height = T)
    first_tip <- vapply(split_T$subtrees, function(s) s$tree$tip.label[1], character(1))
    counts    <- tabulate(tip_founder_i[first_tip], nbins = N0)
    stopifnot(all(counts >= 1))
    data.frame(T0 = T0_i, T = T, gini = ineq::Gini(counts))
  }))

  null_df <- do.call(rbind, lapply(T_seq_i, function(T) {
    p_T <- exp(-lambda_i * (T - T0_i))
    null_draws <- vapply(seq_len(R), function(i) ineq::Gini(1 + rgeom(N0, p_T)), numeric(1))
    ci <- quantile(null_draws, c(0.025, 0.975))
    data.frame(T0 = T0_i, T = T, gini_null_median = median(null_draws), gini_null_lo = ci[[1]], gini_null_hi = ci[[2]])
  }))

  list(obs = obs_df, null = null_df)
}

T0_sweep_results <- lapply(T0_sweep_values, compute_gini_for_T0,
                            tree = bt, total_tips_val = bt_total_tips, T_max = T_max_sens,
                            n_grid = n_grid, R = R_gini_time_sens)

T0_sweep_obs_df  <- do.call(rbind, lapply(T0_sweep_results, `[[`, "obs"))
T0_sweep_null_df <- do.call(rbind, lapply(T0_sweep_results, `[[`, "null"))

write.csv(T0_sweep_obs_df,  file.path(out_dir, "fig_s14_panel_C_gini_over_time_T0sweep.csv"), row.names = FALSE)
write.csv(T0_sweep_null_df, file.path(out_dir, "fig_s14_panel_C_gini_over_time_T0sweep_null.csv"), row.names = FALSE)

# Panel D: excess of observed Gini over the null's 97.5% upper bound, same
# T0 sweep.
T0_sweep_diff_df <- data.frame(T0 = T0_sweep_obs_df$T0, T = T0_sweep_obs_df$T,
                                gini_diff = T0_sweep_obs_df$gini - T0_sweep_null_df$gini_null_hi)
write.csv(T0_sweep_diff_df, file.path(out_dir, "fig_s14_panel_D_gini_diff_to_null_T0sweep.csv"), row.names = FALSE)

#-------------------------------------------------------------------------#
# Panel E: same founder x cell-type heatmap as fig4 panel E above, backbone
# tree instead of full tree -- reuses split_bt/founder_blast/all_types
# already computed above instead of reloading a tree.
bt_all_tips   <- unlist(lapply(split_bt$subtrees, function(s) s$tree$tip.label))
bt_ncells_all <- table(factor(traj_of[bt_all_tips], levels = all_types))
bt_type_order <- names(sort(bt_ncells_all, decreasing = TRUE))

bt_M_all <- t(vapply(split_bt$subtrees, function(s)
  as.numeric(table(factor(traj_of[s$tree$tip.label], levels = all_types))), numeric(length(all_types))))
colnames(bt_M_all) <- all_types

bt_founder_total <- rowSums(bt_M_all)
bt_type_total    <- colSums(bt_M_all)
bt_N_grand       <- sum(bt_M_all)
bt_expected_all  <- outer(bt_founder_total, bt_type_total) / bt_N_grand
bt_log2fc_raw    <- ifelse(bt_M_all > 0, log2(bt_M_all / bt_expected_all), NA_real_)
bt_log2fc_all    <- ifelse(bt_M_all >= MIN_CELLS_G, bt_log2fc_raw, NA_real_)

bt_has_enough_visible <- rowSums(!is.na(bt_log2fc_all)) >= MIN_VISIBLE_TYPES_G
if (!all(bt_has_enough_visible)) {
  cat(sprintf("panel E (backbone): %d founder(s) placed in a terminal block (< %d visible cell types)\n",
              sum(!bt_has_enough_visible), MIN_VISIBLE_TYPES_G))
}
bt_best_type_idx <- rep(NA_integer_, nrow(bt_M_all))
bt_best_type_idx[bt_has_enough_visible] <- apply(bt_log2fc_all[bt_has_enough_visible, , drop = FALSE], 1, which.max)

bt_argmax_type <- rep(NA_character_, nrow(bt_M_all))
bt_argmax_type[bt_has_enough_visible] <- colnames(bt_M_all)[bt_best_type_idx[bt_has_enough_visible]]

bt_maxval <- rep(NA_real_, nrow(bt_M_all))
bt_maxval[bt_has_enough_visible] <- bt_log2fc_all[cbind(which(bt_has_enough_visible), bt_best_type_idx[bt_has_enough_visible])]

bt_founder_ord <- order(match(bt_argmax_type, bt_type_order), -bt_maxval)

bt_col_ids <- sprintf("f%03d", seq_along(bt_founder_ord))
bt_g_long <- data.frame(
  founder   = factor(rep(bt_col_ids, times = ncol(bt_M_all)), levels = bt_col_ids),
  celltype  = factor(rep(colnames(bt_M_all), each = length(bt_founder_ord)), levels = rev(bt_type_order)),
  log2fc    = as.vector(bt_log2fc_all[bt_founder_ord, , drop = FALSE]),
  n_cells   = as.vector(bt_M_all[bt_founder_ord, , drop = FALSE]),
  type_n    = rep(bt_type_total[colnames(bt_M_all)], each = length(bt_founder_ord)),
  founder_n = rep(bt_founder_total[bt_founder_ord], times = ncol(bt_M_all)),
  blast     = to_disp(rep(founder_blast[bt_founder_ord], times = ncol(bt_M_all)))
)
write.csv(bt_g_long, file.path(out_dir_backbone, "fig_s14_panel_E_prevalent_founders_heatmap_backbone.csv"), row.names = FALSE)

#-------------------------------------------------------------------------#
# Panel F: per-cell-type Gini consistency between the full placement tree
# and the backbone tree, one blastomere per facet. A (blastomere, type) pair
# needs >= min_cells_per_type cells and a defined Gini (>= 2 founders) in
# BOTH trees -- below that, Gini is a small-N artifact.
min_cells_per_type <- 10

ft_blast_type_gini_df <- do.call(rbind, lapply(blast_labels, function(b) {
  idx <- which(ft_founder_blast == b)
  M <- t(vapply(ft_split$subtrees[idx], function(s)
    as.numeric(table(factor(traj_of[s$tree$tip.label], levels = g_all_types))), numeric(length(g_all_types))))
  colnames(M) <- g_all_types
  do.call(rbind, lapply(g_all_types, function(ty) {
    counts  <- sort(M[, ty], decreasing = TRUE)
    n_cells <- sum(counts)
    if (n_cells == 0) return(NULL)
    nz <- counts[counts > 0]
    data.frame(blastomere = b, major_trajectory = ty, n_cells = n_cells,
               n_founders = length(nz), gini = if (length(nz) >= 2) ineq::Gini(nz) else NA_real_)
  }))
}))

bt_blast_type_gini_df <- do.call(rbind, lapply(blast_labels, function(b) {
  idx <- which(founder_blast == b)
  M <- t(vapply(split_bt$subtrees[idx], function(s)
    as.numeric(table(factor(traj_of[s$tree$tip.label], levels = all_types))), numeric(length(all_types))))
  colnames(M) <- all_types
  do.call(rbind, lapply(all_types, function(ty) {
    counts  <- sort(M[, ty], decreasing = TRUE)
    n_cells <- sum(counts)
    if (n_cells == 0) return(NULL)
    nz <- counts[counts > 0]
    data.frame(blastomere = b, major_trajectory = ty, n_cells = n_cells,
               n_founders = length(nz), gini = if (length(nz) >= 2) ineq::Gini(nz) else NA_real_)
  }))
}))

type_consistency_all_df <- merge(
  ft_blast_type_gini_df[, c("blastomere", "major_trajectory", "gini", "n_founders", "n_cells")],
  bt_blast_type_gini_df[, c("blastomere", "major_trajectory", "gini", "n_founders", "n_cells")],
  by = c("blastomere", "major_trajectory"), suffixes = c("_full", "_backbone")
)
type_consistency_all_df$min_founders <- pmin(type_consistency_all_df$n_founders_full, type_consistency_all_df$n_founders_backbone)

keep_pair <- type_consistency_all_df$n_cells_full >= min_cells_per_type &
             type_consistency_all_df$n_cells_backbone >= min_cells_per_type &
             !is.na(type_consistency_all_df$gini_full) & !is.na(type_consistency_all_df$gini_backbone)
type_consistency_df <- type_consistency_all_df[keep_pair, ]

n_type_pairs_dropped <- length(union(paste(ft_blast_type_gini_df$blastomere, ft_blast_type_gini_df$major_trajectory),
                                     paste(bt_blast_type_gini_df$blastomere, bt_blast_type_gini_df$major_trajectory))) - nrow(type_consistency_all_df)
cat(sprintf("fig s14: %d (blastomere x type) pairs kept (>= %d cells + defined Gini in both trees); dropped %d not in both, %d filtered\n",
            nrow(type_consistency_df), min_cells_per_type, n_type_pairs_dropped, sum(!keep_pair)))

type_consistency_df$blast <- to_disp(type_consistency_df$blastomere)
write.csv(type_consistency_df, file.path(out_dir, "fig_s14_panel_F_consistency_percelltype.csv"), row.names = FALSE)

# Spearman rho per blastomere between full-tree and backbone-tree Gini,
# summarizing how well type_consistency_df's points track the diagonal.
rho_df <- do.call(rbind, lapply(blast_labels, function(b) {
  sub <- type_consistency_df[type_consistency_df$blastomere == b, ]
  data.frame(blastomere = b, n = nrow(sub),
             rho = if (nrow(sub) >= 3) cor(sub$gini_full, sub$gini_backbone, method = "spearman") else NA_real_)
})) 
rho_df$blast <- to_disp(rho_df$blastomere)
write.csv(rho_df, file.path(out_dir, "fig_s14_panel_F_consistency_percelltype_rho.csv"), row.names = FALSE)

