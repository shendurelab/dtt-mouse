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
