library(castor)
library(ggplot2)
library(patchwork)
library(scales)
library(RColorBrewer) 

#-------------------------------------------------------------------------#
# Figure settings
base_size <- 10
fig_theme <- theme_minimal(base_size = base_size, base_family = "sans") +
  theme(
    plot.title         = element_text(size = base_size, face = "bold", margin = margin(b = 2)),
    plot.subtitle      = element_text(size = base_size - 1.5, color = "grey30", margin = margin(b = 3), lineheight = 0.95),
    axis.title         = element_text(size = base_size),
    axis.text          = element_text(size = base_size - 1, color = "black"),
    legend.text        = element_text(size = base_size - 1.5),
    legend.title        = element_blank(),
    legend.key.size    = unit(2.5, "mm"),
    legend.margin      = margin(t = 0, b = 0),
    legend.box.spacing = unit(1, "mm"),
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_line(linewidth = 0.2, color = "grey88"),
    axis.line          = element_line(linewidth = 0.3, color = "black"),
    axis.ticks         = element_line(linewidth = 0.3, color = "black"),
    plot.margin        = margin(2, 3, 2, 2)
  )

output_dir = "figures/sensitivity_dating/"

#-------------------------------------------------------------------------#
# Load tree


#backbone tree
tree_file = "tree_building/results/3-dated-tree/merged_minB2h_lineage_constrained.nwk"
bt = ape::read.tree(tree_file)

T0_values <-c(5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 10)
total_tips <- length(bt$tip.label)
T_max  <- max(ape::node.depth.edgelength(bt)[seq_along(bt$tip.label)]) - 0.01
n_grid <- 50
R_gini_time <- 1000

compute_gini_for_T0 <- function(T0, bt, total_tips, T_max, n_grid, R_gini_time) {
  split_T0 <- castor::split_tree_at_height(bt, height = T0)
  N0 <- length(split_T0$subtrees)
  
  tip_founder <- integer(length(bt$tip.label))
  names(tip_founder) <- bt$tip.label
  for (i in seq_along(split_T0$subtrees)) {
    tip_founder[split_T0$subtrees[[i]]$tree$tip.label] <- i
  }
  
  dT     <- T_max - T0
  lambda <- log(total_tips / N0) / dT
  T_seq  <- seq(T0, T_max, length.out = n_grid)
  
  obs_df <- do.call(rbind, lapply(T_seq, function(T) {
    split_T <- castor::split_tree_at_height(bt, height = T)
    first_tip <- vapply(split_T$subtrees, function(s) s$tree$tip.label[1], character(1))
    counts <- tabulate(tip_founder[first_tip], nbins = N0)
    stopifnot(all(counts >= 1))
    data.frame(T0 = T0, T = T, gini = ineq::Gini(counts))
  }))
  
  null_df <- do.call(rbind, lapply(T_seq, function(T) {
    p_T <- exp(-lambda * (T - T0))
    null_draws <- vapply(seq_len(R_gini_time), function(i) ineq::Gini(1 + rgeom(N0, p_T)),
                         numeric(1))
    ci <- quantile(null_draws, c(0.025, 0.975))
    data.frame(T0 = T0, T = T, gini_null_median = median(null_draws),
               gini_null_lo = ci[[1]], gini_null_hi = ci[[2]])
  }))
  
  list(obs = obs_df, null = null_df)
}

results <- lapply(T0_values, compute_gini_for_T0,
                  bt = bt, total_tips = total_tips, T_max = T_max,
                  n_grid = n_grid, R_gini_time = R_gini_time)

gini_obs_df  <- do.call(rbind, lapply(results, `[[`, "obs"))
gini_null_df <- do.call(rbind, lapply(results, `[[`, "null"))

write.csv(gini_obs_df,  paste0(output_dir, "sensitivity_founder_gini.csv"),  row.names = FALSE)
write.csv(gini_null_df, paste0(output_dir, "sensitivity_founder_gini_null.csv"), row.names =
            FALSE)


#gini_obs_df = read.csv(paste0(output_dir, "gini_over_time_T0sens.csv"))
#gini_null_df = read.csv(paste0(output_dir, "gini_nulls_over_time_T0sens.csv"))

t0_pal <- setNames(RColorBrewer::brewer.pal(length(T0_values), "Dark2"),
                   as.character(T0_values))

p_d <- ggplot() +
  geom_ribbon(data = gini_null_df,
              aes(x = T, ymin = gini_null_lo, ymax = gini_null_hi, fill = as.factor(T0)),
              alpha = 0.2, color = NA) +
  geom_line(data = gini_obs_df, aes(x = T, y = gini, color = as.factor(T0)), linewidth = 0.6)+
  geom_point(data = gini_obs_df, aes(x = T, y = gini, color = as.factor(T0)), size = 0.6) +
  scale_color_manual(values = t0_pal) +
  scale_fill_manual(values = t0_pal) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Days post-fertilization", y = "Gini coefficient") +
  fig_theme

p_d

gini_diff <- data.frame(gini_diff = gini_obs_df$gini - gini_null_df$gini_null_hi)
gini_diff$T0 <- gini_obs_df$T0
gini_diff$T <- gini_obs_df$T

gini_diff |> 
  group_by(T0) |>
summarise(max_diff = max(gini_diff))

gini_diff |> 
  group_by(T0) |>
  summarise(max_diff = max(gini_diff),
            half_max_diff = max(gini_diff)/2,
            half_max_time = T[which.max(gini_diff >= max(gini_diff)/2)])

ggsave(paste0(output_dir, "sensitivity_founder_gini.pdf"), p_d, height = 3.5, width = 4.5)


p_d2 <- ggplot() +
  geom_line(data = gini_diff, aes(x = T, y = gini_diff, color = as.factor(T0)), linewidth = 0.6)+
  geom_point(data = gini_diff, aes(x = T, y = gini_diff, color = as.factor(T0)), size = 0.6) +
  scale_color_manual(values = t0_pal) +
  scale_fill_manual(values = t0_pal) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Days post-fertilization", y = "Difference in Gini to Null") +
  fig_theme

p_d2

ggsave(paste0(output_dir, "sensitivity_founder_gini_diff.pdf"), p_d2, height = 3.5, width = 4.5)


### 
# Sensitivity with respect to dated tree
## Load trees
# load trees
tree_files = c("tree_building/results/3-dated-tree/merged_minB2h_lineage_constrained.nwk",
               "tree_building/results/3-dated-tree-sensitivity/merged/merged_time_tree_ge7_attempt1_minage_sourced_minB2h_l0.01_sidefrac63-37.nwk",
               "tree_building/results/3-dated-tree-sensitivity-5842/merged/merged_time_tree_ge7_attempt1_minage_sourced_minB2h_l0.01_sidefrac58-42.nwk")

trees = lapply(tree_files, function(x) {
  ape::read.tree(x)
})

names(trees) = c("main", "sensitivity_63-37", "sensitivity_58-42")

T0_values <- 7.0
results_per_tree <- lapply(trees, compute_gini_for_T0,
                  T0 = T0_values, total_tips = total_tips, T_max = T_max,
                  n_grid = n_grid, R_gini_time = R_gini_time)


gini_per_tree_df  <- do.call(rbind, lapply(results_per_tree, `[[`, "obs"))
gini_per_tree_df$tree <- rep(names(trees), each = nrow(gini_per_tree_df) / length(trees))
gini_per_tree_null_df <- do.call(rbind, lapply(results_per_tree, `[[`, "null"))

write.csv(gini_per_tree_df,  paste0(output_dir, "gini_over_time_perTree.csv"),  row.names = FALSE)
write.csv(gini_per_tree_null_df, paste0(output_dir, "gini_nulls_over_time_perTree.csv"), row.names =
            FALSE)

tree_pal <- c("#2B8C2C", "#7B4FA0", "#1B9AAA")

p_d <- ggplot() +
  geom_ribbon(data = gini_per_tree_null_df,
              aes(x = T, ymin = gini_null_lo, ymax = gini_null_hi, fill = "Yule-null"),
              alpha = 0.2, color = NA) +
  geom_line(data = gini_per_tree_df, aes(x = T, y = gini, color = as.factor(tree)), linewidth = 0.6)+
  geom_point(data = gini_per_tree_df, aes(x = T, y = gini, color = as.factor(tree)), size = 0.6) +
  scale_color_manual(values = tree_pal) +
  scale_fill_manual(values = "grey") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Days post-fertilization", y = "Gini coefficient") +
  fig_theme

p_d

#ggsave(paste0("figures/sensitivity/", "sens_2"), p_d, height = 3.5, width = 4.5)

#recreate the analysis per blastomere, e.g. using  some code from here to assign
# a blastomere label
assign_founders <- function(tree, T0, larger_label = "B1", smaller_label = "B2") {
  nt   <- length(tree$tip.label)
  root <- nt + 1L
  kids <- tree$edge[tree$edge[, 1] == root, 2]
  stopifnot(length(kids) == 2)   # bifurcating root by construction
  
  tips_under <- function(node) {
    if (node <= nt) return(tree$tip.label[node])
    castor::get_subtree_at_node(tree, node - nt)$subtree$tip.label
  }
  child_tips  <- lapply(kids, tips_under)
  child_sizes <- vapply(child_tips, length, integer(1))
  larger      <- which.max(child_sizes)
  larger_set  <- child_tips[[larger]]
  
  split      <- castor::split_tree_at_height(tree, height = T0)
  clone_size <- vapply(split$subtrees, function(s) length(s$tree$tip.label), integer(1))
  tips_by_founder <- lapply(split$subtrees, function(s) s$tree$tip.label)
  blastomere <- vapply(split$subtrees, function(s)
    if (s$tree$tip.label[1] %in% larger_set) larger_label else smaller_label, character(1))
  stopifnot(sum(clone_size) == nt)
  
  list(blastomere = blastomere, clone_size = clone_size, tips_by_founder = tips_by_founder)
}

#-------------------------------------------------------------------------#
# Gini over time within each blastomere lineage (A/B1 = larger, B/B2 = smaller),
# for every tree
compute_gini_per_blastomere <- function(tree, T0, T_max, n_grid,
                                        larger_label = "A", smaller_label = "B") {
  founders <- assign_founders(tree, T0, larger_label = larger_label, smaller_label = smaller_label)
  N0 <- length(founders$blastomere)

  tip_founder <- integer(length(tree$tip.label))
  names(tip_founder) <- tree$tip.label
  for (i in seq_along(founders$tips_by_founder)) {
    tip_founder[founders$tips_by_founder[[i]]] <- i
  }

  T_seq <- seq(T0, T_max, length.out = n_grid)

  do.call(rbind, lapply(T_seq, function(T) {
    split_T <- castor::split_tree_at_height(tree, height = T)
    first_tip <- vapply(split_T$subtrees, function(s) s$tree$tip.label[1], character(1))
    counts <- tabulate(tip_founder[first_tip], nbins = N0)
    stopifnot(all(counts >= 1))

    do.call(rbind, lapply(c(larger_label, smaller_label), function(lbl) {
      idx <- which(founders$blastomere == lbl)
      data.frame(T = T, blastomere = lbl, gini = ineq::Gini(counts[idx]))
    }))
  }))
}

gini_blastomere_df <- do.call(rbind, lapply(names(trees), function(nm) {
  df <- compute_gini_per_blastomere(trees[[nm]], T0 = T0_values, T_max = T_max, n_grid = n_grid)
  df$tree <- nm
  df
}))

write.csv(gini_blastomere_df, paste0("figures/sensitivity/", "gini_per_blastomere.csv"), row.names = FALSE)

p_d3 <- ggplot(gini_blastomere_df, aes(x = T, y = gini, color = tree)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 0.6) +
  facet_wrap(~ blastomere, labeller = as_labeller(c(A = "Blastomere A (B1)", B = "Blastomere B (B2)"))) +
  scale_color_manual(values = tree_pal) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Days post-fertilization", y = "Gini coefficient") +
  fig_theme

p_d3

ggsave(paste0("figures/sensitivity/", "gini_per_blastomere.pdf"), p_d3, height = 3.5, width = 7)
