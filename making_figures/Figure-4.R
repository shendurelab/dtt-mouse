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

k_for_frac <- function(sizes, frac) {
  s <- sort(sizes, decreasing = TRUE)
  which(cumsum(s) >= frac * sum(s))[1]
}

# Palette: observed / data = muted steel blue, null + reference lines = neutral grey.
col_obs    <- "#4c72b0"
col_bar    <- "#8a8a8a"   # muted single-hue bars (one series, no identity to encode)
col_null   <- "grey65"    # Yule null bar = neutral reference baseline, not a series
col_cutoff <- "grey45"    # grey dashed reference marker (paper convention)


#-------------------------------------------------------------------------#
#Panel A 
ft_rank_df = read.csv("./figures_data/fig4_panel_A.csv")

k_obs_b <- k_for_frac(ft_rank_df$size, frac = 0.5)


p_a <- ggplot(ft_rank_df, aes(x = rank, y = size)) +
  geom_col(width = 1, fill = col_bar) +
  geom_vline(xintercept = k_obs_b, linetype = "dashed", color = col_cutoff, linewidth = 0.4) +
  annotate("text", x = k_obs_b, y = max(ft_rank_df$size), label = sprintf("top %d", k_obs_b),
           hjust = -0.15, vjust = 1, size = base_size / .pt * 0.85, color = col_cutoff) +
  scale_y_log10(labels = label_comma()) +
  labs(x = "Founder rank", y = "Clone size (log scale)") +
  fig_theme

p_a

#-------------------------------------------------------------------------#
#Panel B
dat_b = read.csv("./figures_data/fig4_panel_B_founders_for_half.csv")


p_b <- ggplot(data.frame(k_null = dat_b$k_null), aes(x = k_null)) +
  geom_histogram(binwidth = 1, boundary = 0.5, fill = col_null, color = "white", linewidth = 0.1) +
  geom_vline(xintercept = k_obs_b, color = col_obs, linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = k_obs_b, y = Inf, label = sprintf("observed = %d", k_obs_b),
           color = col_obs, hjust = -0.08, vjust = 1.5, size = base_size / .pt * 0.85) +
  # atop(), not embedded "\n" -- see panel a's comment.
  labs(
       x = sprintf("# founders for %.0f%%\nof cells (full tree)", 100 * 0.5), y = "count (null reps)") +
  fig_theme

p_b

#-------------------------------------------------------------------------#
#Panel C
gini_df = read.csv("./figures_data/fig4_panel_C_gini_vs_null.csv")

p_c <- ggplot(gini_df, aes(x = gini)) +
  geom_density(fill = "grey75", color = "grey55", alpha = 0.6, linewidth = 0.4) +
  geom_vline(xintercept = gini_df$gini_obs[1], color = col_obs, linetype = "dashed", linewidth = 0.5) +
  
  annotate("text", x = gini_df$gini_obs[1], y = Inf, label = sprintf("observed = %.3f", gini_df$gini_obs[1]),
           color = col_obs, hjust = -0.05, vjust = 1.5, size = base_size / .pt * 0.85) +
  scale_x_continuous(limits = c(0, 1)) +
  # atop(), not embedded "\n" -- see panel a's comment.
  labs(
       x = "Gini coefficient", y = "Density") +
  fig_theme

p_c

#-------------------------------------------------------------------------#
#Panel D

gini_null_time_df = read.csv("figures_data/fig4_panel_D_gini_over_time_null.csv")
gini_fixed_df = read.csv("figures_data/fig4_panel_D_gini_over_time.csv")

p_d <- ggplot(gini_fixed_df, aes(x = T, y = gini)) +
  geom_ribbon(data = gini_null_time_df, aes(x = T, ymin = gini_null_lo, ymax = gini_null_hi),
              inherit.aes = FALSE, fill = col_null, alpha = 0.5) +
  geom_line(data = gini_null_time_df, aes(x = T, y = gini_null_median),
            inherit.aes = FALSE, color = "grey45", linetype = "dashed", linewidth = 0.4) +
  geom_line(color = col_obs, linewidth = 0.6) +
  geom_point(color = col_obs, size = 0.6) +
  
  scale_y_continuous(limits = c(0, 1)) +
  labs(subtitle = "Yule null: median (dashed)\n+ 95% envelope",
       x = "Days post-fertilization", y = "Gini coefficient") +
  fig_theme

p_d

#-------------------------------------------------------------------------#
#Panel F

cons_df_blastnull = read.csv("./figures_data/fig4_panel_F_gini_excess.csv")
# quantile-binned diverging fill scale for diff_hi
pos_diff   <- cons_df_blastnull$diff_hi[cons_df_blastnull$diff_hi > 0]
n_pos_bins <- 4
pos_cuts   <- round(quantile(pos_diff, seq(1, n_pos_bins - 1) / n_pos_bins, na.rm = TRUE), 2)
diff_range_raw <- range(cons_df_blastnull$diff_hi, na.rm = TRUE)
diff_range  <- c(floor(diff_range_raw[1] * 100) / 100, ceiling(diff_range_raw[2] * 100) / 100)
diff_breaks <- sort(unique(c(diff_range[1], 0, pos_cuts, diff_range[2])))
zero_pos    <- scales::rescale(0, from = diff_range)

# y-axis labels carry each type's FULL-TREE cell count (same basis as panel
# F's n= labels, g_ncells_all below) even though row order/filtering above is
# backbone-tree-based -- computed independently here since panel F's block
# runs later in this script.
ft_type_n_all   <- table(factor(traj_of[ft_tree$tip.label], levels = sort(unique(na.omit(traj_of[ft_tree$tip.label])))))
cons_type_label <- setNames(sprintf("%s (n=%s)", cons_type_order, label_comma()(as.numeric(ft_type_n_all[cons_type_order]))), cons_type_order)

p_f_blastnull_heatmap <- ggplot(cons_df_blastnull, aes(x = T, y = major_trajectory, fill = diff_hi)) +
  geom_tile() +
  scale_fill_stepsn(colours = c("#2166ac", "white", "#b2182b"),
                    values = c(0, zero_pos, 1),
                    breaks = diff_breaks, limits = diff_range,
                    name = "Excess gini",
                    
                    guide = guide_coloursteps(even.steps = TRUE, show.limits = TRUE)) +
  scale_y_discrete(limits = rev(cons_type_order), labels = cons_type_label[rev(cons_type_order)]) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_wrap(~blast) +
  labs(
    x = "Days post-fertilization", y = NULL) +
  fig_theme +
  theme(axis.text.y = element_text(size = base_size - 2))

p_f_blastnull_heatmap


#-------------------------------------------------------------------------#
#Panel E founders heatmap
TRAJECTORY_COLORS <- c(
  Neuroectoderm_and_glia            = "#e16b38",
  Intermediate_neuronal_progenitors = "#434195",
  Eye_and_other                     = "#6ab65c",
  Ependymal_cells                   = "#8866a6",
  CNS_neurons                       = "#dec148",
  Mesoderm                          = "#99589c",
  Definitive_erythroid              = "#cd5246",
  Epithelium                        = "#ac9fb4",
  Endothelium                       = "#4ca257",
  Muscle_cells                      = "#d8a8c8",
  Hepatocytes                       = "#3e5f39",
  White_blood_cells                 = "#889ac9",
  Neural_crest_PNS_glia             = "#f3ec7d",
  Primitive_erythroid               = "#e5a9a1",
  Neural_crest_PNS_neurons          = "#b3c595",
  T_cells                           = "#df9c5b",
  Lung_and_airway                   = "#50abca",
  Intestine                         = "#cf3579",
  B_cells                           = "#52aea1",
  Olfactory_sensory_neurons         = "#c03c32",
  Cardiomyocytes                    = "#735898",
  Oligodendrocytes                  = "#a18c5c",
  Mast_cells                        = "#48707b",
  Megakaryocytes                    = "#9b909a",
  Testis_and_adrenal                = "#7b7d65"
)

g_long = read.csv("figures_data/fig4_panel_E_prevalent_founders_heatmap.csv")

# CSV drops factor levels -> rebuild write-time order from columns already in the file
founder_levels  <- unique(g_long$founder[order(as.integer(sub("^f", "", g_long$founder)))])
celltype_levels <- unique(g_long$celltype[order(g_long$type_n)])   # ascending n -> smallest at axis bottom

g_long$founder  <- factor(g_long$founder,  levels = founder_levels)
g_long$celltype <- factor(g_long$celltype, levels = celltype_levels)

# y-axis labels carry each type's total sampled cell count
g_ylabs      <- rev(g_type_order)
g_ycols      <- ifelse(g_ylabs %in% names(TRAJECTORY_COLORS), TRAJECTORY_COLORS[g_ylabs], "black")
g_type_label <- setNames(sprintf("%s (n=%s)", g_type_order, label_comma()(type_total[g_type_order])), g_type_order)

p_g_main <- ggplot(g_long, aes(x = founder, y = celltype, fill = log2fc)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                       na.value = "grey90", name = "log2\n(obs/exp)") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0), labels = g_type_label[rev(g_type_order)]) +
  labs(x = sprintf("All founders (full tree, n = %s)", format(nrow(M_all), big.mark = ",")), y = NULL) +
  theme_classic(base_size = base_size, base_family = "sans") +
  theme(
    axis.text.x     = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y     = element_text(size = base_size - 1, colour = g_ycols),
    axis.line       = element_blank(),
    axis.title.x    = element_text(size = base_size),
    legend.title    = element_text(size = base_size - 1),
    legend.text     = element_text(size = base_size - 1),
    legend.key.size = unit(3, "mm"),
    plot.margin     = margin(3, 5, 3, 3)
  )

# blastomere-of-origin strip, aligned to p_g_main's founder order
g_strip <- unique(g_long[, c("founder", "blast")])
g_strip$blast <- to_disp(g_strip$blast)

p_g_strip <- ggplot(g_strip, aes(x = founder, y = 1, fill = blast)) +
  geom_tile() +
  scale_fill_manual(values = blast_colors, name = NULL) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(legend.text = element_text(size = base_size - 1), legend.key.size = unit(2.5, "mm"))

p_g <- p_g_strip + p_g_main + plot_layout(heights = c(0.06, 1), guides = "collect")

p_g

