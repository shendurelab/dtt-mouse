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

#-------------------------------------------------------------------------#
#Panel A -- blastomere-split Gini(T), main tree vs two root-split-ratio
# sensitivity trees, one blastomere per facet
gini_blast_sens_df <- read.csv("./figures_data/fig_s14_panel_A_gini_per_blastomere_sensitivity.csv")

# named explicitly, not positional -- avoids the CSV round-trip factor-order
# issue (alphabetical re-sort on read.csv) silently remapping colors to trees
tree_pal <- c(main = "#2B8C2C", "sensitivity_58-42" = "#7B4FA0", "sensitivity_63-37" = "#1B9AAA")

p_a <- ggplot(gini_blast_sens_df, aes(x = T, y = gini, color = tree)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 0.6) +
  facet_wrap(~blast) +
  scale_color_manual(values = tree_pal, name = NULL) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Days post-fertilization", y = "Gini coefficient") +
  fig_theme

p_a

#-------------------------------------------------------------------------#
#Panel B -- blastomere-split Gini(T), main tree only, single panel
gini_emergence_main_df <- read.csv("./figures_data/fig_s14_panel_B_gini_emergence.csv")

blast_colors <- c("Blastomere A" = "#3573b9", "Blastomere B" = "#e8863b")

p_b <- ggplot(gini_emergence_main_df, aes(x = T, y = gini, color = blast)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = blast_colors, name = NULL) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Days post-fertilization", y = "Gini coefficient") +
  fig_theme

p_b

#-------------------------------------------------------------------------#
#Panel C -- Gini(T) vs Yule null, backbone tree, one line per candidate T0
t0_obs_df  <- read.csv("./figures_data/fig_s14_panel_C_gini_over_time_T0sweep.csv")
t0_null_df <- read.csv("./figures_data/fig_s14_panel_C_gini_over_time_T0sweep_null.csv")

t0_pal <- setNames(RColorBrewer::brewer.pal(length(unique(t0_obs_df$T0)), "Dark2"),
                    as.character(sort(unique(t0_obs_df$T0))))

p_c <- ggplot() +
  geom_ribbon(data = t0_null_df, aes(x = T, ymin = gini_null_lo, ymax = gini_null_hi, fill = as.factor(T0)),
              alpha = 0.2, color = NA) +
  geom_line(data = t0_obs_df, aes(x = T, y = gini, color = as.factor(T0)), linewidth = 0.6) +
  geom_point(data = t0_obs_df, aes(x = T, y = gini, color = as.factor(T0)), size = 0.6) +
  scale_color_manual(values = t0_pal, name = "T0") +
  scale_fill_manual(values = t0_pal, guide = "none") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Days post-fertilization", y = "Gini coefficient") +
  fig_theme

p_c

#-------------------------------------------------------------------------#
#Panel D -- excess of observed Gini over the null's 97.5% bound, same T0 sweep
t0_diff_df <- read.csv("./figures_data/fig_s14_panel_D_gini_diff_to_null_T0sweep.csv")

# no y-limits here (unlike panel C) -- gini_diff can go negative, source's
# reused [0,1] limit would silently clip those points
p_d <- ggplot(t0_diff_df, aes(x = T, y = gini_diff, color = as.factor(T0))) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 0.6) +
  scale_color_manual(values = t0_pal, name = "T0") +
  labs(x = "Days post-fertilization", y = "Difference in Gini to Null") +
  fig_theme + 
  ylim(0, 1.0)

p_d

#-------------------------------------------------------------------------#
#Panel E -- founder x cell-type heatmap (backbone tree), fill = log2(obs/exp),
# blastomere-of-origin strip above
g_long_bt <- read.csv("./figures_data/fig_s14_panel_E_prevalent_founders_heatmap_backbone.csv")
g_long_bt$founder <- factor(g_long_bt$founder, levels = unique(g_long_bt$founder))

# celltype's row order in the CSV cycles through all_types (alphabetical),
# NOT the count order founders were actually positioned by -- unique() on
# row order (as founder above) would silently keep that wrong order, so
# reconstruct from type_n (each type's total cell count) instead.
e_type_order  <- unique(g_long_bt$celltype[order(-g_long_bt$type_n)])
g_long_bt$celltype <- factor(g_long_bt$celltype, levels = rev(e_type_order))

# labels carry each type's total sampled cell count
e_type_totals <- setNames(g_long_bt$type_n[match(e_type_order, g_long_bt$celltype)], e_type_order)
e_ycols       <- ifelse(e_type_order %in% names(TRAJECTORY_COLORS), TRAJECTORY_COLORS[e_type_order], "black")
e_ylabs       <- setNames(sprintf("%s (n=%s)", e_type_order, label_comma()(e_type_totals)), e_type_order)

p_e_main <- ggplot(g_long_bt, aes(x = founder, y = celltype, fill = log2fc)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                        na.value = "grey90", name = "log2\n(obs/exp)") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0), labels = e_ylabs) +
  labs(x = sprintf("All founders (backbone tree, n = %d)", length(unique(g_long_bt$founder))), y = NULL) +
  theme_classic(base_size = base_size, base_family = "sans") +
  theme(
    axis.text.x     = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y     = element_text(size = base_size - 1, colour = e_ycols),
    axis.line       = element_blank(),
    axis.title.x    = element_text(size = base_size),
    legend.title    = element_text(size = base_size - 1),
    legend.text     = element_text(size = base_size - 1),
    legend.key.size = unit(3, "mm"),
    plot.margin     = margin(3, 5, 3, 3)
  )

e_strip <- unique(g_long_bt[, c("founder", "blast")])

p_e_strip <- ggplot(e_strip, aes(x = founder, y = 1, fill = blast)) +
  geom_tile() +
  scale_fill_manual(values = blast_colors, name = NULL) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(legend.text = element_text(size = base_size - 1), legend.key.size = unit(2.5, "mm"))

p_e <- p_e_strip + p_e_main + plot_layout(heights = c(0.06, 1), guides = "collect")

p_e

#-------------------------------------------------------------------------#
#Panel F -- per-cell-type Gini consistency between full and backbone trees,
# one blastomere per facet
type_consistency_df <- read.csv("./figures_data/fig_s14_panel_F_consistency_percelltype.csv")
rho_df <- read.csv("./figures_data/fig_s14_panel_F_consistency_percelltype_rho.csv")
rho_df$label <- sprintf("n=%d\nrho=%.2f", rho_df$n, rho_df$rho)

p_f <- ggplot(type_consistency_df, aes(x = gini_backbone, y = gini_full)) +
  geom_abline(slope = 1, intercept = 0, color = "grey70", linetype = "dashed", linewidth = 0.4) +
  geom_point(aes(color = major_trajectory, size = min_founders), alpha = 0.85) +
  geom_text(data = rho_df, aes(x = 0.05, y = 0.98, label = label), inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = base_size / .pt * 0.75, color = "grey30") +
  scale_color_manual(values = TRAJECTORY_COLORS, name = "major trajectory") +
  scale_size_continuous(name = "min(n founders)", range = c(1.5, 4)) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~blast) +
  guides(color = guide_legend(order = 1, override.aes = list(size = 2.4)), size = guide_legend(order = 2)) +
  labs(x = "Per-trajectory Gini on backbone tree", y = "Per-trajectory Gini on full tree") +
  fig_theme

p_f
