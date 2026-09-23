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

# how many founders are needed to get to 50% of cells in the tree 
FRAC_B = 0.5

# Palette: observed / data = muted steel blue, null + reference lines = neutral grey.
col_obs    <- "#4c72b0"
col_bar    <- "#8a8a8a"   # muted single-hue bars (one series, no identity to encode)
col_null   <- "grey65"    # Yule null bar = neutral reference baseline, not a series
col_cutoff <- "grey45"    # grey dashed reference marker (paper convention)


#-------------------------------------------------------------------------#
#Panel A 
bt_rank_df = read.csv("./figures_data/fig_s13_panel_A_bt.csv")


k_obs_bt <- k_for_frac(bt_rank_df$size, frac = 0.5)
k_obs_bt

## plus panel for backbone tree
p_as <- ggplot(bt_rank_df, aes(x = rank, y = size)) +
  geom_col(width = 1, fill = col_bar) +
  geom_vline(xintercept = k_obs_bt, linetype = "dashed", color = col_cutoff, linewidth = 0.4) +
  annotate("text", x = k_obs_b_bt, y = max(bt_rank_df$size), label = sprintf("top %d", k_obs_b_bt),
           hjust = -0.15, vjust = 1, size = base_size / .pt * 0.85, color = col_cutoff) +
  scale_y_log10(labels = label_comma()) +
  labs(
       x = "Founder rank", y = "Clone size (log scale)") +
  fig_theme

p_as

#-------------------------------------------------------------------------#
#Panel B
dat_b = read.csv("./figures_data/fig_s13_panel_B_founders_for_half.csv")

p_b2 <- ggplot(data.frame(k_null = dat_b$k_null), aes(x = k_null)) +
  geom_histogram(binwidth = 1, boundary = 0.5, fill = col_null, color = "white", linewidth = 0.1) +
  geom_vline(xintercept = k_obs_b_bt, color = col_obs, linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = k_obs_b_bt, y = Inf, label = sprintf("observed = %d", k_obs_b_bt),
           color = col_obs, hjust = -0.08, vjust = 1.5, size = base_size / .pt * 0.85) +
  # atop(), not embedded "\n" -- see panel a's comment.
  labs(
       x = sprintf("# founders for %.0f%%\nof cells (backbone tree)", 100 * FRAC_B), y = "count (null reps)") +
  fig_theme

p_b2

#-------------------------------------------------------------------------#
#Panel C
gini_df_bt = read.csv("./figures_data/fig_s13_gini_vs_null.csv")


p_c2 <- ggplot(gini_df_bt, aes(x = gini)) +
  geom_density(fill = "grey75", color = "grey55", alpha = 0.6, linewidth = 0.4) +
  geom_vline(xintercept = gini_df_bt$gini_obs[1], color = col_obs, linetype = "dashed", linewidth = 0.5) +
  
  annotate("text", x = gini_df_bt$gini_obs[1], y = Inf, label = sprintf("observed = %.3f", gini_df_bt$gini_obs[1]),
           color = col_obs, hjust = 0.6, vjust = 1.5, size = base_size / .pt * 0.85) +
  scale_x_continuous(limits = c(0, 1)) +
  # atop(), not embedded "\n" -- see panel a's comment.
  labs(
    x = "Gini coefficient", y = "Density") +
  fig_theme

p_c2

#-------------------------------------------------------------------------#
#Panel D -- blastomere-split rank-abundance (full tree)
d_rank_df   <- read.csv("./figures_data/fig_s13_panel_D_rank_abundance.csv")
d_cutoff_df <- read.csv("./figures_data/fig_s13_panel_D_rank_abundance_cutoff.csv")

# geom_text(data=...), not annotate() -- each facet (blastomere) needs its
# own label/position; annotate() would draw the same one on every facet.
p_d <- ggplot(d_rank_df, aes(x = rank, y = size)) +
  geom_col(width = 1, fill = col_bar) +
  geom_vline(data = d_cutoff_df, aes(xintercept = k), linetype = "dashed", color = col_cutoff, linewidth = 0.4) +
  geom_text(data = d_cutoff_df, aes(x = k, y = max(d_rank_df$size), label = sprintf("Top %d", k)),
            inherit.aes = FALSE, hjust = -0.15, vjust = 1, size = base_size / .pt * 0.85, color = col_cutoff) +
  scale_y_log10(labels = label_comma()) +
  facet_wrap(~blast, scales = "free_x") +
  labs(x = "Founder rank", y = "Clone size (log scale)") +
  fig_theme

p_d

#-------------------------------------------------------------------------#
#Panel E -- blastomere-split founders-for-half vs Yule null (full tree)
e_null_df <- read.csv("./figures_data/fig_s13_panel_E_founders_for_half_null.csv")
e_obs_df  <- read.csv("./figures_data/fig_s13_panel_E_founders_for_half_obs.csv")
e_obs_df$label <- sprintf("observed = %d\np %s", e_obs_df$k_obs,
                           ifelse(e_obs_df$p_val == 0, "< 1e-04", sprintf("= %.2g", e_obs_df$p_val)))

p_e <- ggplot(e_null_df, aes(x = k_null)) +
  geom_histogram(binwidth = 1, boundary = 0.5, fill = col_null, color = "white", linewidth = 0.1) +
  geom_vline(data = e_obs_df, aes(xintercept = k_obs), color = col_obs, linetype = "dashed", linewidth = 0.5) +
  geom_text(data = e_obs_df, aes(x = k_obs, y = Inf, label = label),
            inherit.aes = FALSE, color = col_obs, hjust = -0.08, vjust = 1.5, size = base_size / .pt * 0.85) +
  facet_wrap(~blast) +
  labs(x = "# founders needed for 50%\nof E13.5 cells in full tree", y = "Count (null reps)") +
  fig_theme

p_e

#-------------------------------------------------------------------------#
#Panel F -- blastomere-split Gini vs Yule null (backbone tree)
f_null_df <- read.csv("./figures_data/fig_s13_panel_F_gini_vs_null_null.csv")
f_obs_df  <- read.csv("./figures_data/fig_s13_panel_F_gini_vs_null_obs.csv")
f_obs_df$label <- sprintf("obs = %.3f\nz = %.1f\np %s", f_obs_df$gini_obs, f_obs_df$gini_z,
                           ifelse(f_obs_df$p_gini == 0, "< 1e-04", sprintf("= %.2g", f_obs_df$p_gini)))

p_f <- ggplot(f_null_df, aes(x = gini)) +
  geom_density(fill = "grey75", color = "grey55", alpha = 0.6, linewidth = 0.4) +
  geom_vline(data = f_obs_df, aes(xintercept = gini_obs), color = col_obs, linetype = "dashed", linewidth = 0.5) +
  geom_text(data = f_obs_df, aes(x = gini_obs, y = Inf, label = label),
            inherit.aes = FALSE, color = col_obs, hjust = 1.1, vjust = 1.5, size = base_size / .pt * 0.85) +
  facet_wrap(~blast) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = "Gini coefficient", y = "Null density") +
  fig_theme

p_f
