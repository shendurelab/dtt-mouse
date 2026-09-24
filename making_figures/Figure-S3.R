
###############################################################################
### Figure S3. Per-integration bulk lineage-genotype trees (mini trees)
### One midpoint-rooted NJ tree per TAPE integration (11 for embryo 3), each
### paired with a tip-aligned genotype/tape heatmap. Reads the trees
### tree_analysis/mini_trees/fig_s3_build_mini_trees.R writes; no tree-building
### here.

suppressMessages({
  library(ape); library(phangorn)
  library(ggtree); library(ggplot2); library(patchwork); library(tidyr); library(png)
})
source("tree_analysis/mini_trees/mini_tree_lib.R")
source("tree_analysis/mini_trees/mini_tree_colors.R")

GENOTYPES_CSV <- "bulk_tape/tables/DTTz_3_S3.bulk_lineage_genotypes.csv"
TREES_DIR <- "tree_analysis/mini_trees/results/trees"
PANEL_DIR <- "tree_analysis/mini_trees/results/panels"
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

INK       <- "#0b0b0b"    # primary ink -- branches, panel titles
MUTED_INK <- "#898781"    # muted ink -- axis text
SURFACE   <- "white"      # figure background
DPI       <- 300

# ---- 1. load genotypes, split by tapebc, one fixed color map for all panels -
genotypes <- read.csv(GENOTYPES_CSV, stringsAsFactors = FALSE)
color_map <- build_symbol_color_map(genotypes)   # same symbol -> hue in every panel

by_tapebc <- split(genotypes, genotypes$tapebc)
by_tapebc <- by_tapebc[order(-lengths(lapply(by_tapebc, `[[`, "genotype_id")))]

##########################################
### Fig. S3 A-K: per-integration tree + tip-aligned tape heatmap

# panels are pre-rendered to PNG (not just kept as live ggplot objects) because
# patchwork silently corrupts tile fills once more than ~6 panels are nested
# together in one composition -- a patchwork/ggplot composition limit (verified
# in vector PDF output too, not a raster artifact), so step 2 below composites
# from these images instead.
build_panel <- function(tapebc, df) {
  tree <- ape::read.tree(file.path(TREES_DIR, paste0(tapebc, ".nwk")))

  p_tree <- ggtree(tree, color = INK, linewidth = 0.3) +
    geom_tippoint(size = 0.5, color = INK) +
    ggtitle(sprintf("%s (n = %d)", tapebc, nrow(df))) +
    theme(plot.title = element_text(size = 8, face = "bold", color = INK),
          plot.background = element_rect(fill = SURFACE, color = NA))
  tip_y <- p_tree$data[p_tree$data$isTip, c("label", "y")]
  long <- build_tape_long_df(df, tip_y)

  p_tape <- ggplot(long, aes(x = site, y = y, fill = symbol)) +
    geom_tile(color = SURFACE, linewidth = 0.4) +
    scale_fill_manual(values = color_map, guide = "none") +
    scale_y_continuous(limits = range(p_tree$data$y) + c(-0.5, 0.5), expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(axis.text.x = element_text(size = 5.5, color = MUTED_INK, margin = margin(t = 1)),
          plot.background = element_rect(fill = SURFACE, color = NA))

  p_tree + p_tape + plot_layout(widths = c(1.5, 1))
}

panels <- Map(build_panel, names(by_tapebc), by_tapebc)
panel_paths <- character(length(panels))
for (i in seq_along(panels)) {
  tapebc <- names(by_tapebc)[i]
  n_tips <- nrow(by_tapebc[[i]])
  h_mm <- pmin(260, 30 + n_tips * 1.1)   # panel height scales with tip count
  panel_paths[i] <- file.path(PANEL_DIR, paste0(tapebc, ".png"))
  ggsave(panel_paths[i], panels[[i]], width = 100, height = h_mm, units = "mm", dpi = DPI)
}

##########################################
### Fig. S3 legend: symbol -> color, shared by every panel above

named <- names(SYMBOL_COLORS)[names(SYMBOL_COLORS) %in% names(color_map)]
has_other <- OTHER_COLOR %in% color_map[setdiff(names(color_map), c(named, "U"))]
legend_df <- data.frame(
  label = c(named, if (has_other) "other", "unedited"),
  fill  = c(named, if (has_other) "other", "U")
)
legend_df$label <- factor(legend_df$label, levels = legend_df$label)
legend_colors <- setNames(
  c(SYMBOL_COLORS[named], if (has_other) c(other = OTHER_COLOR), U = UNEDITED_COLOR),
  c(named, if (has_other) "other", "U")
)
p_legend <- ggplot(legend_df, aes(x = label, y = 1, fill = fill)) +
  geom_tile(color = SURFACE, linewidth = 0.4, width = 0.9, height = 0.9) +
  geom_text(aes(label = label), y = 0.35, size = 2.4, color = MUTED_INK) +
  scale_fill_manual(values = legend_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 1.6)) +
  coord_cartesian(clip = "off") +
  theme_void() + theme(plot.margin = margin(2, 2, 2, 2))
legend_path <- file.path(PANEL_DIR, "symbol_legend.png")
legend_w_mm <- 14 * nrow(legend_df)
ggsave(legend_path, p_legend, width = legend_w_mm, height = 22, units = "mm", dpi = DPI)

p_legend

##########################################
### Fig. S3: combined grid of all 11 panels

as_panel_grob <- function(path) grid::rasterGrob(png::readPNG(path), interpolate = TRUE)
p <- wrap_plots(lapply(panel_paths, as_panel_grob), ncol = 2) /
  wrap_elements(as_panel_grob(legend_path)) +
  plot_layout(heights = c(30, 2))
p
