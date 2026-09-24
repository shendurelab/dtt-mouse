# mini_tree_colors.R
# One consistent color mapping for edit-symbol -> hue, shared across all 11
# mini-tree panels (a symbol like "AAG" is written by the same TAPE recorder
# vocabulary regardless of which genomic integration/tapebc it lands on, so
# assigning it one fixed color across every panel -- not re-deriving colors
# per panel -- is what makes the symbol identity comparable across the figure).
#
# Colors match the existing manuscript figure's per-site-symbol legend
# (ColorBrewer "Paired" 12-class palette). Confirmed that legend's 12 named
# symbols are exactly 12 of our 13 observed edit symbols; the 13th (CTG,
# ~7.3% of edits -- too common to fold into "other") gets one added hue not
# reused from the 12. Note: "Paired" was not designed to the dataviz skill's
# CVD/contrast validator (several pastels fail it, e.g. light-purple vs
# light-blue under deuteranopia) -- kept anyway for consistency with the
# existing figure, per Sophie's call.
SYMBOL_COLORS <- c(
  GCC  = "#a6cee3",
  AAG  = "#1f78b4",
  CAC  = "#b2df8a",
  GATG = "#33a02c",
  CCC  = "#fb9a99",
  GGC  = "#e31a1c",
  ACG  = "#fdbf6f",
  GAA  = "#ff7f00",
  ACC  = "#cab2d6",
  ACA  = "#6a3d9a",
  CCG  = "#ffff99",
  ACT  = "#b15928",
  CTG  = "#e7298a"    # not in the reference legend -- added, not reused
)
OTHER_COLOR    <- "#999999"   # any observed symbol outside SYMBOL_COLORS (none currently)
UNEDITED_COLOR <- "#e1e0d9"   # 'U' (unedited site) -- an absence, not a data color

source("tree_analysis/mini_trees/mini_tree_lib.R")   # split_genotype_sites()

# Build the symbol -> color map from every genotype in the full (all-tapebc)
# table. Any symbol not in SYMBOL_COLORS folds into OTHER_COLOR (the "a 9th
# series folds into Other" rule) -- a no-op on the current data, since all 13
# observed symbols are named above, but keeps future/other data safe.
build_symbol_color_map <- function(all_genotypes_df) {
  sites <- unlist(lapply(all_genotypes_df$genotype, split_genotype_sites))
  symbols <- setdiff(unique(sites), "U")
  named <- intersect(symbols, names(SYMBOL_COLORS))
  other <- setdiff(symbols, names(SYMBOL_COLORS))
  cmap <- SYMBOL_COLORS[named]
  if (length(other) > 0) cmap <- c(cmap, setNames(rep(OTHER_COLOR, length(other)), other))
  c(cmap, U = UNEDITED_COLOR)
}

# Text color (white or ink) to overlay on a given fill color, chosen by
# relative luminance so an inline label stays legible on both light and dark
# tile fills -- the one case marks-and-anatomy.md allows text inside a
# colored fill. Vectorized over `hex`.
text_color_for_fill <- function(hex) {
  rgb <- grDevices::col2rgb(hex) / 255
  luminance <- 0.299 * rgb["red", ] + 0.587 * rgb["green", ] + 0.114 * rgb["blue", ]
  ifelse(luminance > 0.6, "#0b0b0b", "white")
}
