#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# fig_ltt_paper.R
#
# Lineages-through-time (LTT) of two
# dated trees against the mouse whole-embryo cell-count "ceiling".
#
# Our constrained dating (green) keeps the sampled lineage count at or
# below the physically-possible embryo cell number, removing the ~60x early
# overshoot of the unconstrained dating (blue).
#
# ---------------------------------------------------------------------------

suppressMessages({
  library(ape)
  library(ggplot2)
  library(scales)
})


## ---- 0a. font handling -----------------------------------------------
# Match the paper's clean sans face (Helvetica); render through ragg/cairo so
# it is embedded, degrading gracefully to "sans" + default devices otherwise.
font_family  <- "sans"
has_ragg     <- requireNamespace("ragg", quietly = TRUE)
has_sysfonts <- requireNamespace("systemfonts", quietly = TRUE)
if (has_sysfonts) {
  hv <- tryCatch(systemfonts::match_fonts("Helvetica"), error = function(e) NULL)
  if (!is.null(hv) && nrow(hv) && nzchar(hv$path[1]) && file.exists(hv$path[1]))
    font_family <- "Helvetica"
}
has_cairo_pdf <- capabilities("cairo")

# target point size -> mm units used by annotate("text", ...)
pt <- function(x) x / ggplot2::.pt

## ---- 0b. paths ------------------------------------------------------------
tree_constrained_main <- "tree_building/results/3-dated-tree/merged_minB2h_lineage_constrained.nwk"
tree_constrained_sens <- "tree_building/results/3-dated-tree-sensitivity/merged_minB2h_lineage_constrained_sidefrac63-37.nwk"
tree_unconstrained <- "tree_building/results/3-dated-tree/merged_minB2h_unconstrained.nwk"

PANEL_TAG          <- "B"
fig_stem           <- paste0("panel_", PANEL_TAG, "_ltt_vs_ceiling") 

# Single cell-count source for the whole ceiling: the CX (Kojima/Cao) table,
# E0.5..E13.5 (day, cells_lo, cells_hi). Nothing here mixes in Qiu counts.
counts_full        <- "tree_analysis/supporting_data/cell_counts_full_kojima_cao.csv"


## ---- 1. LTT helper: ltt_at(tree, days)  -----
## (trees ultrametric, root = day-0 zygote, so node depth = embryonic day)
ltt_at <- function(tr, days) {
  ag <- ape::node.depth.edgelength(tr)
  pa <- ag[tr$edge[, 1]]; ch <- ag[tr$edge[, 2]]
  ht <- max(ag)
  vapply(days, function(d) {
    dd <- min(d, ht - 1e-3)
    sum(pa < dd & ch >= dd)
  }, numeric(1))
}

## ---- 2. grid & LTT curves for every tree ------------------------------------------------
grid <- seq(1.5, 13.5, 0.02)   
ltt_b <- ltt_at(read.tree(tree_unconstrained), grid)
ltt_g <- ltt_at(read.tree(tree_constrained_main),   grid)
ltt_s <- ltt_at(read.tree(tree_constrained_sens),   grid)

ltt_b[ltt_b < 1] <- NA
ltt_g[ltt_g < 1] <- NA
ltt_s[ltt_s < 1] <- NA

## ---- 3. embryo cell-count ceiling that was used for constrained dating --------
cnt      <- read.csv(counts_full)                       # day, cells_lo, cells_hi
ceil_all <- data.frame(day = cnt$day, cells = cnt$cells_hi)
ceil_all <- ceil_all[order(ceil_all$day), ]
ceil_line <- exp(approx(ceil_all$day, log(ceil_all$cells), grid)$y)

## ---- 4. tidy frames + ONE unified series key ----------------------------
# All four elements share a single colour+linetype scale so they render as one
# legend. Ceiling is split by whether the counts were used to
# build the min-age constraint (<=E7.5) or not (>E7.5, context only).
lv <- c("unconstrained dating", "constrained dating 50/50", "constrained dating 63/37",
        "cell ceiling, used (≤ E7.5)", "cell ceiling, unused (> E7.5)")

df_ltt <- rbind(
  data.frame(day = grid, count = ltt_b, series = lv[1]),
  data.frame(day = grid, count = ltt_g, series = lv[2]),
  data.frame(day = grid, count = ltt_s, series = lv[3])
  )

df_ceil      <- data.frame(day = grid, count = ceil_line)
df_ceil_used <- transform(df_ceil[df_ceil$day <= 7.5, ], series = lv[4])
df_ceil_free <- transform(df_ceil[df_ceil$day >= 7.5, ], series = lv[5])  # shares E7.5 join
for (d in c("df_ltt","df_ceil_used","df_ceil_free"))
  assign(d, within(get(d), series <- factor(series, levels = lv)))

# landmark dots (no legend entry), coloured by the same split
df_land      <- ceil_all[ceil_all$day >= min(grid) & ceil_all$day <= max(grid), ]
df_land_used <- df_land[df_land$day <= 7.5, ]
df_land_free <- df_land[df_land$day >  7.5, ]

## ---- 5. limits, colours, theme ------------------------------------------
y_top <- 3e7; y_bot <- 8; x_lim <- c(1.5, 13.5)   # ceiling now tops out at ceiling max 13M (E13.5)
col_uncon <- "#2a78d6"; col_con5050 <- "#008300"; col_con6347 <- "#F4C430"; col_ceil <- "#d1495b"; col_free <- "#898781"
pal_col <- c(
  "unconstrained dating"          = col_uncon,
  "constrained dating 50/50"      = col_con5050,
  "constrained dating 63/37"      = col_con6347,
  "cell ceiling, used (≤ E7.5)"   = col_ceil,
  "cell ceiling, unused (> E7.5)" = col_free
)
pal_lty <- setNames(c("solid", "solid", "solid", "solid", "solid"), lv)  # all series drawn as normal solid lines

theme_paper <- theme_classic(base_size = 8) + theme(
  text = element_text(family = font_family, colour = "black"),
  legend.text = element_text(size = 7, margin = margin(l = 2)),
  legend.title = element_blank(),
  legend.background = element_blank(), legend.key = element_blank(),
  plot.margin = margin(t = 6, r = 8, b = 3, l = 3),
  plot.title = element_blank())

lab_end <- df_land[which.max(df_land$day), ]   # E13.5 endpoint (max cell count, CX)

## ---- 6. build panel ------------------------------------------------------
p <- ggplot() +
  # ceiling: grey (unused) first, coral (used) on top so the E7.5 join is clean
  geom_line(data = df_ceil_used, aes(day, count, linetype = series, colour = series), linewidth = 0.7) +
  geom_line(data = df_ceil_free, aes(day, count, linetype = series, colour = series), linewidth = 0.7) +
  
  #geom_line(data = df_ceil_free, aes(day, count), colour = "#898781", linetype = series), linewidth = 0.7) +
  geom_point(data = df_land_free, aes(day, cells), colour = col_free, size = 1.4, show.legend = FALSE) +
  #geom_line(data = df_ceil_used, aes(day, count, colour = "#d1495b", linetype = series), linewidth = 0.7) +
  geom_point(data = df_land_used, aes(day, cells), colour = "#d1495b", size = 1.4, show.legend = FALSE) +
  # LTT curves
  geom_line(data = df_ltt, aes(day, count, colour = series, linetype = series), linewidth = 0.7, na.rm = TRUE) +
  scale_colour_manual(values = pal_col, breaks = lv, name = NULL) +
  scale_linetype_manual(values = pal_lty, breaks = lv, name = NULL) +
  scale_y_log10(limits = c(y_bot, y_top),
                breaks = 10^(1:8),
                labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                expand = c(0, 0)) +
  scale_x_continuous(breaks = c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5),
                     limits = x_lim, expand = c(0.01, 0)) +
  labs(x = "Embryonic day (E)", y = "Number of lineage / embryo cells", tag = PANEL_TAG) +
  theme_paper +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.985, 0.02),
    legend.justification = c(1, 0),
    legend.spacing.y = grid::unit(1, "pt"),
    legend.key.spacing.y = grid::unit(1, "pt"),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.width = grid::unit(13, "pt"),
    legend.key.height = grid::unit(10, "pt"),
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0.01, 0.985)) +
  guides(colour = guide_legend(byrow = TRUE), linetype = guide_legend(byrow = TRUE))

p

## ---- 7b. tidy CSV export --------------------------------------------------
# Row-bind the five dataframes layered into the plot (step 6) into one long
# table so every plotted point is available outside R. `series` is copied
# verbatim from the plot's own legend text (`lv`); the landmark points don't
# get their own legend entry (show.legend=FALSE) but sit on -- and belong to --
# the same used/free ceiling series as the line under them, so they reuse that
# line's label. `kind` (line vs landmark) keeps the two visually distinct even
# though they can share a `series` value.
line_ltt   <- transform(df_ltt,        kind = "line")
line_ceilU <- transform(df_ceil_used,  kind = "line")
line_ceilF <- transform(df_ceil_free,  kind = "line")
land_used  <- data.frame(day = df_land_used$day, count = df_land_used$cells,
                          series = lv[3], kind = "landmark")
land_free  <- data.frame(day = df_land_free$day, count = df_land_free$cells,
                          series = lv[4], kind = "landmark")
tidy_cols <- c("day", "count", "series", "kind")
tidy_long <- rbind(line_ltt[, tidy_cols], line_ceilU[, tidy_cols],
                    line_ceilF[, tidy_cols], land_used[, tidy_cols],
                    land_free[, tidy_cols])
tidy_long$series <- as.character(tidy_long$series)

# Can write csv if helpful
#csv_path <- file.path(out_dir, paste0(fig_stem, ".csv"))
#write.csv(tidy_long, csv_path, row.names = FALSE)

## ---- 8. sanity-check numbers ---------------------------------------------------
i65 <- which.min(abs(grid - 6.5))
cat(sprintf("At E6.5: unconstrained LTT=%g constrained LTT=%g ceiling=%g (x%.1f over)\n",
            ltt_b[i65], ltt_g[i65], ceil_line[i65], ltt_b[i65] / ceil_line[i65]))
