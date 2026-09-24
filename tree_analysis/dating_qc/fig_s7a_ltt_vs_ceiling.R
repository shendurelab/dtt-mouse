#!/usr/bin/env Rscript
# fig_s7a_ltt_vs_ceiling.R
# Fig. S7A: lineages-through-time (LTT) of the constrained vs. unconstrained
# dated trees, against the literature embryo cell-count ceiling. 
#
# Run from the repo root:
#   Rscript tree_analysis/dating_qc/fig_s7a_ltt_vs_ceiling.R

suppressMessages({ library(ape) })

TREE_CONSTRAINED   <- "tree_building/results/3-dated-tree/merged_minB2h_lineage_constrained.nwk"
TREE_UNCONSTRAINED <- "tree_building/results/3-dated-tree/merged_minB2h_unconstrained.nwk"
COUNTS_CSV         <- "tree_building/processed_data/cell_counts_full_kojima_cao.csv"
stopifnot(file.exists(TREE_CONSTRAINED), file.exists(TREE_UNCONSTRAINED), file.exists(COUNTS_CSV))
OUT_CSV <- "figures_data/fig_s7_dating_panel_B_ltt_vs_ceiling.csv"

# Lineages-through-time: number of tree edges crossing each day in `days`.
# An edge (parent age pa, child age ch) crosses day d iff pa < d <= ch
# (half-open: a split exactly at d has already happened by d); a query at/after
# the tree height is clamped just inside it so the tip day returns the full
# tip count instead of 0 (tips sit a hair below the nominal height in fp).
ltt_at <- function(tr, days) {
  ag <- ape::node.depth.edgelength(tr)
  pa <- ag[tr$edge[, 1]]; ch <- ag[tr$edge[, 2]]
  ht <- max(ag)
  vapply(days, function(d) {
    dd <- min(d, ht - 1e-3)
    sum(pa < dd & ch >= dd)
  }, numeric(1))
}

# ---- 1. grid & LTT curves --------------------------------------------------
grid <- seq(1.5, 13.45, 0.02)   # stop at 13.45 to avoid a spurious count drop at the tips
ltt_b <- ltt_at(read.tree(TREE_UNCONSTRAINED), grid)
ltt_g <- ltt_at(read.tree(TREE_CONSTRAINED),   grid)
ltt_b[ltt_b < 1] <- NA
ltt_g[ltt_g < 1] <- NA

# ---- 2. embryo cell-count ceiling ------------------------------------------
cnt      <- read.csv(COUNTS_CSV)                          # day, cells_lo, cells_hi, basis
ceil_all <- data.frame(day = cnt$day, cells = cnt$cells_hi)
ceil_all <- ceil_all[order(ceil_all$day), ]
ceil_line <- exp(approx(ceil_all$day, log(ceil_all$cells), grid)$y)

# ---- 3. tidy frames + ONE unified series key -------------------------------
# Ceiling is split by whether the counts were used to build the min-age
# constraint (<=E7.5) or not (>E7.5, context only).
lv <- c("unconstrained dating", "constrained dating",
        "cell ceiling, used (≤ E7.5)", "cell ceiling, unused (> E7.5)")

df_ltt <- rbind(
  data.frame(day = grid, count = ltt_b, series = lv[1]),
  data.frame(day = grid, count = ltt_g, series = lv[2]))

df_ceil      <- data.frame(day = grid, count = ceil_line)
df_ceil_used <- transform(df_ceil[df_ceil$day <= 7.5, ], series = lv[3])
df_ceil_free <- transform(df_ceil[df_ceil$day >= 7.5, ], series = lv[4])  # shares E7.5 join
for (d in c("df_ltt","df_ceil_used","df_ceil_free"))
  assign(d, within(get(d), series <- factor(series, levels = lv)))

# landmark points (the actual literature day/count rows within the grid range)
df_land      <- ceil_all[ceil_all$day >= min(grid) & ceil_all$day <= max(grid), ]
df_land_used <- df_land[df_land$day <= 7.5, ]
df_land_free <- df_land[df_land$day >  7.5, ]

# ---- 4. tidy CSV export -----------------------------------------------------
# `series` matches the plot's legend text; landmark points share the used/free
# ceiling series of the line under them (kind distinguishes line vs landmark).
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

write.csv(tidy_long, OUT_CSV, row.names = FALSE)
cat("wrote", OUT_CSV, "\n")

# ---- 5. sanity numbers ------------------------------------------------------
i65 <- which.min(abs(grid - 6.5))
cat(sprintf("At E6.5: unconstrained LTT=%g constrained LTT=%g ceiling=%g (x%.1f over)\n",
            ltt_b[i65], ltt_g[i65], ceil_line[i65], ltt_b[i65] / ceil_line[i65]))
