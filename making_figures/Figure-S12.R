#####################################################################
### Figure S12. Held-out-cell placement accuracy.

##################################################################################
### Fig. S12 A-D: placement distance (edges) and time-to-MRCA (days) vs. the best
### available anchor, stratified by retained tape number (A, B), total edit depth
### (C), and cell type (D). Reads the scored table
### tree_analysis/placement_accuracy/fig_s12_placement_accuracy_analysis.R writes;
### does not rerun placement.

library(data.table)
library(ggplot2)
library(cowplot)

scored <- fread("./figures_data/fig_s12_placement_accuracy_scored_dtt.csv.gz")

# every query's true anchor survives in this dataset (no query's whole side of
# the tree got dropped), so "placement" below is never split by anchor-lost
stopifnot(all(scored$true_anchor_survived))


# --- Panel A: median edge distance, by retained TAPE number ---
## A small fraction of the scored table can have NA node_dist (placement
## algorithm did not place these cells), excluded here as in the scored table.
placement_a <- scored[!is.na(node_dist), .(y = as.double(median(node_dist)),
                                           lo = quantile(node_dist, 0.25),
                                           hi = quantile(node_dist, 0.75),
                                           which = "placement"),
                      by = .(group = retained_tapes)]
anchor_a <- scored[!is.na(true_anchor_dist), .(y = as.double(median(true_anchor_dist)),
                                               lo = quantile(true_anchor_dist, 0.25),
                                               hi = quantile(true_anchor_dist, 0.75),
                                               which = "best available\nanchor in backbone"),
                   by = .(group = retained_tapes)]
panel_a <- rbind(placement_a, anchor_a)
panel_a[, which := factor(which, levels = c("placement", "best available\nanchor in backbone"))]

p_A <- ggplot(panel_a, aes(factor(group), y, fill = which)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(0.7), width = 0.2) +
  labs(x = "retained TAPE number", y = "median distance (edges)", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")
p_A

# --- Panel B: median time to MRCA (days), by retained TAPE number ---
# tree is ultrametric, so half the round-trip divergence time is the true time
# back to the MRCA from either lineage

placement_b <- scored[!is.na(mrca_divergence_days), .(y = as.double(median(mrca_divergence_days / 2)),
                                                       lo = quantile(mrca_divergence_days / 2, 0.25),
                                                       hi = quantile(mrca_divergence_days / 2, 0.75),
                                                       which = "placement"),
                      by = .(group = retained_tapes)]
anchor_b <- scored[!is.na(true_anchor_divergence_days), .(y = as.double(median(true_anchor_divergence_days / 2)),
                                                          lo = quantile(true_anchor_divergence_days / 2, 0.25),
                                                          hi = quantile(true_anchor_divergence_days / 2, 0.75),
                                                          which = "best available\nanchor in backbone"),
                   by = .(group = retained_tapes)]
panel_b <- rbind(placement_b, anchor_b)
panel_b[, which := factor(which, levels = c("placement", "best available\nanchor in backbone"))]

p_B <- ggplot(panel_b, aes(factor(group), y, fill = which)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(0.7), width = 0.2) +
  labs(x = "retained TAPE number", y = "median time to MRCA (days)", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")

p_B
# --- Panel C: median edge distance, by total edits per query (binned by 5) ---

scored[, edit_depth_bin := 5 * (total_edits %/% 5)]

placement_c <- scored[!is.na(node_dist), .(y = as.double(median(node_dist)),
                                           lo = quantile(node_dist, 0.25),
                                           hi = quantile(node_dist, 0.75),
                                           which = "placement"),
                      by = .(group = edit_depth_bin)]
anchor_c <- scored[!is.na(true_anchor_dist), .(y = as.double(median(true_anchor_dist)),
                                               lo = quantile(true_anchor_dist, 0.25),
                                               hi = quantile(true_anchor_dist, 0.75),
                                               which = "best available\nanchor in backbone"),
                   by = .(group = edit_depth_bin)]
panel_c <- rbind(placement_c, anchor_c)
panel_c[, which := factor(which, levels = c("placement", "best available\nanchor in backbone"))]

p_C <- ggplot(panel_c, aes(factor(group), y, fill = which)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(0.7), width = 0.2) +
  labs(x = "total edits per query (binned)", y = "median distance (edges)", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")


p_C
# --- Panel D: median edge distance, by major trajectory ---

placement_d <- scored[!is.na(node_dist), .(y = as.double(median(node_dist)),
                                           lo = quantile(node_dist, 0.25),
                                           hi = quantile(node_dist, 0.75),
                                           which = "placement"),
                      by = .(group = major_trajectory)]
anchor_d <- scored[!is.na(true_anchor_dist), .(y = as.double(median(true_anchor_dist)),
                                               lo = quantile(true_anchor_dist, 0.25),
                                               hi = quantile(true_anchor_dist, 0.75),
                                               which = "best available\nanchor in backbone"),
                   by = .(group = major_trajectory)]
panel_d <- rbind(placement_d, anchor_d)
panel_d[, which := factor(which, levels = c("placement", "best available\nanchor in backbone"))]

p_D <- ggplot(panel_d, aes(factor(group), y, fill = which)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(0.7), width = 0.2) +
  labs(x = "major trajectory", y = "median distance (edges)", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

p_A
p_B
p_C
p_D


