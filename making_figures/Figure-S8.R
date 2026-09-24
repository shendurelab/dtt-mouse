

#####################################################
### Figure S8. Reliability of backbone tree topology.

#####################################################################################################
### Fig. S8A: node support over developmental time, counted by distinct tapes edited (solid) vs.
### total edited sites (dashed, Fig. 2K's metric), for blastomere A and B

library(ggplot2)
library(tidyr)
library(dplyr)

dat = read.csv("./figures_data/fig_s8_ancestral_panel_A_support_over_time_indep_tapes.csv")

dat_long <- dat %>%
    pivot_longer(cols = c(pct_ge1, pct_ge2, pct_ge3), names_to = "threshold", values_to = "pct") %>%
    mutate(threshold = recode(threshold, pct_ge1 = "support ≥ 1", pct_ge2 = "support ≥ 2", pct_ge3 = "support ≥ 3"),
           threshold = factor(threshold, levels = c("support ≥ 1", "support ≥ 2", "support ≥ 3")),
           metric = factor(metric, levels = c("total_edits", "indep_tapes")))

p = ggplot(dat_long[!is.na(dat_long$pct),],
           aes(x = day_mid, y = pct, color = blastomere, linetype = metric,
               group = interaction(blastomere, metric))) +
    geom_line() +
    # points only on the "indep_tapes" (solid) line, so the dashed reference stays uncluttered
    geom_point(data = dat_long[!is.na(dat_long$pct) & dat_long$metric == "indep_tapes",]) +
    facet_wrap(~threshold, nrow = 1) +
    scale_x_continuous(breaks = seq(0, 13.5, 2), limits = c(0, 13.4)) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100)) +
    scale_color_manual(values = c("A" = "#4783B5", "B" = "#F78C1E")) +
    scale_linetype_manual(values = c("total_edits" = "dashed", "indep_tapes" = "solid"),
                          labels = c("total_edits" = "Total edits (previous)",
                                     "indep_tapes" = "Independent tapes")) +
    theme_classic() +
    theme(legend.position = "top", strip.background = element_rect(fill = "white")) +
    labs(x = "Developmental time (days)", y = "% of internal nodes with support ≥ threshold",
         color = "Blastomere", linetype = NULL)

p
##################################################################
### Fig. S8B: Leave-one-tape-out validation of reconstructed tree.

library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyverse)

tapebc_color_plate = c(
    "AATAGAAAACGA" = "#A6CEE3",
    "ATATCAAATTGA" = "#1F78B4",
    "CAGCTAACGCCT" = "#B2DF8A",
    "CATATAATCGCA" = "#33A02C",
    "CGGCGAAAAGGT" = "#FB9A99",
    "CGGGGAATTGTA" = "#E31A1C",
    "GTGTAAATCGGC" = "#FDBF6F",
    "TAACGAATGCCG" = "#FF7F00",
    "TCCGGAAGACCC" = "#CAB2D6",
    "TGACTAAAGCGG" = "#6A3D9A",
    "TGGGGAACATAT" = "#B15928"
)


dat = read.csv("./figures_data/panel_heldout_tape_distance_correlation.csv")

dat_long <- dat |>
    pivot_longer(
        cols = c(rho_leave_one_out, rho_random_mean),
        names_to = "type",
        values_to = "rho"
    ) |>
    mutate(
        type = factor(type, 
                      levels = c("rho_leave_one_out", "rho_random_mean"),
                      labels = c("Leave-one-out", "Random (mean)")),
        scenario = factor(scenario, levels = c("random", "clade"))
    )

p = ggplot(dat_long, aes(x = type, y = rho)) +
    geom_boxplot(outlier.shape = NA, fill = "grey90") +
    geom_jitter(aes(color = heldout_tape), width = 0.2, size = 1.2, alpha = 0.8) +
    facet_grid(blastomere ~ scenario) +
    labs(x = NULL, y = expression(rho)) +
    theme_classic() +
    theme(legend.position = "none") +
    scale_color_manual(values = tapebc_color_plate)


