

#####################################################
### Figure S8. Reliability of backbone tree topology.

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


