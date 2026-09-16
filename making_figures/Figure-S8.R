

#####################################################
### Figure S8. Reliability of backbone tree topology.

##################################################################
### Fig. S8B: Leave-one-tape-out validation of reconstructed tree.

library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyverse)

data_path = "/Volumes/f0085ts/work/tapemouse/making_figures/FigSR7"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/panel_heldout_tape_distance_correlation.csv"))

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

ggsave(paste0(save_path, "/FigS8/FigSR7_Leave-one-tape-out.pdf"), p, height = 6, width = 6)
