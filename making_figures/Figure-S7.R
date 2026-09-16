
###########################################################################
### Figure S7. A cell-number ceiling corrects implausibly early node dates, 
### and most internal branches carry multiple independent supporting edits


##########################################################################################################
### Fig. S7A: Number of inferred lineages or anticipated number of cells in embryo over developmental time

library(ggplot2)
library(tidyr)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/figures/fig_S3"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/panel_B_ltt_vs_ceiling.csv"))

### regroup
summary(dat$day[dat$series == "cell ceiling, unused (> E7.5)"])
summary(dat$day[dat$series == "cell ceiling, used (≤ E7.5)"])

series_update = dat$series
series_update[dat$day > 6.5 & dat$series %in% c("cell ceiling, unused (> E7.5)", "cell ceiling, used (≤ E7.5)")] = "cell ceiling, unused (> E6.5)"
series_update[dat$day <= 6.5 & dat$series %in% c("cell ceiling, unused (> E7.5)", "cell ceiling, used (≤ E7.5)")] = "cell ceiling, used (≤ E6.5)"
dat$series = series_update

summary(dat$day[dat$series == "cell ceiling, unused (> E6.5)"])
summary(dat$day[dat$series == "cell ceiling, used (≤ E6.5)"])

dat$log10_count = log10(dat$count)

p = ggplot(dat[dat$log10_count > 1,], aes(x = day, y = log10_count, color = series)) +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = seq(1.5, 13.5, 2)) +
    scale_y_continuous(breaks = 1:7) +
    labs(x = "Time (embryonic day)", y = "Number of inferred lineages or cells in embryo", fill = NULL) +
    theme_classic() +
    scale_color_manual(values = c("unconstrained dating" = "#4076bb", "cell ceiling, used (≤ E6.5)" = "#d14b5c",
                                  "cell ceiling, unused (> E6.5)" = "#8a8881", "constrained dating" = "#108441")) +
    theme(legend.position = "none")

ggsave(paste0(save_path, "/FigS7/FigS3A_update.pdf"), p, height = 5, width = 6.5)


##########################################################################################################
### Fig. S7B: Scatter plot of shifts in dating (constrained minus unconstrained) vs. the unconstrained date

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrastr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/figures/fig_S3"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/panel_C_date_shift.csv"))

p = ggplot(dat, aes(x = old, y = shift, color = blast)) +
    rasterise(geom_point(size = 0.1), dpi = 300) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_x_continuous(breaks = seq(1.5, 13.5, 2)) +
    scale_y_continuous(breaks = seq(0, 4, 1)) +
    labs(x = "Unconstrained node date (embryonic day)", 
         y = "Date shift, new - old (days)", fill = NULL) +
    theme_classic() +
    scale_color_manual(values = c("B1" = "#4783B5", "B2" = "#F78C1E")) +
    theme(legend.position = "none")

ggsave(paste0(save_path, "/FigS7/FigS3B.pdf"), p, height = 5, width = 6.5)


###################################################################################
### Fig. S7C: Shift in node dates, shown separately for blastomere A and B lineages

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrastr)

data_path = "/Users/cxqiu/GitHub/dtt-mouse-analysis/figures/sensitivity"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/node_time_shifts.csv"))

p <- ggplot(dat, aes(x = original_depth, y = depth_diff_hours, color = sensitivity)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "grey70") +
    geom_hline(yintercept = c(-3, 3), linetype = "dashed", color = "black") +
    rasterise(
        geom_point(alpha = 0.15, size = 0.6),
        dpi = 300
    ) +
    facet_wrap(~ blastomere, nrow = 1) +
    scale_color_manual(values = c("diff_58_42" = "#8E4A9E",
                                  "diff_63_37" = "#2CA4A4")) +
    scale_x_continuous(breaks = seq(1.5, 13.5, by = 1)) +
    scale_y_continuous(breaks = c(-6, -3, -1, 0, 1, 3, 6)) +
    labs(x = "Original node time [days]",
         y = "Node time difference [hours]",
         color = "Sensitivity analysis") +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    theme_classic() +
    theme(legend.position = "top",
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "white"))

ggsave(paste0(save_path, "/FigS7/FigSR1_node_time_shifts.pdf"), p, width = 10, height = 5)



#####################################################################################################################
### Fig. S7D: For each internal node, support is quantified as the number of edits accumulated on the incoming branch

library(ggplot2)
library(tidyr)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/figures/fig_S3"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/panel_A_support_histogram.csv"))

dat_frac <- dat %>%
    group_by(blastomere, support) %>%
    summarise(n = n(), .groups = "drop_last") %>%
    mutate(frac = n / sum(n)) %>%
    ungroup()

medians <- dat %>%
    group_by(blastomere) %>%
    summarise(med = median(support), .groups = "drop")

p = ggplot(dat_frac, aes(x = support, y = frac, fill = blastomere)) +
    geom_col(position = position_dodge(preserve = "single")) +
    geom_vline(data = medians,
               aes(xintercept = med, color = blastomere),
               linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
    scale_x_continuous(breaks = seq(0, max(dat_frac$support), by = 5)) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Support", y = "Fraction of nodes", fill = NULL) +
    theme_classic() +
    scale_fill_manual(values = c("B1" = "#4783B5", "B2" = "#F78C1E")) +
    theme(legend.position = "none")

ggsave(paste0(save_path, "/FigS7/FigS3C.pdf"), p, height = 5, width = 6.5)
