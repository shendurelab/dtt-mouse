
############################################################################################
### Figure S15. Tree siblings share annotations far in excess of chance, and reproducibly so


#######################################################################################################
### Fig. 15D: Fold-enrichment of each progenitor & post-mitotic pairing computed separately for A and B


library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/tape_pipeline/figures/v8/heterotypic"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/figSX_captured_divisions_AB_replication.csv"))

library(ggplot2)
library(scales)

d <- subset(dat, tested_in_both == 1)
d <- d[order(d$captured), ]  # captured on top

sp <- cor(d$fold_A, d$fold_B, method = "spearman")
pe <- cor(log2(d$fold_A), log2(d$fold_B), method = "pearson")

d$group <- ifelse(d$captured == 1,
                  sprintf("captured division (%d)", sum(d$captured == 1)),
                  sprintf("tested in both halves (%d)", sum(d$captured == 0)))

p = ggplot(d, aes(fold_A, fold_B, color = group, size = group)) +
    geom_abline(slope = 1, linetype = "dashed", color = "grey55") +
    geom_hline(yintercept = 3, linetype = "dotted", color = "grey55") +
    geom_vline(xintercept = 3, linetype = "dotted", color = "grey55") +
    geom_point(alpha = 0.8) +
    scale_x_log10(limits = c(0.6, NA),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(limits = c(0.6, NA),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_color_manual(values = c("#C0503C", "#87A9CB")) +
    scale_size_manual(values = c(3, 1.6), guide = "none") +
    coord_fixed() +
    labs(x = "fold enrichment, blastomere A (B1)",
         y = "fold enrichment, blastomere B (B2)",
         color = NULL,
         title = sprintf("Spearman = %.2f    Pearson(log2) = %.2f", sp, pe)) +
    theme_classic()

ggsave(paste0(save_path, "/FigS15/FigS7C.pdf"), p, height = 5, width = 6)



