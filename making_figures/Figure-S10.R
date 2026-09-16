
###############################################################################
### Figure S10. Query cells are placed at short edit distances from their 
### anchors and recover cell-type concordance nearly as well as backbone cells

##########################################
### Fig. S10G: Terminal topology by origin

library(ggplot2)
library(tidyr)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/tape_pipeline/figures/v8"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/figS5EF_placement_concordance.csv"))

plot_dat <- dat %>%
    filter(tier != "overall", side != "all") %>%
    mutate(
        bar  = paste(side, origin, sep = "_"),
        tier = factor(tier, levels = c("loner", "polytomy", "cherry"))
    )

plot_dat$bar <- factor(plot_dat$bar, levels = c("B1_backbone_orig", "B2_backbone_orig",
                                                "B1_backbone", "B2_backbone",
                                                "B1_placed", "B2_placed"))

p <- ggplot(plot_dat, aes(x = bar, y = pct_of_origin, fill = tier)) +
    geom_col(width = 0.75) +
    geom_text(aes(label = sprintf("%.0f", pct_of_origin)),
              position = position_stack(vjust = 0.5),
              color = "black", size = 3.5) +
    scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 101),
                       expand = c(0, 0)) +
    labs(x = NULL, y = "Percent of origin", fill = "Tier") +
    theme_classic(base_size = 12) +
    scale_fill_manual(values = c("loner" = "#B84A45", "polytomy" = "#DDA13A", "cherry" = "#55A868")) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1),
          panel.grid.major.x = element_blank())

ggsave(paste0(save_path, "/FigS10/FigS5E.pdf"), p, height = 4, width = 7)


############################################################################################################
### Fig. S10H: For each cell, the fraction of its equally related nearest neighbors that share its cell type

library(ggplot2)
library(tidyr)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/tape_pipeline/figures/v8"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/figS5EF_placement_concordance.csv"))

plot_dat <- dat %>%
    filter(origin != "backbone_orig", side != "all") %>%
    mutate(
        group = paste(origin, side),
        group = factor(group, levels = c("backbone B1", "backbone B2",
                                         "placed B1",   "placed B2")),
        tier  = factor(tier,  levels = c("cherry", "polytomy", "loner", "overall"))
    )

p = ggplot(plot_dat, aes(x = tier, y = concordance, fill = group)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    geom_text(aes(label = sprintf("%.1f", fold)),
              position = position_dodge(width = 0.8),
              vjust = -0.3, size = 3) +
    geom_hline(yintercept = 5.7, linetype = "dashed", color = "grey40") +
    scale_fill_manual(values = c(
        "backbone B1" = "#4783B5",   # B1 blue
        "backbone B2" = "#F78C1E",   # B2 orange
        "placed B1"   = "#B8D3E6",   # B1 blue, lighter
        "placed B2"   = "#FCD9B0"    # B2 orange, lighter
    )) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = "% nearest neighbor shares cell type", fill = NULL) +
    theme_classic(base_size = 12) +
    theme(panel.grid.major.x = element_blank())

ggsave(paste0(save_path, "/FigS10/FigS5F.pdf"), p, height = 4, width = 8)



