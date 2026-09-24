
###############################################################################
### Figure S10. Query cells are placed at short edit distances from their
### anchors and recover cell-type concordance nearly as well as backbone cells

library(ggplot2)
library(dplyr)

BLASTOMERE_COL <- c(B1 = "#2a78d6", B2 = "#e07b39")
BLASTOMERE_LAB <- c(B1 = "Blastomere A", B2 = "Blastomere B")

qc <- read.csv("./figures_data/fig_s10_placement_qc_combined.csv")
qc$lab <- BLASTOMERE_LAB[qc$side]

hist_panel <- function(dat, key, xlabel, integer = TRUE) {
    meds <- dat %>% group_by(side, lab) %>% summarise(med = median(.data[[key]]), .groups = "drop")
    p <- ggplot(dat, aes(x = .data[[key]], y = after_stat(count / sum(count)), fill = lab)) +
        geom_histogram(aes(y = after_stat(density * width)), position = "identity",
                        alpha = 0.45, binwidth = if (integer) 1 else NULL, bins = if (integer) NULL else 60) +
        geom_vline(data = meds, aes(xintercept = med, color = lab), linetype = "dashed", linewidth = 0.6) +
        scale_fill_manual(values = setNames(BLASTOMERE_COL[names(BLASTOMERE_LAB)], BLASTOMERE_LAB)) +
        scale_color_manual(values = setNames(BLASTOMERE_COL[names(BLASTOMERE_LAB)], BLASTOMERE_LAB)) +
        labs(x = xlabel, y = "Fraction of query cells", fill = NULL, color = NULL) +
        theme_classic(base_size = 12) +
        theme(legend.position = "top")
    p
}

##########################################
### Fig. S10B: shared sequential edits with anchor cell

p <- hist_panel(qc, "shared_edits", "Shared sequential edits with anchor cell", integer = TRUE)
p

##########################################
### Fig. S10C: edits contradicting sequential order

conf <- read.csv("./figures_data/fig_s10_panel_C_irrev_conflicts.csv")
conf$lab <- BLASTOMERE_LAB[conf$side]
conf$cat <- cut(conf$n_conflicts, breaks = c(-Inf, 0, 1, 2, 3, 4, Inf),
                labels = c("0", "1", "2", "3", "4", "≥5"))
conf_agg <- conf %>%
    group_by(side, lab, cat) %>%
    summarise(n = sum(n_cells), .groups = "drop") %>%
    group_by(side) %>%
    mutate(frac = n / sum(n)) %>%
    ungroup()

p <- ggplot(conf_agg, aes(x = cat, y = frac, fill = lab)) +
    geom_col(position = position_dodge(0.7), width = 0.6) +
    scale_fill_manual(values = setNames(BLASTOMERE_COL[names(BLASTOMERE_LAB)], BLASTOMERE_LAB)) +
    labs(x = "Edits contradicting sequential order", y = "Fraction of query cells", fill = NULL) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top")
p

##########################################
### Fig. S10D: private edits per query

p <- hist_panel(qc, "new_edits", "Private edits per query", integer = TRUE)
p

##########################################
### Fig. S10E: implied divergence time

p <- hist_panel(qc, "pendant_days", "Implied divergence time (days)", integer = FALSE)
p

##########################################
### Fig. S10G: Terminal topology by origin

library(ggplot2)
library(tidyr)
library(dplyr)

dat = read.csv("./figures_data/figS5EF_placement_concordance.csv")

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


p
############################################################################################################
### Fig. S10H: For each cell, the fraction of its equally related nearest neighbors that share its cell type

library(ggplot2)
library(tidyr)
library(dplyr)

dat = read.csv("./figures_data/figS5EF_placement_concordance.csv")

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



