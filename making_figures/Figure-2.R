
################################################################################
### Figure 2. Single-nucleus RNA-seq of a DNA Typewriter embryo enables 
### reconstruction of a time-calibrated, zygote-rooted, tip-annotated phylogeny.

###############################################################################
### Fig. 2C: Histogram of the number of edited monomers per cell in embryo #3, 
### for the 4,371 cells with tape genotypes at all 11 tape integrations.

library(ggplot2)
library(tidyr)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/tape_pipeline/tables"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/e3v8.saturation_sites_written.csv"), skip = 2)

p = ggplot(dat, aes(x = sites_written, y = pct_of_cells)) +
    geom_bar(stat="identity") +
    labs(x = "# of sites written", y = "% of cells", fill = NULL) +
    theme_classic()

ggsave(paste0(save_path, "/Fig2/Fig2_saturation_sites_written.pdf"), p, height =4, width = 6)

x = rep(dat$sites_written, times = dat$n_cells)
print(paste0(mean(x), " +/- ", sd(x)))


###################################################################################################
### Fig. 2D: Cell-level rarefaction analysis of lineage-genotype diversity (Hurlbert interpolation)

library(ggplot2)
library(tidyr)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/tape_pipeline/tables"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/e3v8.lineage_rarefaction_curve.csv"))

p = ggplot(data = filter(dat, region == "observed"), aes(x = cells_sampled, y = expected_genotypes)) +
#    geom_line(linewidth = 1) +
    geom_point() +
    labs(x = "# of genome equivalents sampled", y = "Expected distinct lineage genotypes", fill = NULL) +
    theme_classic()

ggsave(paste0(save_path, "/Fig2/Fig2_lineage_rarefaction_curve.pdf"), p, height = 4, width = 6)


###############################################################################################
### Fig. 2E: Log2-scaled percentages of cells per cell type in blastomere A versus blastomere B

library(ggplot2)
library(tidyr)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/tape_pipeline/tables"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"
work_path = "/Volumes/f0085ts/work/tapemouse/tree_analysis"

dat = read.table(paste0(data_path, "/v8/e3v8.routing_labels.tsv.gz"), sep='\t', header=T)
celltype = read.table(paste0(work_path, "/cell_metadata.v8.txt"), sep='\t', header=T)
major_trajectory_celltype = read.table(paste0(work_path, "/major_trajectory_celltype_table.txt"), sep='\t', header=T)

all_celltypes = unique(celltype$celltype[celltype$cell_id %in% dat$cell[dat$blastomere %in% c("B1","B2")]])
### 141 cell types

cell_num_1 = celltype %>% 
    filter(cell_id %in% dat$cell[dat$blastomere == "B1"]) %>% 
    group_by(celltype) %>% 
    tally() %>% 
    complete(celltype = all_celltypes, fill = list(n = 0)) %>%
    mutate(total_n = sum(n)) %>%
    mutate(log2_frac = log2(100*(n/total_n)+1)) %>%
    select(celltype, B1_log2_frac = log2_frac)

cell_num_2 = celltype %>% 
    filter(cell_id %in% dat$cell[dat$blastomere == "B2"]) %>% 
    group_by(celltype) %>% 
    tally() %>% 
    complete(celltype = all_celltypes, fill = list(n = 0)) %>%
    mutate(total_n = sum(n)) %>%
    mutate(log2_frac = log2(100*(n/total_n)+1)) %>%
    select(celltype, B2_log2_frac = log2_frac)

df = cell_num_1 %>% left_join(cell_num_2, by = "celltype") %>%
    left_join(major_trajectory_celltype, by = "celltype")

p = ggplot(df, aes(x = B1_log2_frac, y = B2_log2_frac, color = major_trajectory)) +
    geom_point(size = 3) +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    labs(x = "Log2[Fraction (%) + 1], B1 blastomere", y = "Log2[Fraction (%) + 1], B2 blastomere") +
    scale_color_manual(values=major_trajectory_color_plate)

ggsave(paste0(save_path, "/Fig2/Fig2_celltype_frac_two_blastomere.pdf"), p, height = 5, width = 5)

fit = cor.test(df$B1_log2_frac, df$B2_log2_frac, method = "spearman")
print(fit$estimate) ### 0.9931551
print(fit$p.value)  ### 1.617634e-131


######################################################################
### Fig. 2I: Mean edits accumulated per cell versus developmental time

library(ggplot2)
library(tidyr)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/figures/fig_3"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/panel_F_temporal_resolution.csv"))

dat$rate[dat$architecture == "empirical_11x6" & dat$day == 0] = 13.33

p = ggplot(dat, aes(x = day, y = rate, color = architecture)) +
    geom_line(size = 2) +
    scale_x_continuous(breaks = c(seq(1.5, 13.5, 2))) +
    scale_y_continuous(breaks = seq(0, 16, 2), limits = c(0, 16)) +
    labs(x = "Time (embryonic day)", y = "Edits per cell per day", fill = NULL) +
    theme_classic() +
    scale_color_manual(values = c("sequential_1x66_flat" = "#4076bb", "non_sequential_1x66" = "#d14b5c",
                                  "empirical_11x6" = "black", "sequential_11x6_constant_rate" = "#108441")) +
    theme(legend.position = "none")

ggsave(paste0(save_path, "/Fig2/Fig2I.pdf"), p, height = 4, width = 5.5)


####################################################################################################
### Fig. 2J: Observed tape editing rate across developmental time for blastomere A and blastomere B.

library(ggplot2)
library(tidyr)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/figures/fig_3"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/panel_D_editing_rate.csv"))

p <- ggplot(dat[!is.na(dat$rate),], aes(x = day_mid, y = rate, color = blastomere)) +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = c(seq(1.5, 13.5, 2))) +
    scale_y_continuous(breaks = seq(0, 16, 2), limits = c(0, 16)) +
    scale_color_manual(values = c("A" = "#4783B5", "B" = "#F78C1E")) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = "Time (days)", y = "Editing rate (edits/day)", fill = NULL)

ggsave(paste0(save_path, "/Fig2/Fig2J.pdf"), p, height = 4, width = 5)


###################################################################################
### Fig. 2K: Internal node support across developmental time for blastomere A and B

library(ggplot2)
library(tidyr)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/figures/fig_3"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/panel_E_support_over_time.csv"))

# Reshape to long format so both pct columns can be plotted as separate lines
dat_long <- dat %>%
    pivot_longer(cols = c(pct_ge1, pct_ge2),
                 names_to = "threshold", values_to = "pct") %>%
    mutate(threshold = recode(threshold,
                              pct_ge1 = ">1",
                              pct_ge2 = ">2"))

p <- ggplot(dat_long[!is.na(dat_long$pct),], aes(x = day_mid, y = pct,
                          color = blastomere, linetype = threshold,
                          group = interaction(blastomere, threshold))) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = c(seq(1.5, 13.5, 2))) +
    scale_color_manual(values = c("A" = "#4783B5", "B" = "#F78C1E")) +
    scale_linetype_manual(values = c(">1" = "solid", ">2" = "dotted")) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = "Time (days)", y = "% of internal nodes",
         color = "Blastomere", linetype = "Threshold")

ggsave(paste0(save_path, "/Fig2/Fig2K.pdf"), p, height = 4, width = 5)



