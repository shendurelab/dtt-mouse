
################################################################################
### Figure 2. Single-nucleus RNA-seq of a DNA Typewriter embryo enables 
### reconstruction of a time-calibrated, zygote-rooted, tip-annotated phylogeny.


###############################################################################
### Fig. 2A: 2D UMAP co-embedding of 1.58M scRNA-seq profiles from embryo #3 and 
### 1.65M from 7 timepoints of our high-temporal-resolution single-cell atlas of mouse development

### scRNAseq_processing/step4_integration_atlas_E1275_E1425.R



###############################################################################
### Fig. 2B: Left: log2 percentage of cells per cell type (169 cell types from panel A), 
### each point colored by major trajectory. Right: Spearman correlation of cell type proportions, embryo #3 vs. each wildtype atlas timepoint (E12.75–E14.25).

### Please see the script: scRNAseq_processing/step4_integration_atlas_E1275_E1425.R



###############################################################################
### Fig. 2C: Histogram of the number of edited monomers per cell in embryo #3, 
### for the 4,371 cells with tape genotypes at all 11 tape integrations.

library(ggplot2)
library(tidyr)
library(dplyr)

dat = read.csv("./figures_data/e3v8.saturation_sites_written.csv", skip = 2)

p = ggplot(dat, aes(x = sites_written, y = pct_of_cells)) +
    geom_bar(stat="identity") +
    labs(x = "# of sites written", y = "% of cells", fill = NULL) +
    theme_classic()


###################################################################################################
### Fig. 2D: Cell-level rarefaction analysis of lineage-genotype diversity (Hurlbert interpolation)

library(ggplot2)
library(tidyr)
library(dplyr)

dat = read.csv("./figures_data/e3v8.lineage_rarefaction_curve.csv")

p = ggplot(data = filter(dat, region == "observed"), aes(x = cells_sampled, y = expected_genotypes)) +
    geom_point() +
    labs(x = "# of genome equivalents sampled", y = "Expected distinct lineage genotypes", fill = NULL) +
    theme_classic()


###############################################################################################
### Fig. 2E: Log2-scaled percentages of cells per cell type in blastomere A versus blastomere B

library(ggplot2)
library(tidyr)
library(dplyr)

major_trajectory_color_plate = c(
    "Neuroectoderm_and_glia"             = "#f96100",
    "Intermediate_neuronal_progenitors"  = "#2e0ab7",
    "Eye_and_other"                      = "#00d450",
    "Ependymal_cells"                    = "#b75bff",
    "CNS_neurons"                        = "#e5c000",
    "Mesoderm"                           = "#bb46c5",
    "Definitive_erythroid"               = "#dc453e",
    "Epithelium"                         = "#af9fb6",
    "Endothelium"                        = "#00a34e",
    "Muscle_cells"                       = "#ffa1f5",
    "Hepatocytes"                        = "#185700",
    "White_blood_cells"                  = "#7ca0ff",
    "Neural_crest_PNS_glia"              = "#fff167",
    "Adipocytes"                         = "#7f3e39",
    "Primitive_erythroid"                = "#ffa9a1",
    "Neural_crest_PNS_neurons"           = "#b5ce92",
    "T_cells"                            = "#ff9d47",
    "Lung_and_airway"                    = "#02b0d1",
    "Intestine"                          = "#ff007a",
    "B_cells"                            = "#01b7a6",
    "Olfactory_sensory_neurons"          = "#e6230b",
    "Cardiomyocytes"                     = "#643e8c",
    "Oligodendrocytes"                   = "#916e00",
    "Mast_cells"                         = "#005361",
    "Megakaryocytes"                     = "#3f283d",
    "Testis_and_adrenal"                 = "#585d3b"
)

dat = read.table("./support_data/e3v8.routing_labels.tsv.gz", sep='\t', header=T)
celltype = read.table(gzfile("./support_data/cell_metadata.v8.txt.gz"), sep='\t', header=T)
major_trajectory_celltype = read.table("./support_data/major_trajectory_celltype_table.txt", sep='\t', header=T)

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



######################################################################
### Fig. 2I: Mean edits accumulated per cell versus developmental time

library(ggplot2)
library(tidyr)
library(dplyr)

dat = read.csv("./figures_data/panel_F_temporal_resolution.csv")

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



####################################################################################################
### Fig. 2J: Observed tape editing rate across developmental time for blastomere A and blastomere B.

library(ggplot2)
library(tidyr)
library(dplyr)

dat = read.csv("./figures_data/fig2_ancestral_panel_D_editing_rate.csv")

p <- ggplot(dat[!is.na(dat$rate),], aes(x = day_mid, y = rate, color = blastomere)) +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = c(seq(1.5, 13.5, 2))) +
    scale_y_continuous(breaks = seq(0, 16, 2), limits = c(0, 16)) +
    scale_color_manual(values = c("A" = "#4783B5", "B" = "#F78C1E")) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = "Time (days)", y = "Editing rate (edits/day)", fill = NULL)
p

###################################################################################
### Fig. 2K: Internal node support across developmental time for blastomere A and B

library(ggplot2)
library(tidyr)
library(dplyr)

dat = read.csv("figures_data/fig2_ancestral_panel_E_support_over_time.csv")

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
p


