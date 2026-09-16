
##############################################################################
### Figure 1. Pronuclear zygotic injection (PNI) of DNA Typewriter components 
### identifies an embryo with robust levels of sequential editing by E13.5.

#############################################################################################################
### Fig. 1B: Number of distinct tape-BC sequences (unique integration barcodes; blue) 
### and the estimated total number of genomic integrations after accounting for multi-copy barcodes (orange).

library(ggplot2)
library(tidyr)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/tape_pipeline/tables"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/panel1_barcode_copies.csv"))
dat = dat[dat$embryo %in% paste0("embryo ", c(2,3,6)),]

df = dat %>% group_by(embryo) %>% 
    summarize(num_tapebc = n_distinct(tapebc), sum_copies = sum(copies))

df_long <- df %>% 
    pivot_longer(cols = c(num_tapebc, sum_copies), 
                 names_to = "metric", 
                 values_to = "value")

p = ggplot(df_long, aes(x = embryo, y = value, fill = metric)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = value), 
              position = position_dodge(width = 0.9), 
              vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = c("num_tapebc" = "#4C72B0", 
                                 "sum_copies" = "#DD8452"),
                      labels = c("num_tapebc" = "# TapeBc",
                                 "sum_copies" = "# Copy")) +
    labs(x = NULL, y = "Count", fill = NULL) +
    theme_classic() +
    theme(legend.position = "top")

ggsave(paste0(save_path, "/Fig1/Fig1_integration_barcode_summary.pdf"), p, height = 4, width = 4.5)


#############################################################################################################
### Fig. 1C: Per-monomer editing rates for each embryo, expressed as the fraction of resolved reads edited 
### at each of the six monomers (colored 1-6), collapsed with respect to independent integrations.

library(ggplot2)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/tape_pipeline/tables"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/panel2_persite_by_embryo.csv"))
dat = dat[dat$embryo %in% paste0("embryo ", c(2,3,6)),]
dat$site = factor(paste0("site_", dat$site), levels = paste0("site_", 1:6))

p = ggplot(dat, aes(x = embryo, y = edit_rate, fill = site)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = NULL, y = "Editing rate", fill = NULL) +
    theme_classic() +
    theme(legend.position = "top") +
    scale_fill_manual(values = site_color_plate) +
    guides(fill = guide_legend(nrow = 1))

ggsave(paste0(save_path, "/Fig1/Fig1_per_site_editing_rate_per_embryo.pdf"), p, height = 4, width = 4.5)


####################################################################
### Fig. 1D: Per-monomer editing rate in embryo #3, shown separately 
### for each of 11 single-copy tape-BC integrations.

library(ggplot2)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/tape_pipeline/tables"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/panel3_persite_by_integration.csv"))
dat = dat[dat$embryo == "embryo 3",]
dat$site = factor(dat$site, levels = 1:6)

p = ggplot(dat, aes(x = site, y = edit_rate, group = tapebc, color = tapebc)) +
    geom_line() +
    geom_point() +
    labs(x = "Site", y = "Editing rate") +
    scale_color_manual(values = tapebc_color_plate) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                       limits = c(0, 1)) +
    theme_classic()

ggsave(paste0(save_path, "/Fig1/Fig1_per_site_editing_rate_per_integration_embryo3.pdf"), p, height = 4, width = 6)



####################################################################
### Fig. 1E: Frequency of each insertional symbol at each of the six 
### monomers in embryo #3, aggregated across integrations.

library(ggplot2)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/tape_pipeline/tables"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/DTTz_3_S3.site_insertion_freq.csv")) 

top_ins = dat %>% filter(freq_within_site >= 0.005) %>% 
    pull(insertion)
top_ins = unique(top_ins)

dat = dat %>%
    mutate(insertion = if_else(insertion %in% top_ins, insertion, "other"),
           insertion = factor(insertion, levels = c(top_ins, "other")))

p = ggplot(dat, aes(x = factor(site), y = freq_within_site, fill = insertion)) +
    geom_col() +
    theme_classic() +
    scale_fill_manual(values = edits_color_plate) +
    labs(x = "Site (position in TAPE)", y = "Fraction of edits", fill = NULL) 

ggsave(paste0(save_path, "/Fig1/Fig1_site_insertion_freq.pdf"), p, height = 3, width = 4)

### The 5 symbols that are present at ≥0.5% frequency only at site-1 and/or site-2 are marked with an asterisk
site_1 = dat %>% filter(site %in% c(1,2), freq_within_site >= 0.005) %>% pull(insertion) %>% unique()
site_2 = dat %>% filter(site %in% c(3:6), freq_within_site >= 0.005) %>% pull(insertion) %>% unique()
x = site_1[!site_1 %in% site_2]

