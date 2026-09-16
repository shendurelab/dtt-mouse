

############################################################################
### Figure S17. Clade co-occurrence analyses yield consistent patterns 
### when conducted independently on subtrees defined by blastomeres A and B


################################################################################################
### Fig. 17C: Co-occurrence enrichment at recent vs. ancient clade scales is highly reproducible

library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyverse)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/tape_pipeline/figures/v8/coupling_depth"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/c18_s8c_k15_vs_k200_AB.csv"))

df_A = dat[,c("log2_enr_A_K200", "log2_enr_A_K15", "category")]
colnames(df_A) = c("log2_enr_K200", "log2_enr_K15", "category")
df_A$blastomere = "blastomere_A"

df_B = dat[,c("log2_enr_B_K200", "log2_enr_B_K15", "category")]
colnames(df_B) = c("log2_enr_K200", "log2_enr_K15", "category")
df_B$blastomere = "blastomere_B"

df = rbind(df_A, df_B)

p = ggplot() +
    geom_point(data = df, aes(x = log2_enr_K200, y = log2_enr_K15), color = "grey70", size = 1) +
    geom_point(data = df[df$category != "neither",], aes(x = log2_enr_K200, y = log2_enr_K15), color = "white", size = 2) +
    geom_point(data = df[df$category != "neither",], aes(x = log2_enr_K200, y = log2_enr_K15, color = category), size = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    facet_wrap(~ blastomere, nrow = 1) +
    labs(x = "log2 enrichment, K = 200", y = "log2 enrichment, K = 15", fill = NULL) +
    scale_x_continuous(breaks = seq(-6, 5, by = 2), limits = c(-6, 5)) +
    scale_y_continuous(breaks = seq(-6, 5, by = 2), limits = c(-6, 5)) +
    theme_classic() +
    scale_color_manual(
        values = c(K15_only  = "#F28E2B",
                   K200_only = "#1F77B4",
                   both_K      = "#7A1F5C"))

ggsave(paste0(save_path, "/FigS17/FigS8C.pdf"), p, height = 5, width = 11)



#######################################################################
### Fig. 17D: For each of the 130 significantly coupled cell-type pairs

library(ggplot2)
library(tidyr)
library(dplyr)
library(viridis)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/tape_pipeline/figures/v8/coupling_depth"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

df = read.csv(paste0(data_path, "/figSXa_coupling_depth_AB.csv"))

df_count = df %>% group_by(coupling_depth_A_E, coupling_depth_B_E) %>% tally()
df_count$n = factor(df_count$n)

p = ggplot() +
    geom_point(data = df_count, aes(x = coupling_depth_A_E, y = coupling_depth_B_E), color = "grey70", size = 3) +
    geom_point(data = df_count, aes(x = coupling_depth_A_E, y = coupling_depth_B_E, color = n), size = 2.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    labs(x = "coupling depth, blastomere A", y = "coupling depth, blastomere B", fill = NULL) +
    scale_x_continuous(breaks = seq(7, 13, by = 1), limits = c(7, 13.5)) +
    scale_y_continuous(breaks = seq(7, 13, by = 1), limits = c(7, 13.5)) +
    theme_classic() +
    scale_color_viridis(discrete=TRUE) 


ggsave(paste0(save_path, "/FigS17/FigS8D.pdf"), p, height = 5, width = 6)

