

#############################################################################################
### Figure S4. Robustness of recording in DNA Typewriter embryos from additional PNI sessions


##################################################
### Fig. S4A: Number of distinct tape-BC sequences

library(ggplot2)
library(tidyr)
library(dplyr)

dat = read.csv("./figures_data/figSR5b_integrations.csv")

df_long <- dat %>% 
    pivot_longer(cols = c(n_tapebc, n_integrations), 
                 names_to = "metric", 
                 values_to = "value")

df_long$embryo = factor(df_long$embryo, levels = c("embryo #2", "embryo #3", "E14.5", "E18.5 #1", "E18.5 #2", "E19.5"))

p = ggplot(df_long, aes(x = embryo, y = value, fill = metric)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = value), 
              position = position_dodge(width = 0.9), 
              vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = c("n_tapebc" = "#4C72B0", 
                                 "n_integrations" = "#DD8452"),
                      labels = c("num_tapebc" = "# TapeBc",
                                 "n_integrations" = "# Copy")) +
    labs(x = NULL, y = "Count", fill = NULL) +
    theme_classic() +
    theme(legend.position = "top")


#######################################################
### Fig. S4B: Per-monomer editing rates for each embryo

library(ggplot2)
library(dplyr)

site_color_plate = c(
    "site_1" = "#440154",
    "site_2" = "#414487",
    "site_3" = "#2a788e",
    "site_4" = "#22a884",
    "site_5" = "#7ad151",
    "site_6" = "#ec1b37"
)

dat = read.csv("./figures_data/figSR5c_persite_editrate.csv")
dat$site = factor(paste0("site_", dat$site), levels = paste0("site_", 1:6))
dat$embryo = factor(dat$embryo, levels = c("embryo #2", "embryo #3", "E14.5", "E18.5 #1", "E18.5 #2", "E19.5"))

p = ggplot(dat, aes(x = embryo, y = edit_rate, fill = site)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = NULL, y = "Editing rate", fill = NULL) +
    theme_classic() +
    theme(legend.position = "top") +
    scale_fill_manual(values = site_color_plate) +
    guides(fill = guide_legend(nrow = 1))



###############################################################################
### Fig. S4C: Rarefaction analysis of recorded diversity in each of six embryos

library(ggplot2)
library(dplyr)

dat = read.csv("./figures_data/figSR5d_rarefaction_curve.csv")
dat$embryo = factor(dat$embryo, levels = c("embryo #2", "embryo #3", "E14.5", "E18.5 #1", "E18.5 #2", "E19.5"))

p = ggplot(data = filter(dat, region == "observed"), 
           aes(x = cells_sampled, y = expected_lineages, color = tapebc)) +
    geom_point(size = 0.3) +
    facet_wrap(~ embryo, nrow = 2, ncol = 3) +
    labs(x = "# of genome equivalents sampled", 
         y = "Expected distinct lineage genotypes", 
         fill = NULL) +
    theme_classic() +
    theme(legend.position = "none")


