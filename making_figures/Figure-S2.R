
####################################################################
### Figure S2. Patterns and diversity of tape genotypes in embryo #3

#####################################################################
### Fig. S2C: Rarefaction analysis of recorded diversity in embryo #3

library(ggplot2)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/tape_pipeline/tables"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/DTTz_3_S3.rarefaction_curve.csv"))

p = ggplot(data = filter(dat, region == "observed"), aes(x = cells_sampled, y = expected_lineages, color = tapebc)) +
    geom_point() +
    scale_color_manual(values = tapebc_color_plate) +
    labs(x = "# of genome equivalents sampled", y = "Expected distinct lineage genotypes", fill = NULL) +
    theme_classic()

ggsave(paste0(save_path, "/FigS2/FigS1_rarefaction_curve.pdf"), p, height = 3, width = 6)


