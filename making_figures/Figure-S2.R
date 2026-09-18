
####################################################################
### Figure S2. Patterns and diversity of tape genotypes in embryo #3

#####################################################################
### Fig. S2C: Rarefaction analysis of recorded diversity in embryo #3

library(ggplot2)
library(dplyr)

tapebc_color_plate = c(
    "AATAGAAAACGA" = "#A6CEE3",
    "ATATCAAATTGA" = "#1F78B4",
    "CAGCTAACGCCT" = "#B2DF8A",
    "CATATAATCGCA" = "#33A02C",
    "CGGCGAAAAGGT" = "#FB9A99",
    "CGGGGAATTGTA" = "#E31A1C",
    "GTGTAAATCGGC" = "#FDBF6F",
    "TAACGAATGCCG" = "#FF7F00",
    "TCCGGAAGACCC" = "#CAB2D6",
    "TGACTAAAGCGG" = "#6A3D9A",
    "TGGGGAACATAT" = "#B15928"
)

dat = read.csv("./figures_data/DTTz_3_S3.rarefaction_curve.csv")

p = ggplot(data = filter(dat, region == "observed"), aes(x = cells_sampled, y = expected_lineages, color = tapebc)) +
    geom_point() +
    scale_color_manual(values = tapebc_color_plate) +
    labs(x = "# of genome equivalents sampled", y = "Expected distinct lineage genotypes", fill = NULL) +
    theme_classic()

