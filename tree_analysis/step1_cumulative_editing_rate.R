
############################################################################################
### This script tests whether post-prefix editing rates are heterogeneous across cell types,
### and whether any such heterogeneity is explained by expression of PEmax or SAMHD1. It 
### generates the three panels shown in Fig. S4.

################################################
### Chengxiang Qiu
### Aug-23, 2026

print("Loading packages for regular data analysis, e.g. dplyr")
suppressMessages(library(Matrix))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))

print("Loading packages for plotting, e.g. ggplot2")
suppressMessages(library(ggplot2))
suppressMessages(library(plotly))
suppressMessages(library(htmlwidgets))
suppressMessages(library(gridExtra))
suppressMessages(library(gplots))
suppressMessages(library(patchwork))
suppressMessages(library(viridis))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggrepel))


#################################################################
### Section-1: Are editing rates heterogeneous across cell types?

### data can be found at Github: https://github.com/seidels/mouse_sprint/tree/v6-e3v5v6-merged/tape_pipeline/tables/v8
work_path = "/net/shendure/vol2/projects/cxqiu/work/tapemouse"

### cell-meta table
cell_meta = read.table(paste0(work_path, "/tree_analysis/cell_metadata.v8.txt"), header=T, sep="\t")

### pre-fix edits
prefix_edits = read.table(paste0(work_path, "/tree_analysis/edit_rate/e3.blastomere_defining_prefixes.tsv"), header=T)

### subset cells with blastomere A and B
dat_A = read.table(paste0(work_path, "/tree_analysis/e3v8.B1_tape_consensus.tsv.gz"), header = TRUE)
dat_B = read.table(paste0(work_path, "/tree_analysis/e3v8.B2_tape_consensus.tsv.gz"), header = TRUE)

celltype_A = cell_meta %>% filter(cell_id %in% dat_A$cell_id) %>%
    group_by(celltype) %>% tally() %>% mutate(frac = 100*n/nrow(dat_A)) %>% filter(frac >= 0.01)

celltype_B = cell_meta %>% filter(cell_id %in% dat_B$cell_id) %>%
    group_by(celltype) %>% tally() %>% mutate(frac = 100*n/nrow(dat_B)) %>% filter(frac >= 0.01)

common_celltype = intersect(celltype_A$celltype, celltype_B$celltype)
### n = 102 abundant cell types (>= 0.01%) in both A and B

### analyze each blastomere
per_celltype = list()
for(blastomere in c("A", "B")){
  print(blastomere)
  if (blastomere == "A"){
    dat = read.table(paste0(work_path, "/tree_analysis/e3v8.B1_tape_consensus.tsv.gz"), header=T, sep="\t") ### n = 889,849 cells
  } else {
    dat = read.table(paste0(work_path, "/tree_analysis/e3v8.B2_tape_consensus.tsv.gz"), header=T, sep="\t") ### n = 653,335 cells
  }
    
  long = dat[,c(1:12)] %>%
  pivot_longer(-cell_id, names_to = "integration", values_to = "chain") %>%
  filter(!is.na(chain)) %>% left_join(prefix_edits[prefix_edits$blastomere == blastomere,] %>% select(integration, n_defining_sites), by = "integration")

  per_tape_per_cell = long %>%
    rowwise() %>%
    mutate(
      sites_edited = sum(strsplit(chain, "\\|")[[1]] != "U") - n_defining_sites
    ) %>%
    ungroup() %>%
    filter(n_defining_sites != 6, sites_edited >= 0)
  ### Removing genotypes that were already fully edited during the prefix, as well as a handful of genotypes whose edits don't match the prefix pattern.

  per_tape = per_tape_per_cell %>%
    left_join(cell_meta[,c("cell_id", "celltype")], by = "cell_id") %>%
    group_by(celltype, integration) %>% 
    summarize(mean_per_tape = mean(sites_edited)/13.5)

  per_celltype[[blastomere]] = per_tape %>%
    group_by(celltype) %>%
    summarize(mean_per_celltype = sum(mean_per_tape))
}

df = per_celltype[["A"]] %>% select(celltype, A_edits = mean_per_celltype) %>% 
  inner_join(per_celltype[["B"]] %>% select(celltype, B_edits = mean_per_celltype), by = "celltype") %>%
  filter(celltype %in% common_celltype)

major_trajectory_celltype = read.table(paste0(work_path, "/tree_analysis/major_trajectory_celltype_table.txt"), header=T, sep="\t")
df = df %>% left_join(major_trajectory_celltype, by = "celltype")

df$category = rep("Other", nrow(df))
df$category[df$major_trajectory %in% c("B_cells", "Definitive_erythroid", "Megakaryocytes", 
                                       "T_cells", "Mast_cells", "Primitive_erythroid", "White_blood_cells")] = "Blood"

p = ggplot() +
  geom_point(data = df %>% filter(category != "Blood"),
             aes(x = A_edits, y = B_edits),
             color = "grey70", size = 3) +
  geom_point(data = df %>% filter(category == "Blood"),
             aes(x = A_edits, y = B_edits),
             color = "red", size = 3) +
  geom_text_repel(data = df %>% filter(category == "Blood"),
                  aes(x = A_edits, y = B_edits, label = celltype),
                  color = "black", size = 3,
                  box.padding = 0.4, max.overlaps = Inf) +
  scale_x_continuous(breaks = seq(0, 3, by = 0.5)) +
  scale_y_continuous(breaks = seq(0, 3, by = 0.5)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "Post-prefix editing rate (blastomere A)",
       y = "Post-prefix editing rate (blastomere B)")

ggsave(paste0(work_path, "/tree_analysis/edit_rate/FigS4_post_prefix_edit_rate.pdf"), p, height = 5, width = 5)


fit = cor.test(df$A_edits, df$B_edits, method = "spearman")
print(fit$estimate) ### 0.92

n <- sum(complete.cases(df$A_edits, df$B_edits))
rho <- fit$estimate
t_stat <- rho * sqrt((n - 2) / (1 - rho^2))
p_exact <- 2 * pt(-abs(t_stat), df = n - 2)
p_exact
### < 1e-42

saveRDS(df, paste0(work_path, "/tree_analysis/edit_rate/post_prefix_edit_rate.rds"))

### n = 102 abundant cell types (>0.01%) overlapped between A and B

print(max(df$A_edits)/min(df$A_edits)) ### 1.607212
print(max(df$B_edits)/min(df$B_edits)) ### 2.060734



#################################################################
### Section-2: Does PEmax expression correlated to editing rates?

experiment_list = c("experiment1_20260618_seq4_AD", 
                    "experiment1_20260618_seq4_EH",
                    "experiment1_20260618_seq5_IL",
                    "experiment1_20260618_seq5_MP",
                    "experiment1_20260618_seq6_QT",
                    "experiment1_20260618_seq6_UW",
                    "experiment2_20260713_seq2_XY")

batch_num = 8

### cell-meta table
cell_meta = read.table(paste0(work_path, "/tree_analysis/cell_metadata.v8.txt"), header=T, sep="\t")

PEmax_count = NULL
for(experiment_id in experiment_list){
  print(experiment_id)
    for(batch_id in 1:batch_num){
    x = read.csv(paste0(work_path, "/data_processing/", experiment_id, "/nobackup/output_", batch_id, "/PEmax_count.txt"))
    if (experiment_id == "experiment2_20260713_seq2_XY"){
      x$cell_id = paste0("exp2_", x$cell_id)
    } else {
      x$cell_id = paste0("exp1_", x$cell_id)
    }
    x = x[x$cell_id %in% cell_meta$cell_id,]
    x$fpkm_dtomato = x$n_dtomato / 8.005 / (x$n_total / 1e6)
    PEmax_count = rbind(PEmax_count, x)
  }
}

saveRDS(PEmax_count, paste0(work_path, "/tree_analysis/edit_rate/PEmax_count.rds"))

df = readRDS(paste0(work_path, "/tree_analysis/edit_rate/post_prefix_edit_rate.rds"))

dat = read.table(paste0(work_path, "/tree_analysis/e3v8.B1_tape_consensus.tsv.gz"), header=T, sep="\t")
cell_meta_A = cell_meta %>% filter(cell_id %in% dat$cell_id, celltype %in% df$celltype) %>% left_join(PEmax_count[,c("cell_id", "fpkm_dtomato")], by = "cell_id")
PEmax_A = cell_meta_A %>% group_by(celltype) %>% summarize(A_PEmax = mean(log2(fpkm_dtomato + 1)))

dat = read.table(paste0(work_path, "/tree_analysis/e3v8.B2_tape_consensus.tsv.gz"), header=T, sep="\t")
cell_meta_B = cell_meta %>% filter(cell_id %in% dat$cell_id, celltype %in% df$celltype) %>% left_join(PEmax_count[,c("cell_id", "fpkm_dtomato")], by = "cell_id")
PEmax_B = cell_meta_B %>% group_by(celltype) %>% summarize(B_PEmax = mean(log2(fpkm_dtomato + 1)))

df = df %>% left_join(PEmax_A %>% select(celltype, A_PEmax), by = "celltype") %>% 
  left_join(PEmax_B %>% select(celltype, B_PEmax), by = "celltype")




############# BOXPLOT #########################

print(wilcox.test(df$A_PEmax[df$category == "Blood"], df$A_PEmax[df$category != "Blood"]))
### p-value = 5.72e-06

p1 = ggplot(data = df, aes(x = factor(category), y = A_PEmax, fill = category)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(width = 0.2) +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    scale_fill_manual(values=c("Blood" = "red", "other" = "grey60")) +
    labs(x = "", y = "Mean Log2(PEmax expression FPKM), blastomere A")


print(wilcox.test(df$B_PEmax[df$category == "Blood"], df$B_PEmax[df$category != "Blood"]))
### p-value = 1.312e-06

p2 = ggplot(data = df, aes(x = factor(category), y = B_PEmax, fill = category)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(width = 0.2) +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    scale_fill_manual(values=c("Blood" = "red", "other" = "grey60")) +
    labs(x = "", y = "Mean Log2(PEmax expression FPKM), blastomere B")

ggsave(paste0(work_path, "/tree_analysis/edit_rate/FigS4_PEmax_exp.pdf"), p1 + p2, height = 5, width = 10)




df_long <- df %>%
  pivot_longer(cols = c(A_edits, B_edits, A_PEmax, B_PEmax),
               names_to = c("replicate", ".value"),
               names_sep = "_")

fit = cor.test(df_long$edits[df_long$category == "Blood"], df_long$PEmax[df_long$category == "Blood"], method = "spearman")
print(fit$estimate) ### -0.05035577
print(fit$p.value)  ### 0.7988483

fit = cor.test(df_long$edits[df_long$category != "Blood"], df_long$PEmax[df_long$category != "Blood"], method = "spearman")
print(fit$estimate) ### -0.2521572
print(fit$p.value)  ### 0.000761667


p1 = ggplot() +
  geom_point(data = df_long %>% filter(category == "Blood"),
             aes(x = edits, y = PEmax),
             color = "red", size = 3) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "Edits",
       y = "PEmax")

p2 = ggplot() +
  geom_point(data = df_long %>% filter(category != "Blood"),
             aes(x = edits, y = PEmax),
             color = "grey70", size = 3) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "Edits",
       y = "PEmax")

ggsave(paste0(work_path, "/tree_analysis/edit_rate/FigS4_PEmax_exp_split_blood.pdf"), p1 + p2, height = 5, width = 10)


#######################################################################################
### Section-3: Do selected gene expressions (e.g., Samhd1) correlated to editing rates?

library(Seurat)

df_gene = readRDS(paste0(work_path, "/df_gene.rds"))
mouse_genes <- c("Samhd1", "Mlh1", "Pms2", "Msh2", "Msh3",
                 "Exo1", "Fen1", "Lig1", "Rrm1", "Rrm2")
gene_list = df_gene[df_gene$gene_short_name %in% mouse_genes,]

batch_num = 8

experiment_list = c("experiment1_20260618_seq4_AD", 
                    "experiment1_20260618_seq4_EH",
                    "experiment1_20260618_seq5_IL",
                    "experiment1_20260618_seq5_MP",
                    "experiment1_20260618_seq6_QT",
                    "experiment1_20260618_seq6_UW",
                    "experiment2_20260713_seq2_XY")

gene_count_all = NULL
for(experiment_id in experiment_list){
  print(experiment_id)
  obj = readRDS(paste0(work_path, "/data_analysis/", experiment_id, "/obj.rds"))
  obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  gene_count = GetAssayData(obj, slot = "data")[gene_list$gene_ID,]
  rownames(gene_count) = as.vector(gene_list$gene_short_name)
  gene_count_all = cbind(gene_count_all, gene_count)
}

saveRDS(gene_count_all, paste0(work_path, "/tree_analysis/edit_rate/selected_gene_exp.rds"))


### cell-meta table
cell_meta = read.table(paste0(work_path, "/tree_analysis/cell_metadata.v8.txt"), header=T, sep="\t")

df = readRDS(paste0(work_path, "/tree_analysis/edit_rate/post_prefix_edit_rate.rds"))

dat = read.table(paste0(work_path, "/tree_analysis/e3v8.B1_tape_consensus.tsv.gz"), header=T, sep="\t")
cell_meta_A = cell_meta %>% filter(cell_id %in% dat$cell_id, celltype %in% df$celltype)

dat = read.table(paste0(work_path, "/tree_analysis/e3v8.B2_tape_consensus.tsv.gz"), header=T, sep="\t")
cell_meta_B = cell_meta %>% filter(cell_id %in% dat$cell_id, celltype %in% df$celltype)

i = 2
  gene_name <- rownames(gene_count_all)[i]
  dat <- data.frame(cell_id = colnames(gene_count_all),
                    gene_exp = as.vector(gene_count_all[i,])) %>%
    left_join(cell_meta_A[,c("cell_id", "celltype")], by = "cell_id") %>%
    filter(!is.na(celltype)) %>%
    group_by(celltype) %>%
    summarize(mean_exp = mean(gene_exp))
  
  df_x <- df %>% select(category, celltype, edits = A_edits) %>% left_join(dat, by = "celltype")
  df_x_1 = df_x
  
  print(wilcox.test(df_x$mean_exp[df$category == "Blood"], df_x$mean_exp[df$category != "Blood"]))
  ### p-value = 0.000379

  p1 = ggplot(data = df_x, aes(x = factor(category), y = mean_exp, fill = category)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(width = 0.2) +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    scale_fill_manual(values=c("Blood" = "red", "other" = "grey60")) +
    labs(x = "", y = "Mean normalized Samhd1 expression, blastomere A")


  dat <- data.frame(cell_id = colnames(gene_count_all),
                    gene_exp = as.vector(gene_count_all[i,])) %>%
    left_join(cell_meta_B[,c("cell_id", "celltype")], by = "cell_id") %>%
    filter(!is.na(celltype)) %>%
    group_by(celltype) %>%
    summarize(mean_exp = mean(gene_exp))
  
  df_x <- df %>% select(category, celltype, edits = B_edits) %>% left_join(dat, by = "celltype")
  df_x_2 = df_x
  
  print(wilcox.test(df_x$mean_exp[df$category == "Blood"], df_x$mean_exp[df$category != "Blood"]))
  ### p-value = 0.0003519

  p2 = ggplot(data = df_x, aes(x = factor(category), y = mean_exp, fill = category)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(width = 0.2) +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    scale_fill_manual(values=c("Blood" = "red", "other" = "grey60")) +
    labs(x = "", y = "Mean normalized Samhd1 expression, blastomere B")

ggsave(paste0(work_path, "/tree_analysis/edit_rate/FigS4_Samhd1_exp.pdf"), p1 + p2, height = 5, width = 10)





df_long <- rbind(df_x_1, df_x_2)

fit = cor.test(df_long$edits[df_long$category == "Blood"], df_long$mean_exp[df_long$category == "Blood"], method = "spearman")
print(fit$estimate) ### -0.07881773
print(fit$p.value)  ### 0.6892061

fit = cor.test(df_long$edits[df_long$category != "Blood"], df_long$mean_exp[df_long$category != "Blood"], method = "spearman")
print(fit$estimate) ### 0.2100242
print(fit$p.value)  ### 0.005222277


p1 = ggplot() +
  geom_point(data = df_long %>% filter(category == "Blood"),
             aes(x = edits, y = mean_exp),
             color = "red", size = 3) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "Edits",
       y = "Samhd1")

p2 = ggplot() +
  geom_point(data = df_long %>% filter(category != "Blood"),
             aes(x = edits, y = mean_exp),
             color = "grey70", size = 3) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "Edits",
       y = "Samhd1")

ggsave(paste0(work_path, "/tree_analysis/edit_rate/FigS4_Samhd1_exp_split_blood.pdf"), p1 + p2, height = 5, width = 10)















