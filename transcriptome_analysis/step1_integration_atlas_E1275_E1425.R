

#######################################################################
### UMAP co-embedding of 1.58M scRNA-seq profiles from embryo #3 and
### 1.65M from 7 timepoints of our high-temporal-resolution single-cell
### atlas of mouse development (E12.75 - E14.25 in 6-hr increments)

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import os, sys
import gc

work_path = '/net/shendure/vol2/projects/cxqiu/work/tapemouse'

experiment_list = ["experiment1_20260618_seq4_AD", 
                   "experiment1_20260618_seq4_EH",
                   "experiment1_20260618_seq5_IL",
                   "experiment1_20260618_seq5_MP",
                   "experiment1_20260618_seq6_QT",
                   "experiment1_20260618_seq6_UW",
                   "experiment2_20260713_seq2_XY"]

adatas = []
for experiment_id in experiment_list:
    print(experiment_id)
    a = sc.read_h5ad(f"{work_path}/data_analysis/{experiment_id}/adata.h5ad")
    adatas.append(a)

adata = ad.concat(adatas, axis=0)

del adatas
gc.collect()

mouse_gene = pd.read_csv("/net/gs/vol1/home/cxqiu/work/tome/code/mouse.v37.geneID.txt", sep="\t", index_col=4)
adata.var = mouse_gene.loc[adata.var_names]

# exclude sex + mito chromosomes, only keep lncRNA and protein_coding
exclude_chrom = ['chrX', 'chrY', 'chrM']
keep_type = ['lncRNA', 'protein_coding']
adata = adata[:, ~adata.var['chr'].isin(exclude_chrom) & adata.var['gene_type'].isin(keep_type)].copy()


# Load jax data and align genes
day_list = ["E12.75","E13.0","E13.25","E13.5","E13.75","E14.0","E14.25"]

adatas = []
for day_id in day_list:
    print(day_id)
    a = sc.read_h5ad(f"/net/shendure/vol2/projects/cxqiu/JAX_rna_mm39/gene_count/adata.{day_id}.h5ad")
    adatas.append(a)

adata_jax = ad.concat(adatas, axis=0)

del adatas
gc.collect()


# Merge two datasets by common genes
common = adata.var_names.intersection(adata_jax.var_names)
adata     = adata[:, common].copy()
adata_jax = adata_jax[:, common].copy()
print(f"Common genes: {len(common)}")

# Strip var to empty (matching jax) then concat
adata.var = adata.var[[]]
adata = ad.concat([adata, adata_jax], label="dataset", keys=["tapemouse", "jax"])
print(adata.shape, adata.obs['dataset'].value_counts().to_dict())
adata = adata.copy()

del adata_jax
gc.collect()

sc.pp.normalize_total(adata, target_sum=1e4)
print("Done normalization by total counts ...")

sc.pp.log1p(adata)
print("Done log transformation ...")

sc.pp.highly_variable_genes(adata, n_top_genes=2500)
print("Done finding highly variable genes ...")

adata = adata[:, adata.var.highly_variable]
print("Done filtering in highly variable genes ...")

sc.pp.scale(adata, max_value=10)
print("Done scaling data ...")

sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
print("Done performing PCA ...")

sc.pp.neighbors(adata, n_neighbors=50, n_pcs=50)
print("Done computing neighborhood graph ...")

sc.tl.umap(adata, min_dist=0.1, n_components=3)
adata.obs['UMAP_1'] = list(adata.obsm['X_umap'][:,0])
adata.obs['UMAP_2'] = list(adata.obsm['X_umap'][:,1])
adata.obs['UMAP_3'] = list(adata.obsm['X_umap'][:,2])

sc.tl.umap(adata, min_dist=0.1, n_components=2)
adata.obs['UMAP_2d_1'] = list(adata.obsm['X_umap'][:,0])
adata.obs['UMAP_2d_2'] = list(adata.obsm['X_umap'][:,1])
print("Done UMAP ...")

adata.write(f"{work_path}/transcriptome_analysis/adata_integration.h5ad", compression="gzip")

adata.obs.to_csv(f"{work_path}/transcriptome_analysis/adata_integration.obs.csv")
pd.DataFrame(adata.obsm['X_pca']).to_csv(f"{work_path}/transcriptome_analysis/adata_integration.pca.csv")


#########################################
### Performing knn to transferring labels

import pandas as pd
import numpy as np
from annoy import AnnoyIndex

work_path = '/net/shendure/vol2/projects/cxqiu/work/tapemouse'

OBS_PATH = f"{work_path}/transcriptome_analysis/adata_integration.obs.csv"
PCA_PATH = f"{work_path}/transcriptome_analysis/adata_integration.pca.csv"
N_NEIGHBORS = 20
N_TREES = 150
METRIC = "euclidean"

obs = pd.read_csv(OBS_PATH, index_col=0)
pca = pd.read_csv(PCA_PATH, index_col=0)

pca_values = pca.values
n_dims = pca_values.shape[1]

is_jax = (obs["dataset"] == "jax").values
is_tapemouse = (obs["dataset"] == "tapemouse").values

jax_ids = obs.index[is_jax].to_numpy()
tapemouse_ids = obs.index[is_tapemouse].to_numpy()

jax_pca = pca_values[is_jax]
tapemouse_pca = pca_values[is_tapemouse]

index = AnnoyIndex(n_dims, METRIC)
for i, vec in enumerate(jax_pca):
    index.add_item(i, vec)

index.build(N_TREES)

neighbor_idx = np.zeros((len(tapemouse_ids), N_NEIGHBORS), dtype=np.int64)
for ti, vec in enumerate(tapemouse_pca):
    neighbor_idx[ti] = index.get_nns_by_vector(vec, N_NEIGHBORS)

df = pd.DataFrame(neighbor_idx, columns=[f"neighbor_{i+1}" for i in range(N_NEIGHBORS)])
df.insert(0, "tapemouse_id", tapemouse_ids)
df.to_csv(f"{work_path}/transcriptome_analysis/adata_integration.knn.csv", index=False)
 
pd.Series(jax_ids).to_csv(f"{work_path}/transcriptome_analysis/adata_integration.knn_index.csv", index=True, header=["jax_id"])



########################
### Plotting the 3D UMAP

source("~/work/scripts/utils.R")
work_path = "/net/shendure/vol2/projects/cxqiu/work/tapemouse"
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/tapemouse"

pd = read.csv(paste0(work_path, "/transcriptome_analysis/adata_integration.obs.csv"), row.names=1)
pd$cell_id = rownames(pd)
pd_1 = pd[pd$dataset == 'jax',]
pd_2 = pd[pd$dataset == 'tapemouse',]

pd_jax = readRDS("/net/shendure/vol2/projects/cxqiu/JAX_rna_mm39/pd.rds")
pd_1_x = pd_1 %>% left_join(pd_jax, by = "cell_id") %>% as.data.frame()
rownames(pd_1_x) = pd_1_x$cell_id

knn = read.csv(paste0(work_path, "/transcriptome_analysis/adata_integration.knn.csv"), row.names=1)
knn = knn + 1
knn_id = read.csv(paste0(work_path, "/transcriptome_analysis/adata_integration.knn_index.csv"), row.names=1)
pd_1_x = pd_1_x[as.vector(knn_id$jax_id),]

tmp = NULL
for(i in 1:ncol(knn)){
    tmp = cbind(tmp, pd_1_x$major_trajectory[knn[,i]])
}
major_trajectory = apply(tmp, 1, function(x) names(which.max(table(x))))

tmp = NULL
for(i in 1:ncol(knn)){
    tmp = cbind(tmp, pd_1_x$celltype[knn[,i]])
}
celltype = apply(tmp, 1, function(x) names(which.max(table(x))))

pd_2_x = data.frame(cell_id = rownames(knn), major_trajectory = major_trajectory, celltype = celltype)
pd_2_x = pd_2 %>% left_join(pd_2_x, by = "cell_id")
pd_2_x$SampleName = "Embryo_3"

pd_1$major_trajectory = pd_1_x$major_trajectory
pd_1$celltype = pd_1_x$celltype
pd_1$SampleName = pd_1_x$SampleName
pd_1$day = pd_1_x$day

pd_2$major_trajectory = pd_2_x$major_trajectory
pd_2$celltype = pd_2_x$celltype
pd_2$SampleName = pd_2_x$SampleName
pd_2$day = "E13.5"

pd = rbind(pd_1, pd_2)
saveRDS(pd, paste0(work_path, "/transcriptome_analysis/adata_integration.obs.rds"))

set.seed(2016)
pd_sub = pd %>% group_by(dataset) %>% slice_sample(n = 150000)

fig = plot_ly(pd_sub, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~dataset)
saveWidget(fig, paste0(save_path, "/integration_dataset.html"), selfcontained = FALSE, libdir = "tmp")

fig = plot_ly(pd_sub[pd_sub$dataset == "tapemouse",], x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~celltype)
saveWidget(fig, paste0(save_path, "/tapemouse_celltype.html"), selfcontained = FALSE, libdir = "tmp")

fig = plot_ly(pd_sub[pd_sub$dataset == "tapemouse",], x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~major_trajectory, colors = major_trajectory_color_plate)
saveWidget(fig, paste0(save_path, "/tapemouse_major_trajectory.html"), selfcontained = FALSE, libdir = "tmp")

pd_out = pd[pd$dataset == "tapemouse", c("cell_id", "major_trajectory", "celltype", "UMAP_1", "UMAP_2", "UMAP_3", "UMAP_2d_1", "UMAP_2d_2")]

write.table(pd_out, paste0(save_path, "/cell_metadata.v8.txt"), row.names=F, quote=F, sep='\t')


https://shendure-web.gs.washington.edu/content/members/cxqiu/private/nobackup/tapemouse/cell_metadata.v8.txt



########################
### Plotting the 2D UMAP

source("~/work/scripts/utils.R")
work_path = "/net/shendure/vol2/projects/cxqiu/work/tapemouse"
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/tapemouse"

pd = readRDS(paste0(work_path, "/transcriptome_analysis/adata_integration.obs.rds"))

p = ggplot() +
    geom_point(data = pd, aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.15, color = "black") +
    geom_point(data = pd, aes(x = UMAP_2d_1, y = UMAP_2d_2, color = major_trajectory), size=0.1) +
    theme_void() +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=major_trajectory_color_plate)
ggsave("~/share/Fig2_umap_major_trajectory.png", p, dpi = 300, height = 5, width = 5)

p = ggplot() +
    geom_point(data = pd %>% filter(dataset == "jax"), aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.15, color = "grey80") +
    geom_point(data = pd %>% filter(dataset == "tapemouse"), aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.15, color = "red") +
    theme_void() +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5))
ggsave("~/share/Fig2_umap_tapemouse.png", p, dpi = 300, height = 5, width = 5)

p = ggplot() +
    geom_point(data = pd %>% filter(dataset == "tapemouse"), aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.15, color = "grey80") +
    geom_point(data = pd %>% filter(dataset == "jax"), aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.15, color = "blue") +
    theme_void() +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5))
ggsave("~/share/Fig2_umap_jax.png", p, dpi = 300, height = 5, width = 5)




#################################################################
### Compare cell type compositions between new data and JAX E13.5

pd = readRDS(paste0(work_path, "/transcriptome_analysis/adata_integration.obs.rds"))

major_trajectory_celltype = pd %>% filter(dataset == "jax") %>% group_by(major_trajectory, celltype) %>% tally() %>% select(-n)

library(tidyr)

all_celltypes <- unique(pd$celltype)

cell_num_1 = pd %>% 
    filter(dataset == "tapemouse") %>% 
    group_by(celltype) %>% 
    tally() %>% 
    complete(celltype = all_celltypes, fill = list(n = 0)) %>%
    mutate(total_n = sum(n)) %>%
    mutate(log2_frac = log2(100*(n/total_n)+1)) %>%
    select(celltype, new_log2_frac = log2_frac)

df = NULL
cor_res = NULL
for(i in names(table(pd$day[pd$dataset == "jax"]))){
    x = pd %>% 
        filter(dataset == "jax", day == i) %>% 
        group_by(celltype) %>% 
        tally() %>% 
        complete(celltype = all_celltypes, fill = list(n = 0)) %>%
        mutate(day = i, total_n = sum(n)) %>%
        mutate(log2_frac = log2(100*(n/total_n)+1)) %>%
        select(celltype, day, old_log2_frac = log2_frac) %>%
        left_join(cell_num_1, by = "celltype")
    fit = cor.test(x$new_log2_frac, x$old_log2_frac, method = "spearman")
    cor_res = rbind(cor_res, data.frame(day = i, corr = fit$estimate, pval = -log10(fit$p.val)))
    df = rbind(df, x)
}

df = df %>% left_join(major_trajectory_celltype, by = "celltype")

df_x = df %>% filter(day == "E13.5")
fit = cor.test(df_x$new_log2_frac, df_x$old_log2_frac, method = "spearman")
print(fit$estimate) ### 0.94
print(fit$p.val) ### < 1e-78

p = ggplot(df %>% filter(day == "E13.5"), aes(x = new_log2_frac, y = old_log2_frac, color = major_trajectory)) +
  geom_point(size = 3) +
  theme_classic(base_size = 12) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
  labs(x = "Log2[Fraction (%) + 1] in this study", y = "Log2[Fraction (%) + 1] in JAX E13.5") +
  scale_color_manual(values=major_trajectory_color_plate)

ggsave("~/share/Fig2_celltype_frac.pdf", p, height = 5, width = 5)


library(ggrepel)

p <- ggplot(cor_res, aes(x = day, y = corr)) +
  geom_line(aes(group = 1), color = "grey50") +
  geom_point(aes(color = day), size = 3) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color = "black")) +
  labs(y = "Spearman correlation coefficient") +
  scale_color_viridis(discrete=TRUE)

ggsave("~/share/Fig2_celltype_frac_2.pdf", p, height = 5, width = 4)




###########################################################################################
### Comparing cell-type-compositions between backbone tree vs. placed cells vs. E13.5 atlas

source("~/work/scripts/utils.R")
work_path = "/net/shendure/vol2/projects/cxqiu/work/tapemouse"
library(ape)
library(tidyr)

pd = readRDS(paste0(work_path, "/transcriptome_analysis/adata_integration.obs.rds"))
all_celltypes = unique(pd$celltype)

major_trajectory_celltype_table = read.table(paste0(work_path, "/tree_analysis/major_trajectory_celltype_table.txt"), sep="\t", header=T)


### subset cells with blastomere A and B
dat_A = read.table(paste0(work_path, "/tree_analysis/e3v8.B1_tape_consensus.tsv.gz"), header = TRUE)
dat_B = read.table(paste0(work_path, "/tree_analysis/e3v8.B2_tape_consensus.tsv.gz"), header = TRUE)

### E13.5 atlas
pd_E135 = pd %>% filter(day == "E13.5", dataset == "jax")

### backbone
tree_backbone = read.tree(paste0(work_path, "/tree_analysis/merged_minB2h_lineage_constrained.nwk"))

pd_backbone_A = pd %>% filter(dataset == "tapemouse", 
    cell_id %in% tree_backbone$tip.label, cell_id %in% dat_A$cell_id)

pd_backbone_B = pd %>% filter(dataset == "tapemouse", 
    cell_id %in% tree_backbone$tip.label, cell_id %in% dat_B$cell_id)

### full-tree
tree_full = read.tree(paste0(work_path, "/tree_analysis/merged_full_placed.nwk"))

pd_placed_A = pd %>% filter(dataset == "tapemouse", 
    cell_id %in% tree_full$tip.label, cell_id %in% dat_A$cell_id, !cell_id %in% tree_backbone$tip.label)

pd_placed_B = pd %>% filter(dataset == "tapemouse", 
    cell_id %in% tree_full$tip.label, cell_id %in% dat_B$cell_id, !cell_id %in% tree_backbone$tip.label)


### comparision and making plot

get_frac <- function(pd, all_celltypes, min_frac = 0.01) {
  pd %>%
    group_by(celltype) %>%
    tally() %>%
    complete(celltype = all_celltypes, fill = list(n = 0)) %>%
    mutate(total_n = sum(n),
           frac = 100 * (n / total_n),
           log2_frac = log2(frac + 1)) %>%
    filter(frac >= min_frac) %>%
    select(celltype, log2_frac)
}

plot_celltype_cor <- function(pd_x, pd_y,
                              x_label, y_label,
                              all_celltypes,
                              celltype_table = major_trajectory_celltype_table,
                              color_plate = major_trajectory_color_plate,
                              min_frac = 0.01) {

  df <- get_frac(pd_x, all_celltypes, min_frac) %>%
    rename(x_log2_frac = log2_frac) %>%
    inner_join(get_frac(pd_y, all_celltypes, min_frac) %>%
                 rename(y_log2_frac = log2_frac),
               by = "celltype") %>%
    left_join(celltype_table, by = "celltype")

  fit     <- cor.test(df$x_log2_frac, df$y_log2_frac, method = "spearman")
  rho_lab <- round(unname(fit$estimate), 2)

  # bound for tiny p, exact value otherwise
  ttl <- if (fit$p.value < 1e-3) {
    bquote("Spearman's" ~ rho == .(rho_lab) * "," ~
             italic(p) < 10^.(ceiling(log10(fit$p.value))))
  } else {
    bquote("Spearman's" ~ rho == .(rho_lab) * "," ~
             italic(p) == .(signif(fit$p.value, 2)))
  }

  ggplot(df, aes(x = x_log2_frac, y = y_log2_frac, color = major_trajectory)) +
    geom_point(size = 3) +
    theme_classic(base_size = 10) +
    theme(legend.position = "none",
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black"),
          plot.title  = element_text(hjust = 0.5, size = 12)) +
    labs(x = paste0("Log2(% of cells per cell type), ", x_label),
         y = paste0("Log2(% of cells per cell type), ", y_label),
         title = ttl) +
    scale_color_manual(values = color_plate)
}

p1 <- plot_celltype_cor(pd_backbone_A, pd_placed_A,
                        "backbone-blastomere-A", "placed-blastomere-A",
                        all_celltypes)

p2 <- plot_celltype_cor(pd_backbone_B, pd_placed_B,
                        "backbone-blastomere-B", "placed-blastomere-B",
                        all_celltypes)

p3 <- plot_celltype_cor(pd_backbone_A, pd_E135,
                        "backbone-blastomere-A", "atlas @E13.5",
                        all_celltypes)

p4 <- plot_celltype_cor(pd_backbone_B, pd_E135,
                        "backbone-blastomere-B", "atlas @E13.5",
                        all_celltypes)

p5 <- plot_celltype_cor(pd_placed_A, pd_E135,
                        "placed-blastomere-A", "atlas @E13.5",
                        all_celltypes)

p6 <- plot_celltype_cor(pd_placed_B, pd_E135,
                        "placed-blastomere-B", "atlas @E13.5",
                        all_celltypes)

pp <- (p1 | p3 | p5) / (p2 | p4 | p6) &
  theme(plot.margin = margin(t = 20, r = 20, b = 20, l = 20))

ggsave("~/share/celltype_cor_grid.pdf", pp, width = 15, height = 10)






