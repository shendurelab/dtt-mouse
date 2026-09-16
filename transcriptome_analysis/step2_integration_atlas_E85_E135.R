
#######################################################################
### UMAP co-embedding of 1.58M scRNA-seq profiles from embryo #3 and
### ~5M from 21 timepoints of our high-temporal-resolution single-cell
### atlas of mouse development (E8.5 - E13.5 in 6-hr increments)


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
day_list = ["E8.5", "E8.75", "E9.0", "E9.25", "E9.5", "E9.75", 
"E10.0","E10.25","E10.5","E10.75",
"E11.0","E11.25","E11.5","E11.75",
"E12.0","E12.25","E12.5","E12.75",
"E13.0","E13.25","E13.5"]

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

adata.write(f"{work_path}/transcriptome_analysis/adata_integration_early.h5ad", compression="gzip")

adata.obs.to_csv(f"{work_path}/transcriptome_analysis/adata_integration_early.obs.csv")
pd.DataFrame(adata.obsm['X_pca']).to_csv(f"{work_path}/transcriptome_analysis/adata_integration_early.pca.csv")


########################
### Plotting the 3D UMAP

source("~/work/scripts/utils.R")
work_path = "/net/shendure/vol8/projects/cxqiu/work/tapemouse"
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/tapemouse"

pd = read.csv(paste0(work_path, "/transcriptome_analysis/adata_integration_early.obs.csv"), row.names=1)
pd$cell_id = rownames(pd)
pd_1 = pd[pd$dataset == 'jax',]
pd_2 = pd[pd$dataset == 'tapemouse',]

pd_jax = readRDS("/net/shendure/vol2/projects/cxqiu/JAX_rna_mm39/pd.rds")
rownames(pd_jax) = pd_jax$cell_id = paste0(pd_jax$experiment_id, "_", pd_jax$cell_id)
pd_1_x = pd_1 %>% left_join(pd_jax, by = "cell_id") %>% as.data.frame()
rownames(pd_1_x) = pd_1_x$cell_id

pd_new = readRDS(paste0(work_path, "/transcriptome_analysis/adata_integration.obs.rds"))
pd_2_x = pd_2 %>% left_join(pd_new, by = "cell_id") %>% as.data.frame()
rownames(pd_2_x) = pd_2_x$cell_id

pd_1$major_trajectory = pd_1_x$major_trajectory
pd_1$celltype = pd_1_x$celltype
pd_1$RT_group = pd_1_x$RT_group
pd_1$day = pd_1_x$day

pd_2$major_trajectory = pd_2_x$major_trajectory
pd_2$celltype = pd_2_x$celltype
pd_2$RT_group = pd_2_x$RT_group
pd_2$day = "E13.5"

pd = rbind(pd_1, pd_2)
saveRDS(pd, paste0(work_path, "/transcriptome_analysis/adata_integration_early.obs.rds"))

set.seed(2016)
pd_1 = pd %>% filter(dataset == "jax") %>% group_by(day) %>% slice_sample(n = 10000)
pd_2 = pd %>% filter(dataset == "tapemouse"); pd_2$day = "TapeMouse:E13.5"
pd_sub = rbind(pd_1, pd_2)

day_color_map = c("#5E4FA2", "#476AAE", "#3287BB", "#4BA4B0", "#66C1A4", "#87CFA4", "#AADCA3", "#C8E89D", "#E4F498", "#F2F9AA",
"#FDFDBC", "#FEEFA4", "#FDDF8A", "#FDC675", "#FCAD61", "#F78C51", "#F36C43", "#E35548", "#D43D4E", "#B81E47", "#9E0142", "#808080")

names(day_color_map) = c("E8.5", "E8.75", "E9.0", "E9.25", "E9.5", "E9.75", "E10.0", "E10.25", "E10.5", "E10.75",
"E11.0", "E11.25", "E11.5", "E11.75", "E12.0", "E12.25", "E12.5", "E12.75", "E13.0", "E13.25", "E13.5", "TapeMouse:E13.5")

pd_sub$day = factor(pd_sub$day, levels = names(day_color_map))
fig = plot_ly(pd_sub, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~day, colors = day_color_map)
saveWidget(fig, paste0(save_path, "/integration_E8_to_E13.5_day.html"), selfcontained = FALSE, libdir = "tmp")

fig = plot_ly(pd_sub, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~major_trajectory, colors = major_trajectory_color_plate)
saveWidget(fig, paste0(save_path, "/integration_E8_to_E13.5_major_trajectory.html"), selfcontained = FALSE, libdir = "tmp")

fig = plot_ly(pd_sub, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~celltype)
saveWidget(fig, paste0(save_path, "/integration_E8_to_E13.5_celltype.html"), selfcontained = FALSE, libdir = "tmp")





