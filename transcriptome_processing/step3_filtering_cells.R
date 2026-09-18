
#########################
### main cell cluster ###
#########################


mouse_gene = read.table("mouse.v37.geneID.txt", header=T, sep="\t", as.is=T)

experiment_id = "experiment1_20260618_seq4_AD"

fd = readRDS(paste0(work_path, "/df_gene.rds"))
rownames(fd) = fd$gene_ID
pd = readRDS(paste0(work_path, "/data_analysis/", experiment_id, "/df_cell.rds"))

pd$log2_umi = log2(pd$UMI_count)
pd$EXON_pct = 100 * pd$all_exon / (pd$all_exon + pd$all_intron)
print(nrow(pd))

### calculate MT_pct and Ribo_pct per cell
MT_gene = as.vector(fd[grep("^mt-",fd$gene_short_name),]$gene_ID)
Rpl_gene = as.vector(fd[grep("^Rpl",fd$gene_short_name),]$gene_ID)
Mrpl_gene = as.vector(fd[grep("^Mrpl",fd$gene_short_name),]$gene_ID)
Rps_gene = as.vector(fd[grep("^Rps",fd$gene_short_name),]$gene_ID)
Mrps_gene = as.vector(fd[grep("^Mrps",fd$gene_short_name),]$gene_ID)
RIBO_gene = c(Rpl_gene, Mrpl_gene, Rps_gene, Mrps_gene)

pd_tmp = NULL
batch_num = 8
for(i in 1:batch_num){
    print(paste0("processing:",i,"/",batch_num))
    count_i = readRDS(paste0(work_path, "/data_analysis/", experiment_id, "/gene_count_", i, ".rds"))
    pd_i = pd[colnames(count_i),]
    pd_i$MT_pct = 100 * Matrix::colSums(count_i[fd$gene_ID %in% MT_gene, ])/Matrix::colSums(count_i)
    pd_i$RIBO_pct = 100 * Matrix::colSums(count_i[fd$gene_ID %in% RIBO_gene, ])/Matrix::colSums(count_i)
    pd_tmp = rbind(pd_tmp, pd_i)
}
pd_tmp = pd_tmp[rownames(pd),]
pd$MT_pct = as.vector(pd_tmp$MT_pct)
pd$RIBO_pct = as.vector(pd_tmp$RIBO_pct)

pd_doublets_DEG = readRDS(paste0(work_path, "/data_analysis/", experiment_id,  "/doublet_cluster_2/res_doubelt_DEG.rds"))
pd$doublets_DEG = pd$cell_id %in% pd_doublets_DEG$cell_id[pd_doublets_DEG$doublet_DEG]

print(sum(pd$detected_doublets | pd$doublet_cluster | pd$doublets_DEG))
pd = pd[!(pd$detected_doublets | pd$doublet_cluster | pd$doublets_DEG),]
pd$detected_doublets = pd$doublet_cluster = pd$doublets_DEG = NULL
print(nrow(pd))
saveRDS(pd, paste0(work_path, "/data_analysis/", experiment_id, "/pd.rds"))

x_tmp = pd$log2_umi[pd$EXON_pct <= 85]
x1 = mean(x_tmp) - sd(x_tmp)
x2 = mean(x_tmp) + 2*sd(x_tmp)

pd = pd[pd$log2_umi >= x1 & pd$log2_umi <= x2 & pd$EXON_pct <=85,]
print(nrow(pd))

keep = pd$doublet_score <= 0.1 &
    pd$RIBO_pct <= 5 &
    pd$MT_pct <= 5
pd = pd[keep,]

print(nrow(pd))
print(median(pd$UMI_count))
print(median(pd$gene_count))


count = NULL
for(i in 1:batch_num){
    print(paste0("processing:",i,"/",batch_num))
    count_i = readRDS(paste0(work_path, "/data_analysis/", experiment_id, "/gene_count_", i, ".rds"))
    count_i = count_i[, colnames(count_i) %in% rownames(pd)]
    count = cbind(count, count_i)
}
print(sum(!colnames(count) %in% rownames(pd)))
count = count[,rownames(pd)]

rownames(pd) = colnames(count) = pd$cell_id = 
    paste0("exp1_", colnames(count))

saveRDS(pd, paste0(work_path, "/data_analysis/", experiment_id, "/pd_filter.rds"))

writeMM(t(count), paste0(work_path, "/data_analysis/", experiment_id, "/h5ad/gene_count.mtx"))
write.table(rownames(count), paste0(work_path, "/data_analysis/", experiment_id, "/h5ad/df_gene.csv"), row.names=F, col.names=F, quote=F, sep=',')
write.table(colnames(count), paste0(work_path, "/data_analysis/", experiment_id, "/h5ad/df_cell.csv"), row.names=F, col.names=F, quote=F, sep=',')

rm doublet_cluster_2/obj_all.rds
rm doublet_cluster_2/cluster_*_pd.rds



#############################
### using Python to save h5ad

import scanpy as sc
import pandas as pd
import numpy as np
import os, sys


experiment_id = "experiment1_20260618_seq4_AD"

adata = sc.read_mtx(f"{work_path}/data_analysis/{experiment_id}/h5ad/gene_count.mtx")
adata.obs.index = pd.read_csv(f"{work_path}/data_analysis/{experiment_id}/h5ad/df_cell.csv", index_col = 0, header=None).index.rename(None)
adata.var.index = pd.read_csv(f"{work_path}/data_analysis/{experiment_id}/h5ad/df_gene.csv", index_col = 0, header=None).index.rename(None)

adata.write(f"{work_path}/data_analysis/{experiment_id}/adata.h5ad", compression="gzip")













