
#########################
### main cell cluster ###
#########################

work_path = "/net/shendure/vol2/projects/cxqiu/work/tapemouse"
source("~/work/scripts/utils.R")

mouse_gene = read.table("/net/gs/vol1/home/cxqiu/work/tome/code/mouse.v37.geneID.txt", header=T, sep="\t", as.is=T)

experiment_id = "experiment1_20260618_seq4_AD"

fd = readRDS(paste0(work_path, "/df_gene.rds"))
rownames(fd) = fd$gene_ID
pd = readRDS(paste0(work_path, "/data_analysis/", experiment_id, "/df_cell.rds"))

pd$log2_umi = log2(pd$UMI_count)
pd$EXON_pct = 100 * pd$all_exon / (pd$all_exon + pd$all_intron)
print(nrow(pd))
### n = 297,290 cells

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
### n = 23029; 7.7%
pd = pd[!(pd$detected_doublets | pd$doublet_cluster | pd$doublets_DEG),]
pd$detected_doublets = pd$doublet_cluster = pd$doublets_DEG = NULL
print(nrow(pd))
### 274261 cells
saveRDS(pd, paste0(work_path, "/data_analysis/", experiment_id, "/pd.rds"))

### first we removed cells with exon% > 85%
### then we identified the peak of the histgram and then set the log2_UMI cutoff as 
### [mean(x) - sd(x), mean(x) + 2*sd(x)]

x_tmp = pd$log2_umi[pd$EXON_pct <= 85]
x1 = mean(x_tmp) - sd(x_tmp)
x2 = mean(x_tmp) + 2*sd(x_tmp)

pd = pd[pd$log2_umi >= x1 & pd$log2_umi <= x2 & pd$EXON_pct <=85,]
print(nrow(pd))
### n = 222532 cells

### In the big dataset, we did those filtering as well
keep = pd$doublet_score <= 0.1 &
    pd$RIBO_pct <= 5 &
    pd$MT_pct <= 5
pd = pd[keep,]

print(nrow(pd))
print(median(pd$UMI_count))
print(median(pd$gene_count))
### n = 218654 cells
### median(UMI_count) = 1169
### median(gene_count) = 795

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

### output mtx and csv to generate h5ad used for analyzing using Scanpy/Python
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

work_path = '/net/shendure/vol2/projects/cxqiu/work/tapemouse'

experiment_id = "experiment1_20260618_seq4_AD"

adata = sc.read_mtx(f"{work_path}/data_analysis/{experiment_id}/h5ad/gene_count.mtx")
adata.obs.index = pd.read_csv(f"{work_path}/data_analysis/{experiment_id}/h5ad/df_cell.csv", index_col = 0, header=None).index.rename(None)
adata.var.index = pd.read_csv(f"{work_path}/data_analysis/{experiment_id}/h5ad/df_gene.csv", index_col = 0, header=None).index.rename(None)

adata.write(f"{work_path}/data_analysis/{experiment_id}/adata.h5ad", compression="gzip")








###############################
### check cutoff of UMI #######
###############################

setwd("~/work/tapemouse/data_analysis/experiment1_20260618_seq4_AD")
library(dplyr)
library(ggplot2)
library(gridExtra) 

pd = readRDS("pd.rds")
print(dim(pd))

p1 = ggplot(pd, aes(log2_umi)) + geom_histogram(binwidth = 0.1) 
p2 = ggplot(pd, aes(RIBO_pct)) + geom_histogram(binwidth = 0.1)
p3 = ggplot(pd, aes(MT_pct)) + geom_histogram(binwidth = 0.1) 
p4 = ggplot(pd, aes(EXON_pct)) + geom_histogram(binwidth = 0.1) + geom_vline(xintercept = 85) 

pdf("quality_summary.pdf",6,5)
grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2) 
dev.off()

pd %>% 
    filter(EXON_pct <= 85) %>%
    ggplot(aes(log2_umi)) + geom_histogram(binwidth = 0.1) 


pd_sub = pd[pd$EXON_pct <= 85,]

x1 = mean(pd_sub$log2_umi) - sd(pd_sub$log2_umi)
x2 = mean(pd_sub$log2_umi) + 2*sd(pd_sub$log2_umi)
hist(pd_sub$log2_umi, 500); abline(v = x1); abline(v = x2)

sum(pd$log2_umi >= x1 & pd$log2_umi <= x2 & pd$EXON_pct <=85)
### 222532

pd_sub = pd[pd$log2_umi >= x1 & pd$log2_umi <= x2 & pd$EXON_pct <=85,]
median(pd_sub$UMI_count) ### 1184
median(pd_sub$gene_count) ### 803

p = pd %>% 
    filter(EXON_pct <= 85) %>%
    ggplot(aes(log2_umi)) + 
    geom_histogram(binwidth = 0.05) +
    geom_vline(xintercept = x1, colour = "#BB0000") +
    geom_vline(xintercept = x2, colour = "blue") +
    labs(x = "", y = "", title = "experiment1_20260618_seq4_AD") +
    theme_classic(base_size = 15) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
pdf("Histgram_log2umi.pdf", 5, 3)
p
dev.off()







