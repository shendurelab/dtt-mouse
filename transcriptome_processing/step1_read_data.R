
#############################################################
### Processing the TapeMouse experiment1_20260618_seq4_AD
### Please contact: CX Qiu (Chengxiang.Qiu@dartmouth.edu)

##########################################################################################
### Step-1: read the raw data from processing pipeline, roughly filtering UMI < 100 nuclei

library(Matrix)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(viridis)
library(patchwork)

work_path = "/net/shendure/vol2/projects/cxqiu/work/tapemouse"

mouse_gene = read.table("/net/gs/vol1/home/cxqiu/work/tome/code/mouse.v37.geneID.txt", header=T, sep="\t", as.is=T)

experiment_id = "experiment1_20260618_seq4_AD"

batch_num = 8

### how many reads in this experiment (after UMI attach)
read_num_fastq = NULL
for(cnt in 1:batch_num){
    print(cnt)
    read_num_cnt = read.table(paste0(work_path, "/data_processing/", experiment_id, "/nobackup/output_", cnt, "/read_num_UMI_attach.txt"),as.is=T)
    read_num_fastq = rbind(read_num_fastq, read_num_cnt)
}
print(sum(read_num_fastq$V1))
### 1,752,100,203

### summary the duplication rate
read_num = NULL
for(cnt in 1:batch_num){
  print(cnt)
  read_num_cnt = read.table(paste0(work_path, "/data_processing/", experiment_id, "/nobackup/output_", cnt, "/read_num.txt"),as.is=T)
  read_num = rbind(read_num, read_num_cnt)
}

print(summary(1 - read_num$V2/read_num$V1)) 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.1017  0.6384  0.6615  0.6587  0.6891  0.7480

print(sum(read_num$V1))
### 1,618,837,684

print(sum(read_num$V2))
### 533,319,833

### saving and performing doublets removing on individual batches if necessary

df_cell_merge = NULL

for(cnt in 1:batch_num){
  print(cnt)
  load(paste0(work_path, "/data_processing/", experiment_id, "/nobackup/output_", cnt, "/report/sci_summary.RData"))
  colnames(df_gene) = c("gene_ID", "gene_type", "gene_short_name")

  keep = df_gene$gene_ID != "PEmax_dTomato"
  df_gene = df_gene[keep,]
  gene_count = gene_count[keep,]

  print(sum(rownames(gene_count) != df_gene$gene_ID))
  print(sum(colnames(gene_count) != df_cell$sample))
  
  rownames(df_gene) = rownames(gene_count) = df_gene$gene_ID =
      unlist(lapply(rownames(df_gene), function(x) strsplit(x,"[.]")[[1]][1]))

  df_cell$UMI_count = Matrix::colSums(gene_count)
  gene_count_copy = gene_count
  gene_count_copy@x[gene_count_copy@x > 0] = 1
  df_cell$gene_count = Matrix::colSums(gene_count_copy)

  RT_barcode = read.csv(paste0(work_path, "/data_processing/", experiment_id, "/RTsamplesheet.csv"), as.is=T, header=T)
  RT_barcode = RT_barcode[RT_barcode$RTprimer == "oligodT",]
  RT_barcode$RTprimer = NULL
  
  df_cell$RTwell = sub(".*_RT-", "", df_cell$sample)
  df_cell = df_cell %>% rename(cell_id = sample) %>% left_join(RT_barcode, by = "RTwell")
  print(sum(is.na(df_cell$SampleName)))

  keep = df_cell$SampleName == "DTT_Z_858923_E2"
  df_cell = df_cell[keep,]
  gene_count = gene_count[,keep]
  
  df_cell$experiment = "experiment1_20260618"
  df_cell$SampleName = "Embryo_3"
  
  df_gene_x = df_gene %>% left_join(mouse_gene, by = "gene_ID")
  df_gene$chr = as.vector(df_gene_x$chr)
  gene_keep = df_gene$chr %in% paste0("chr", c(1:19, "M", "X", "Y"))
  
  keep = df_cell$UMI_count >= 100 & 
    df_cell$gene_count >= 50 & 
    df_cell$unmatched_rate < 0.4
  df_cell = df_cell[keep,]
  df_gene = df_gene[gene_keep,]
  gene_count = gene_count[gene_keep, keep]
  rownames(gene_count) = as.vector(df_gene$gene_ID)
  colnames(gene_count) = as.vector(df_cell$cell_id)
  rownames(df_cell) = as.vector(df_cell$cell_id)
  print(dim(gene_count))
  
  saveRDS(gene_count, paste0(work_path, "/data_analysis/", experiment_id, "/gene_count_", cnt, ".rds"))
  df_cell_merge = rbind(df_cell_merge, df_cell)  
}

print(nrow(df_cell_merge))
### nrow(df_cell) = 297,290

print(median(df_cell_merge$UMI_count))
### 1,052

print(median(df_cell_merge$gene_count))
### 732

print(sum(is.na(df_cell_merge$SampleName)))
### 0

saveRDS(df_cell_merge, paste0(work_path, "/data_analysis/", experiment_id, "/df_cell.rds"))


################################################
### output dataset used for running scrublet ###
################################################

scrublet_bacth = 4
df_cell_merge$split_batch = sample(c(1:scrublet_bacth), nrow(df_cell_merge), replace = T)
print(table(df_cell_merge$split_batch))

gene_count_1 = NULL; df_cell_1 = NULL
gene_count_2 = NULL; df_cell_2 = NULL
gene_count_3 = NULL; df_cell_3 = NULL
gene_count_4 = NULL; df_cell_4 = NULL

for(cnt in 1:batch_num){
  gene_count = readRDS(paste0(work_path, "/data_analysis/", experiment_id, "/gene_count_", cnt, ".rds"))
  df_cell = df_cell_merge[colnames(gene_count),]

  index_1 = df_cell$split_batch == 1
  gene_count_1 = cbind(gene_count_1, gene_count[,index_1])
  df_cell_1 = rbind(df_cell_1, df_cell[index_1,])
  
  index_2 = df_cell$split_batch == 2
  gene_count_2 = cbind(gene_count_2, gene_count[,index_2])
  df_cell_2 = rbind(df_cell_2, df_cell[index_2,])
  
  index_3 = df_cell$split_batch == 3
  gene_count_3 = cbind(gene_count_3, gene_count[,index_3])
  df_cell_3 = rbind(df_cell_3, df_cell[index_3,])
  
  index_4 = df_cell$split_batch == 4
  gene_count_4 = cbind(gene_count_4, gene_count[,index_4])
  df_cell_4 = rbind(df_cell_4, df_cell[index_4,])
}

print(sum(colnames(gene_count_1) != rownames(df_cell_1)))
print(sum(colnames(gene_count_2) != rownames(df_cell_2)))
print(sum(colnames(gene_count_3) != rownames(df_cell_3)))
print(sum(colnames(gene_count_4) != rownames(df_cell_4)))

writeMM(t(gene_count_1), paste0(work_path, "/data_analysis/", experiment_id, "/gene_count_1.mtx"))
writeMM(t(gene_count_2), paste0(work_path, "/data_analysis/", experiment_id, "/gene_count_2.mtx"))
writeMM(t(gene_count_3), paste0(work_path, "/data_analysis/", experiment_id, "/gene_count_3.mtx"))
writeMM(t(gene_count_4), paste0(work_path, "/data_analysis/", experiment_id, "/gene_count_4.mtx"))

write.csv(df_cell_1, paste0(work_path, "/data_analysis/", experiment_id, "/df_cell_1.csv"))
write.csv(df_cell_2, paste0(work_path, "/data_analysis/", experiment_id, "/df_cell_2.csv"))
write.csv(df_cell_3, paste0(work_path, "/data_analysis/", experiment_id, "/df_cell_3.csv"))
write.csv(df_cell_4, paste0(work_path, "/data_analysis/", experiment_id, "/df_cell_4.csv"))

write.csv(df_gene, paste0(work_path, "/data_analysis/", experiment_id, "/df_gene.csv"))

### run scrublet using python to detect doublets

###########################################
### after running scrublet using python ###
###########################################

df_cell = readRDS(paste0(work_path, "/data_analysis/", experiment_id, "/df_cell.rds"))
scrublet_bacth = 4

doublet_scores_observed_cells = NULL
doublet_scores_simulated_doublets = NULL
df = NULL
for(i in 1:scrublet_bacth){
  print(i)
  df_i = read.csv(paste0(work_path, "/data_analysis/", experiment_id, "/df_cell_", i, ".csv"), header=T, row.names=1, as.is=T)
  
  doublet_scores_observed_cells_i = read.csv(paste0(work_path, "/data_analysis/", experiment_id, "/doublet_scores_observed_cells_", i, ".csv"), header=F)
  df_i$doublet_score = as.vector(doublet_scores_observed_cells_i$V1)
  df = rbind(df, df_i)
  
  doublet_scores_simulated_doublets_i = read.csv(paste0(work_path, "/data_analysis/", experiment_id, "/doublet_scores_simulated_doublets_", i, ".csv"), header=F)
  doublet_scores_simulated_doublets = rbind(doublet_scores_simulated_doublets, doublet_scores_simulated_doublets_i)
}

df = df[rownames(df_cell),]
print(sum(rownames(df) != rownames(df_cell)))

df_cell$doublet_score = as.vector(df$doublet_score)
df_cell$detected_doublets = df_cell$doublet_score > 0.2

### sum(df_cell$detected_doublets)/nrow(df_cell) = 0.03942951

###############################################################
### checking if sub-clusters include over 15% doublet cells ###
###############################################################

global = read.csv(paste0(work_path, "/data_analysis/", experiment_id, "/doublet_cluster/global.csv"), header=T)
main_cluster_list = sort(as.vector(unique(global$louvain)))

res = NULL

for(i in 1:length(main_cluster_list)){
  print(paste0(i, "/", length(main_cluster_list)))
  dat = read.csv(paste0(work_path, "/data_analysis/", experiment_id, "/doublet_cluster/adata.obs.louvain_", (i-1), ".csv"), header=T)
  print(nrow(dat))
  dat$louvain = as.vector(paste0("cluster_", dat$louvain))
  dat = dat %>%
    left_join(df_cell[,c("cell_id", "detected_doublets", "doublet_score")], by = "cell_id")
  
  tmp2 = dat %>%
    group_by(louvain) %>%
    tally() %>%
    dplyr::rename(n_sum = n)
  
  tmp1 = dat %>%
    filter(detected_doublets == "TRUE") %>%
    group_by(louvain) %>%
    tally() %>%
    left_join(tmp2, by = "louvain") %>%
    mutate(frac = n/n_sum) %>%
    filter(frac > 0.15)
  
  dat$doublet_cluster = dat$louvain %in% as.vector(tmp1$louvain) 
  
  p1 = ggplot(dat, aes(umap_1, umap_2, color = louvain)) + geom_point() + theme(legend.position="none") 
  p2 = ggplot(dat, aes(umap_1, umap_2, color = doublet_cluster)) + geom_point()
  p3 = ggplot(dat, aes(umap_1, umap_2, color = detected_doublets)) + geom_point()
  p4 = ggplot(dat, aes(umap_1, umap_2, color = doublet_score)) + geom_point() + scale_color_viridis(option = "plasma")
  p5 = ggplot(dat, aes(umap_1, umap_2, color = UMI_count)) + geom_point() + scale_color_viridis(option = "plasma")
  p6 = ggplot(dat, aes(umap_1, umap_2, color = gene_count)) + geom_point() + scale_color_viridis(option = "plasma")

  combined_plot = (p1 | p2) / (p3 | p4) / (p5 | p6)
  ggsave(paste0(work_path, "/data_analysis/", experiment_id, "/doublet_cluster/adata.obs.louvain_", (i-1), ".png"), combined_plot, 
         width = 12, height = 18, dpi = 300)
  
  dat[,!colnames(dat) %in% c("umap_1", "umap_2")]
  dat$main_louvain = (i-1)
  
  res = rbind(res, dat)
}

rownames(res) = as.vector(res$cell_id)
res = res[rownames(df_cell),]
df_cell$doublet_cluster = res$doublet_cluster

### sum(df_cell$detected_doublets | df_cell$doublet_cluster) = 19363
### sum(df_cell$detected_doublets | df_cell$doublet_cluster)/nrow(df_cell) = 0.06513169
saveRDS(df_cell, paste0(work_path, "/data_analysis/", experiment_id, "/df_cell.rds"))

### mv doublet_scores_observed_cells_*.csv doublet_cluster/
### mv doublet_scores_simulated_doublets_*.csv doublet_cluster/
### rm df_cell_*.csv gene_count_*.mtx df_gene.csv



