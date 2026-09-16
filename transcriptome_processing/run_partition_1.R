
### This step might not be necessary, but still worth to take a look if additional doublets

work_path = "/net/shendure/vol2/projects/cxqiu/work/tapemouse"
source("~/work/scripts/utils.R")

args = commandArgs(trailingOnly=TRUE)
experiment_id = args[1]

batch_num = 8

folder = paste0(work_path, "/data_analysis/", experiment_id, '/doublet_cluster_2')
if (!dir.exists(folder)) {dir.create(folder, recursive = TRUE)}

fd = readRDS(paste0(work_path, "/df_gene.rds"))
rownames(fd) = fd$gene_ID
fd = fd[(fd$gene_type %in% c('protein_coding', 'lncRNA')) & (!fd$chr %in% c('chrX', 'chrY')),]

count_all = NULL
pd_all = NULL
for(i in 1:batch_num){
    count = readRDS(paste0(work_path, "/data_analysis/", experiment_id, "/gene_count_", i, ".rds"))
    pd = readRDS(paste0(work_path, "/data_analysis/", experiment_id, "/df_cell.rds"))
    pd = pd[colnames(count),]
    pd = pd[!(pd$detected_doublets | pd$doublet_cluster),]
    count = count[rownames(fd), rownames(pd)]
    print(paste0(i, "/", dim(count)))
    pd_all = rbind(pd_all, pd)
    count_all = cbind(count_all, count)
}
count = count_all
pd = pd_all

cell_keep = pd$gene_count >= 100
count = count[,cell_keep]
pd = pd[cell_keep,]

count_binary = count
count_binary@x[count_binary@x > 0] = 1
gene_keep = rowSums(count_binary) >= 10
count = count[gene_keep,]

obj = CreateSeuratObject(count, meta.data = pd)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2500)
obj = ScaleData(object = obj, verbose = FALSE)
obj = RunPCA(object = obj, npcs = 30, verbose = FALSE)
obj = RunUMAP(object = obj, reduction = "pca", dims = 1:30, min.dist = 0.01, n.neighbors = 50, n.components = 2)
res = doIdentifyPartitions(Embeddings(obj, reduction = "umap"), data.frame(obj[[]]))
obj$my_partition = res[[1]]
obj$my_cluster = res[[2]]
print(table(obj$my_partition))

p = DimPlot(obj, group.by = "my_partition", label = T, raster=FALSE) + NoLegend()
ggsave(paste0(work_path, "/data_analysis/", experiment_id, '/doublet_cluster_2/partition.png'), p, dpi = 300)

saveRDS(obj, paste0(work_path, "/data_analysis/", experiment_id, '/doublet_cluster_2/obj_all.rds'))


#####################################################################
### calculating top200 differential expressed genes for each partition

Idents(obj) = as.vector(obj$my_partition)
obj_sub = subset(obj, downsample = 2500)
DEG = FindAllMarkers(obj_sub, only.pos = T)
top_DEG = DEG %>% group_by(cluster) %>% 
    filter(p_val_adj < 0.05) %>%
    slice_max(order_by = avg_log2FC, n = 200)
saveRDS(top_DEG, paste0(work_path, "/data_analysis/", experiment_id, '/doublet_cluster_2/top_DEG.rds'))

