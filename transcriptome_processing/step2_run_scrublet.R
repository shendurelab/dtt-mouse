
###########################################################################################
### This step might not be necessary, but still worth to take a look if additional doublets

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



########################################################################################
### we take the top200 genes of each partition to perform subclustering on each partition


args = commandArgs(trailingOnly=TRUE)
experiment_id = args[2]

obj_all = readRDS(paste0(work_path, "/data_analysis/", experiment_id, '/doublet_cluster_2/obj_all.rds'))

cluster_list = c(1:max(obj_all$my_partition))

fd = readRDS(paste0(work_path, "/df_gene.rds"))
rownames(fd) = fd$gene_ID
fd = fd[rownames(obj_all),]

top_DEG = readRDS(paste0(work_path, "/data_analysis/", experiment_id, '/doublet_cluster_2/top_DEG.rds'))
gene_used = as.vector(unique(top_DEG$gene))

cluster_i = cluster_list[as.numeric(args[1])]

    print(cluster_i)
    obj = subset(obj_all, subset = my_partition == cluster_i)
    count = GetAssayData(obj, slot = "counts"); meta_data = data.frame(obj[[]])
    
    if(ncol(count) <= 30){
        saveRDS(meta_data, paste0(work_path, "/data_analysis/", experiment_id, '/doublet_cluster_2/cluster_', cluster_i, '_pd.rds'))
    } else {
        obj = CreateSeuratObject(count, meta.data = meta_data)
        obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
        obj = ScaleData(object = obj, features = gene_used, verbose = FALSE)
        obj = RunPCA(object = obj, npcs = 30, features = gene_used, verbose = FALSE)
        obj = RunUMAP(object = obj, reduction = "pca", dims = 1:30, min.dist = 0.1, n.neighbors = 30, n.components = 2)
        obj = FindNeighbors(object = obj, dims = 1:2, reduction = "umap")
        obj = FindClusters(object = obj, resolution = 0.3)
        pd = data.frame(obj[[]])
        
        p1 = DimPlot(obj, group.by = "seurat_clusters", label = T, raster=FALSE, label.size = 3) + NoLegend()

        cds = new_cell_data_set(count, cell_metadata = meta_data, gene_metadata = fd)
        reducedDims(cds)$UMAP = as.matrix(Embeddings(obj, reduction = "umap"))
        p2 = plot_cells(cds, color_cells_by = "doublet_score", cell_size = 0.6) + NoLegend()
        
        medians = pd %>% group_by(seurat_clusters) %>% summarize(Q1 = quantile(doublet_score, 0.25), Q2 = quantile(doublet_score, 0.5), Q3 = quantile(doublet_score, 0.75)) %>%
            mutate(fill_color = ifelse(Q1 > 0.025 | Q2 > 0.05 | Q3 > 0.075, "red", "gray"))
        pd2 = pd %>% left_join(medians, by = "seurat_clusters")
        
        p3 = ggplot(pd2, aes(seurat_clusters, doublet_score, fill = fill_color)) + geom_boxplot() + 
            coord_flip() + scale_fill_identity() + theme(axis.text.y = element_text(size = 8))
        
        p = (p1 / p2) | p3
        ggsave(paste0(work_path, "/data_analysis/", experiment_id, "/doublet_cluster_2/cluster_", cluster_i, ".png"), p, dpi=300, height=10, width=10)
    
        saveRDS(pd, paste0(work_path, "/data_analysis/", experiment_id, '/doublet_cluster_2/cluster_', cluster_i, '_pd.rds'))
    }





################################################################################
### removing potential doublets based on the top marker gene of each cluster ###
################################################################################



experiment_id = "experiment1_20260618_seq4_AD"

pd_all = NULL

for(i in 1:13){
    print(i)
    exclude_cluster_list = c(-1)
    if(i == 3){exclude_cluster_list = c(13)}
    if(i == 5){exclude_cluster_list = c(8,22,24)}
    if(i == 6){exclude_cluster_list = c(2,8,9)}
    if(i == 7){exclude_cluster_list = c(3,7,8,14:19,21:23)}
    if(i == 8){exclude_cluster_list = c(20)}
    if(i == 9){exclude_cluster_list = c(10,11)}
    if(i == 10){exclude_cluster_list = c(6)}
    if(i == 12){exclude_cluster_list = c(1,2)}

    pd_tmp = readRDS(paste0(work_path, "/data_analysis/", experiment_id, "/doublet_cluster_2/cluster_", i, "_pd.rds"))
    if(!'seurat_clusters' %in% colnames(pd_tmp)) {
        pd_tmp$doublet_DEG = FALSE
    } else {
        pd_tmp$doublet_DEG = pd_tmp$seurat_clusters %in% exclude_cluster_list
    }
    pd_tmp = pd_tmp[,c("cell_id", "doublet_DEG", "doublet_score")]
    pd_all = rbind(pd_all, pd_tmp)
}

saveRDS(pd_all[,c("cell_id", "doublet_DEG")], paste0(work_path, "/data_analysis/", experiment_id,  "/doublet_cluster_2/res_doubelt_DEG.rds"))


