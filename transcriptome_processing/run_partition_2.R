
########################################################################################
### we take the top200 genes of each partition to perform subclustering on each partition

work_path = "/net/shendure/vol2/projects/cxqiu/work/tapemouse"
source("~/work/scripts/utils.R")

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


