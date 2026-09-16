
################################################################################
### removing potential doublets based on the top marker gene of each cluster ###
################################################################################

work_path = "/net/shendure/vol2/projects/cxqiu/work/tapemouse"
source("~/work/scripts/utils.R")

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


