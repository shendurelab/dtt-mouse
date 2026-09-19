
#######################################################################
### Compare gene expression between wildtype E13.5 and tapemouse E13.5

### Datasets

### 1) Transcriptome data of emrbyo #3 in this study:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/tapemouse/adata.{experiment_id}.h5a

### 2) Transcriptome data of single-cell atlas of mouse development:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax/mm39_version/adata.{day_id}.h5ad
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax/mm39_version/pd.rds
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax/mm39_version/df_gene.rds

### 3) Other support data:
### https://github.com/shendurelab/dtt-mouse/tree/main/support_data
### mouse.v37.geneID.txt
### cell_metadata.v8.txt
### insertion_sites.txt


mouse_gene = read.table("mouse.v37.geneID.txt", header=T, sep="\t")

pd_all = readRDS(paste0(work_path, "/adata_integration.obs.rds"))

pd_jax = readRDS(paste0(work_path, "/pd.rds"))
pd_jax = pd_jax[,c("cell_id", "UMI_count", "gene_count")]
pd_tapemouse = read.table(paste0(work_path, "/cell_metadata.v8.txt"), header=T)
pd_tapemouse = pd_tapemouse[,c("cell_id", "UMI_count", "gene_count")]
pd_x = rbind(pd_jax, pd_tapemouse)

pd_all = pd_all %>% filter(day == "E13.5") %>% left_join(pd_x, by = "cell_id")

set.seed(1)

lo <- max(min(pd_all$UMI_count[pd_all$dataset == "jax"]),
          min(pd_all$UMI_count[pd_all$dataset == "tapemouse"]))
hi <- min(max(pd_all$UMI_count[pd_all$dataset == "jax"]),
          max(pd_all$UMI_count[pd_all$dataset == "tapemouse"]))

pd_clip <- pd_all %>% filter(UMI_count >= lo, UMI_count <= hi)

pd_clip$bin <- cut(pd_clip$UMI_count, breaks = 40, include.lowest = TRUE)

pd_all_sub <- pd_clip %>%
  group_by(bin) %>%
  filter(n_distinct(dataset) == 2) %>%
  group_by(bin, dataset) %>%
  mutate(bin_n = n()) %>%
  group_by(bin) %>%
  mutate(n_take = min(bin_n)) %>%
  group_by(bin, dataset) %>%
  group_modify(~ slice_sample(.x, n = .x$n_take[1])) %>%
  ungroup() %>%
  select(-bin, -bin_n, -n_take)

pd_all_sub %>% group_by(dataset) %>%
  summarise(median_umi = median(UMI_count),
            median_gene = median(gene_count),
            n = n())

summary(pd_all_sub$UMI_count[pd_all_sub$dataset == "jax"])
summary(pd_all_sub$UMI_count[pd_all_sub$dataset == "tapemouse"])

summary(pd_all_sub$gene_count[pd_all_sub$dataset == "jax"])
summary(pd_all_sub$gene_count[pd_all_sub$dataset == "tapemouse"])


   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    803    1439    1988    2189    2769    5577
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    803    1438    1988    2187    2766    5577
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    306     980    1263    1325    1614    3284
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    478     943    1208    1271    1543    3180

pd_all_sub$log2_umi = log2(pd_all_sub$UMI_count)

p = ggplot(pd_all_sub, aes(log2_umi)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  facet_wrap(~ dataset, ncol = 1) +
  labs(x = "Log2(UMI count)", y = "# of cells") +
  theme_classic()
ggsave("log2_umi_count.pdf", p, width = 4, height = 5)


pd_tmp = pd_all_sub %>%
    group_by(dataset, celltype) %>% tally() %>%
    filter(n >= 500)

common_celltype = intersect(pd_tmp$celltype[pd_tmp$dataset == "jax"],
                            pd_tmp$celltype[pd_tmp$dataset == "tapemouse"])

set.seed(1234)
pd = pd_all_sub %>% filter(celltype %in% common_celltype) %>%
    group_by(dataset, celltype) %>% slice_sample(n = 500) %>% as.data.frame()

write.csv(pd[,c("cell_id","dataset","celltype")], paste0(work_path, "/pd_E13.5_500_cells_per_celltype.csv"), row.names=F)


### subset the gene expression data

experiment_list = c("experiment1_20260618_seq4_AD", 
                    "experiment1_20260618_seq4_EH",
                    "experiment1_20260618_seq5_IL",
                    "experiment1_20260618_seq5_MP",
                    "experiment1_20260618_seq6_QT",
                    "experiment1_20260618_seq6_UW",
                    "experiment2_20260713_seq2_XY")

df_gene = readRDS(paste0(work_path, "/df_gene.rds"))
df_gene_sub = df_gene[df_gene$gene_type %in% c("protein_coding") & df_gene$chr %in% paste0('chr', 1:19),]

gene_count = NULL
for(experiment_id in experiment_list){
  print(experiment_id)
  obj_i = readRDS(paste0(work_path, "/obj.rds"))
  gene_count_i = GetAssayData(obj_i, slot = "counts")
  gene_count_i = gene_count_i[df_gene_sub$gene_ID, colnames(gene_count_i) %in% pd$cell_id]
  gene_count = cbind(gene_count, gene_count_i)
}

gene_count_i = readRDS(paste0(work_path, "/gene_count_E13.5.rds"))
gene_count_i = gene_count_i[df_gene_sub$gene_ID, colnames(gene_count_i) %in% pd$cell_id]

gene_count = cbind(gene_count, gene_count_i)
rownames(pd) = pd$cell_id
pd = pd[colnames(gene_count),]

obj = CreateSeuratObject(gene_count, meta.data = pd)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
Idents(obj) = obj$dataset

saveRDS(obj, paste0(work_path, "/obj.rds"))

summary(obj$nCount_RNA[obj$dataset == "jax"])
summary(obj$nCount_RNA[obj$dataset == "tapemouse"])

summary(obj$nFeature_RNA[obj$dataset == "jax"])
summary(obj$nFeature_RNA[obj$dataset == "tapemouse"])

   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    647    1267    1760    1923    2439    4970
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    687    1208    1649    1821    2266    4986
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    459     890    1144    1195    1454    2704
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    487     812    1030    1093    1313    2634

res_list = list()
for(celltype_i in common_celltype){
    print(celltype_i)
    obj_i = subset(obj, subset = celltype == celltype_i)
    res_i = FindMarkers(obj_i, ident.1 = "jax", ident.2 = "tapemouse") %>%
                     tibble::rownames_to_column("gene_ID") %>%
                     mutate(celltype = celltype_i, fdr = p.adjust(p_val, "BH"))
    res_list[[celltype_i]] = res_i
}

res = do.call(rbind, res_list)
rownames(res) = NULL
res = res %>% left_join(mouse_gene[,c("gene_ID", "gene_short_name")], by = "gene_ID")
res$up_down = if_else(res$avg_log2FC > 0, "jax", "tapemouse")
saveRDS(res, paste0(work_path, "/DEG.rds"))

res_x = res %>% filter(p_val_adj < 0.05) %>% group_by(gene_short_name, up_down) %>% tally() %>% arrange(-n) %>% filter(n >= 10)

res_x$gene_short_name = factor(res_x$gene_short_name, levels = as.vector(res_x$gene_short_name))

p = res_x %>%
  ggplot(aes(gene_short_name, n, fill = up_down)) +
  geom_col() +
  scale_fill_manual(values = c(jax = "firebrick", tapemouse = "steelblue")) +
  labs(x = NULL, y = "# cell types", fill = "Higher in") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90))
ggsave("DE.pdf", p, width = 6, height = 5)



read_gmt <- function(path) {
  lines <- readLines(path)
  setNames(
    lapply(strsplit(lines, "\t"), function(x) x[-c(1, 2)]),
    sapply(strsplit(lines, "\t"), `[`, 1)
  )
}
programs = read_gmt(paste0(work_path, "/gene_sets_mouse.gmt"))

programs_update = list()
for(program_i in names(programs)){
    print(program_i)
    programs_update[[program_i]] = mouse_gene$gene_ID[mouse_gene$gene_short_name %in% programs[[program_i]] &
               mouse_gene$gene_ID %in% rownames(obj)]
}


obj <- AddModuleScore(obj, features = programs_update, name = "MS_")
ms_cols <- paste0("MS_", seq_along(programs_update))
new_cols <- paste0(names(programs_update), "_MS")
colnames(obj@meta.data)[match(ms_cols, colnames(obj@meta.data))] <- new_cols

saveRDS(obj, paste0(work_path, "/obj.rds"))

res_list <- list()
for (celltype_i in common_celltype) {
  print(celltype_i)
  df <- obj@meta.data %>%
    filter(celltype == celltype_i) %>%
    select(dataset, ends_with("_MS"))
  
  res_i <- lapply(names(programs_update), function(p) {
    col  <- paste0(p, "_MS")
    jax  <- df[[col]][df$dataset == "jax"]
    tape <- df[[col]][df$dataset == "tapemouse"]
    data.frame(
      celltype = celltype_i,
      program  = p,
      diff     = mean(jax) - mean(tape),
      cohens_d = (mean(jax) - mean(tape)) / sqrt((var(jax) + var(tape)) / 2),
      p_val    = wilcox.test(jax, tape)$p.value
    )
  }) %>% bind_rows()
  
  res_list[[celltype_i]] <- res_i
}

res <- do.call(rbind, res_list)
rownames(res) <- NULL
res <- res %>%
  group_by(program) %>%
  mutate(p_val_adj = p.adjust(p_val, "BH")) %>%
  ungroup() %>%
  mutate(up_down = if_else(cohens_d > 0, "jax", "tapemouse"))
saveRDS(res, paste0(work_path, "/DEpathway.rds"))

res_x <- res %>%
  filter(p_val_adj < 0.05) %>%
  group_by(program, up_down) %>%
  tally() %>%
  arrange(-n)


res_x$program_x = paste0(res_x$up_down, "_", res_x$program)
res_x$program_x = factor(res_x$program_x, levels = unique(res_x$program_x))
p = res_x %>%
  ggplot(aes(program_x, n, fill = up_down)) +
  geom_col() +
  scale_fill_manual(values = c(jax = "firebrick", tapemouse = "steelblue")) +
  scale_x_discrete(labels = setNames(res_x$program, res_x$program_x)) +
  labs(x = NULL, y = "# cell types", fill = "Higher in") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("DEpathway.pdf", p, width = 6, height = 5)



####################################################################################
### Step-2: checking candidate genes which are overlapped with TapeWriter insertions


mouse_gene = read.table("mouse.v37.geneID.txt", header=T, sep="\t")

insertion_sites = read.table(paste0(work_path, "/insertion_sites.txt"))
ins = insertion_sites %>%
  separate(V1, into = c("chr", "start", "end"), sep = "[:-]", convert = TRUE, remove = FALSE) %>%
  mutate(center = (start + end) %/% 2)

mouse_gene = mouse_gene %>%
  mutate(TSS = ifelse(strand == "+", start, end))

window = 50000
hits = mouse_gene %>%
  inner_join(ins, by = "chr", relationship = "many-to-many", suffix = c(".gene", ".ins")) %>%
  mutate(
    region_start = pmin(TSS - window, start.gene),
    region_end   = pmax(TSS + window, end.gene)
  ) %>%
  filter(center >= region_start & center <= region_end) %>%
  mutate(distance = ifelse(strand == "+", center - TSS, TSS - center)) %>%
  select(gene_ID, gene_short_name, gene_type, chr, TSS, strand,
         insertion = V1, distance) %>%
  arrange(chr, TSS)

gene_list = unique(hits$gene_ID)
### n = 91 genes

pd = read.csv(paste0(work_path, "/pd_E13.5_500_cells_per_celltype.csv"))
common_celltype = unique(pd$celltype)

experiment_list = c("experiment1_20260618_seq4_AD", 
                    "experiment1_20260618_seq4_EH",
                    "experiment1_20260618_seq5_IL",
                    "experiment1_20260618_seq5_MP",
                    "experiment1_20260618_seq6_QT",
                    "experiment1_20260618_seq6_UW",
                    "experiment2_20260713_seq2_XY")

gene_count = NULL
for(experiment_id in experiment_list){
  print(experiment_id)
  obj_i = readRDS(paste0(work_path, "/data_analysis/", experiment_id, "/obj.rds"))
  gene_count_i = GetAssayData(obj_i, slot = "counts")
  gene_count_i = gene_count_i[,colnames(gene_count_i) %in% pd$cell_id]
  gene_count = cbind(gene_count, gene_count_i)
}

gene_count_i = readRDS(paste0(work_path, "/gene_count_E13.5.rds"))
gene_count_i = gene_count_i[,colnames(gene_count_i) %in% pd$cell_id]

gene_count = cbind(gene_count, gene_count_i)

rownames(pd) = pd$cell_id
pd = pd[colnames(gene_count),]

obj = CreateSeuratObject(gene_count, meta.data = pd)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
Idents(obj) = obj$dataset

saveRDS(obj, paste0(work_path, "/obj_full.rds"))


res_list = list()
for(celltype_i in common_celltype){
    print(celltype_i)
    obj_i = subset(obj, subset = celltype == celltype_i)
    res_i = FindMarkers(obj_i, ident.1 = "jax", ident.2 = "tapemouse", logfc.threshold = 0, min.pct = 0, features = gene_list) %>%
                     tibble::rownames_to_column("gene_ID") %>%
                     mutate(celltype = celltype_i, fdr = p.adjust(p_val, "BH"))
    res_list[[celltype_i]] = res_i
}

res = do.call(rbind, res_list)
rownames(res) = NULL

res = res %>% left_join(mouse_gene[,c("gene_ID", "gene_short_name")], by = "gene_ID")
res$up_down = if_else(res$avg_log2FC > 0, "jax", "tapemouse")
saveRDS(res, paste0(work_path, "/DEG_genes_overlap_insertions.rds"))

res_sig = res %>% filter(fdr < 0.05) %>% 
    group_by(gene_short_name, up_down) %>% tally() %>% arrange(-n)

res_sig$gene_short_name_x = paste0(res_sig$up_down, "_", res_sig$gene_short_name)
res_sig$gene_short_name_x = factor(res_sig$gene_short_name_x, levels = unique(res_sig$gene_short_name_x))
p = res_sig %>%
  ggplot(aes(gene_short_name_x, n, fill = up_down)) +
  geom_col() +
  scale_fill_manual(values = c(jax = "firebrick", tapemouse = "steelblue")) +
  scale_x_discrete(labels = setNames(res_sig$gene_short_name, res_sig$gene_short_name_x)) +
  labs(x = NULL, y = "# cell types", fill = "Higher in") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("DEG_genes_overlap_insertions.pdf", p, width = 6, height = 5)


