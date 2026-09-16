
#######################################################
### Do sibling cells share the same cell type identity?

source("~/work/scripts/utils.R")
library(dplyr)
library(tidyr)
library(ape)
library(ggplot2)
library(ggrepel)

work_path <- "/net/shendure/vol2/projects/cxqiu/work/tapemouse"

major_trajectory_celltype_table = read.table(paste0(work_path, "/tree_analysis/major_trajectory_celltype_table.txt"), header=T, sep="\t")

# ---- Function 1: extract sibling pairs and annotate ----
analyze_clade <- function(tree, cell_meta, clade_tips, clade_name) {
  subtree <- keep.tip(tree, clade_tips)
  
  tip_parents <- data.frame(
    tip_id   = 1:Ntip(subtree),
    tip_name = subtree$tip.label,
    parent   = subtree$edge[match(1:Ntip(subtree), subtree$edge[, 2]), 1]
  )
  
  sibling_pairs <- tip_parents %>%
    group_by(parent) %>%
    filter(n() >= 2) %>%
    do(as.data.frame(t(combn(.$tip_name, 2)))) %>%
    ungroup() %>%
    rename(cell_A = V1, cell_B = V2)
  
  cat(clade_name, ":", Ntip(subtree), "tips,", nrow(sibling_pairs), "sibling pairs\n")
  
  sibling_celltypes <- sibling_pairs %>%
    left_join(cell_meta %>% select(cell_id, celltype), by = c("cell_A" = "cell_id")) %>%
    rename(celltype_A = celltype) %>%
    left_join(cell_meta %>% select(cell_id, celltype), by = c("cell_B" = "cell_id")) %>%
    rename(celltype_B = celltype)
  
  list(subtree = subtree,
       sibling_pairs = sibling_pairs,
       sibling_celltypes = sibling_celltypes)
}

# ---- Function 2: compute fold change (sibling enrichment) ----
compute_fold_change <- function(subtree, sibling_pairs, sibling_celltypes,
                                cell_meta, strata = NULL,
                                n_perm = 100, seed = 1) {
  set.seed(seed)
  all_tips <- subtree$tip.label

  tip_celltype <- cell_meta$celltype[match(all_tips, cell_meta$cell_id)]
  names(tip_celltype) <- all_tips

  if (is.null(strata)) {
    tip_strata <- rep("all", length(all_tips))
  } else {
    tip_strata <- as.character(strata[[2]][match(all_tips, strata[[1]])])
    tip_strata[is.na(tip_strata)] <- "__unassigned__"
  }
  strata_idx <- split(seq_along(all_tips), tip_strata)

  obs_counts <- sibling_celltypes %>%
    filter(celltype_A == celltype_B, !is.na(celltype_A)) %>%
    count(celltype_A, name = "n_obs") %>%
    rename(celltype = celltype_A)

  null_counts_list <- replicate(n_perm, {
    shuffled_ct <- tip_celltype
    for (idx in strata_idx) {
      shuffled_ct[idx] <- tip_celltype[idx[sample(length(idx))]]
    }
    names(shuffled_ct) <- all_tips
    ct_a <- shuffled_ct[sibling_pairs$cell_A]
    ct_b <- shuffled_ct[sibling_pairs$cell_B]
    table(ct_a[ct_a == ct_b & !is.na(ct_a) & !is.na(ct_b)])
  }, simplify = FALSE)

  all_celltypes <- unique(cell_meta$celltype)
  null_mean <- sapply(all_celltypes, function(ct) {
    mean(sapply(null_counts_list, function(x)
      if (ct %in% names(x)) as.numeric(x[ct]) else 0))
  })

  data.frame(celltype = all_celltypes, n_null_mean = null_mean) %>%
    left_join(obs_counts, by = "celltype") %>%
    mutate(n_obs = ifelse(is.na(n_obs), 0, n_obs),
           fold_change = (n_obs + 1) / (n_null_mean + 1),
           log2_fc = log2(fold_change)) %>%
    arrange(desc(log2_fc))
}



##############################
### Step-1: Major trajectories

cell_meta <- read.table(paste0(work_path, "/tree_analysis/cell_metadata.v8.txt"),
                        header = TRUE, sep = "\t")


cell_meta <- cell_meta %>% select(cell_id, celltype = major_trajectory) %>% as.data.frame()



#################
### BACKBONE TREE

tree <- read.tree(paste0(work_path, "/tree_analysis/merged_minB2h_lineage_constrained.nwk"))

B1 <- read.table(paste0(work_path, "/tree_analysis/e3v8.B1_tape_consensus.tsv.gz"), header = TRUE)
B2 <- read.table(paste0(work_path, "/tree_analysis/e3v8.B2_tape_consensus.tsv.gz"), header = TRUE)

tree_tips_B1 <- tree$tip.label[tree$tip.label %in% B1$cell_id]
tree_tips_B2 <- tree$tip.label[tree$tip.label %in% B2$cell_id]

res_all <- analyze_clade(tree, cell_meta, c(tree_tips_B1, tree_tips_B2), "All")
fc_all <- compute_fold_change(res_all$subtree, res_all$sibling_pairs, res_all$sibling_celltypes, cell_meta, strata = NULL)
saveRDS(fc_all, paste0(work_path, "/tree_analysis/sibling_cells/fc_major_trajectory_backbone_tree_reshuff_global_all.rds"))

### summary
sum(res_all$sibling_celltypes$celltype_A == res_all$sibling_celltypes$celltype_B)/nrow(res_all$sibling_celltypes)
### 64% of tree-sibling pairs share a cell type
sum(fc_all$n_obs)/sum(fc_all$n_null_mean)
### an 3.3-fold enrichment over a permuted null.

### summary (strict-two-tip-cherry)
sum(res_all$sibling_celltypes$celltype_A == res_all$sibling_celltypes$celltype_B)/nrow(res_all$sibling_celltypes)
### 64% of tree-sibling pairs share a cell type
sum(fc_all$n_obs)/sum(fc_all$n_null_mean)
### an 3.2-fold enrichment over a permuted null.

### correlation between full tree and backbone tree
fc_all_full = readRDS(paste0(work_path, "/tree_analysis/sibling_cells/fc_major_trajectory_full_tree_reshuff_global_all.rds"))
fc_all_backbone = readRDS(paste0(work_path, "/tree_analysis/sibling_cells/fc_major_trajectory_backbone_tree_reshuff_global_all.rds"))

df = fc_all_full %>% select(celltype, full_log2_fc = fold_change) %>%
    left_join(fc_all_backbone %>% select(celltype, backbone_log2_fc = fold_change))

cor_test <- cor.test(df$full_log2_fc, df$backbone_log2_fc,
                     method = "spearman")
cat("\nSpearman correlation between full and backbone:", round(cor_test$estimate, 3),
    ", p =", format.pval(cor_test$p.value, digits = 3), "\n")

### Spearman correlation between full and backbone: 0.883 , p = 2.22e-06



##################################################################
### Repeat this analysis, but using a clade-restricted permutation

split = castor::split_tree_at_height(tree, height = 7)

tip_table <- data.frame(
  tip   = tree$tip.label,
  clade = split$clade2subtree[seq_along(tree$tip.label)]
)
### 3848 clones clades at E7.0, of mean size 170 and median size 32

res_all <- analyze_clade(tree, cell_meta, c(tree_tips_B1, tree_tips_B2), "All")
fc_all <- compute_fold_change(res_all$subtree, res_all$sibling_pairs, res_all$sibling_celltypes, cell_meta, strata = tip_table)
saveRDS(fc_all, paste0(work_path, "/tree_analysis/sibling_cells/fc_major_trajectory_backbone_tree_reshuff_clade_all.rds"))


sum(fc_all$n_obs)/sum(fc_all$n_null_mean)
### an 1.9-fold enrichment over a permuted null.










spearman_exact_p <- function(x, y) {
  ok  <- complete.cases(x, y)
  x   <- x[ok]; y <- y[ok]
  n   <- length(x)
  fit <- suppressWarnings(cor.test(x, y, method = "spearman"))
  rho <- unname(fit$estimate)

  t     <- rho * sqrt((n - 2) / (1 - rho^2))
  log10p <- (log(2) + pt(abs(t), df = n - 2, lower.tail = FALSE, log.p = TRUE)) / log(10)

  list(n       = n,
       rho     = rho,
       p       = fit$p.value,          # 0 if underflowed
       log10_p = log10p,               # always usable
       label   = sprintf("rho = %.2f, p ~ 1e%.0f", rho, log10p))
}



######################
### Step-2: Cell types

cell_meta <- read.table(paste0(work_path, "/tree_analysis/cell_metadata.v8.txt"),
                        header = TRUE, sep = "\t")



#################
### BACKBONE TREE

tree <- read.tree(paste0(work_path, "/tree_analysis/merged_minB2h_lineage_constrained.nwk"))

B1 <- read.table(paste0(work_path, "/tree_analysis/e3v8.B1_tape_consensus.tsv.gz"), header = TRUE)
B2 <- read.table(paste0(work_path, "/tree_analysis/e3v8.B2_tape_consensus.tsv.gz"), header = TRUE)

tree_tips_B1 <- tree$tip.label[tree$tip.label %in% B1$cell_id]
tree_tips_B2 <- tree$tip.label[tree$tip.label %in% B2$cell_id]

res_all <- analyze_clade(tree, cell_meta, c(tree_tips_B1, tree_tips_B2), "All")
res_B1 <- analyze_clade(tree, cell_meta, tree_tips_B1, "B1")
res_B2 <- analyze_clade(tree, cell_meta, tree_tips_B2, "B2")

fc_all <- compute_fold_change(res_all$subtree, res_all$sibling_pairs, res_all$sibling_celltypes, cell_meta, strata = NULL)
fc_B1 <- compute_fold_change(res_B1$subtree, res_B1$sibling_pairs, res_B1$sibling_celltypes, cell_meta, strata = NULL)
fc_B2 <- compute_fold_change(res_B2$subtree, res_B2$sibling_pairs, res_B2$sibling_celltypes, cell_meta, strata = NULL)


saveRDS(fc_all, paste0(work_path, "/tree_analysis/sibling_cells/fc_celltype_backbone_tree_reshuff_global_all.rds"))
saveRDS(fc_B1, paste0(work_path, "/tree_analysis/sibling_cells/fc_celltype_backbone_tree_reshuff_global_B1.rds"))
saveRDS(fc_B2, paste0(work_path, "/tree_analysis/sibling_cells/fc_celltype_backbone_tree_reshuff_global_B2.rds"))


### summary
sum(res_all$sibling_celltypes$celltype_A == res_all$sibling_celltypes$celltype_B)/nrow(res_all$sibling_celltypes)
### 41% of tree-sibling pairs share a cell type
sum(fc_all$n_obs)/sum(fc_all$n_null_mean)
### an 7.9-fold enrichment over a permuted null.



### summary (strict-two-tip-cherry)
sum(res_all$sibling_celltypes$celltype_A == res_all$sibling_celltypes$celltype_B)/nrow(res_all$sibling_celltypes)
### 40% of tree-sibling pairs share a cell type
sum(fc_all$n_obs)/sum(fc_all$n_null_mean)
### an 7.7-fold enrichment over a permuted null.


x = read.csv(paste0(work_path, "/tree_analysis/sibling_cells/fc_celltype_full_tree_reshuff_global.csv"))
celltype_common = unique(x$celltype)

table_out = rbind(fc_B1 %>% filter(celltype %in% celltype_common) %>% mutate(blastomere = "Blastomere_A"),
  fc_B2 %>% filter(celltype %in% celltype_common) %>% mutate(blastomere = "Blastomere_B"),
  fc_all %>% filter(celltype %in% celltype_common) %>% mutate(blastomere = "Both"))

write.csv(table_out, paste0(work_path, "/tree_analysis/sibling_cells/fc_celltype_backbone_tree_reshuff_global.csv"), row.names=F, quote=F)



### correlation between full tree and backbone tree
fc_all_full = readRDS(paste0(work_path, "/tree_analysis/sibling_cells/fc_celltype_full_tree_reshuff_global_all.rds"))
fc_all_backbone = readRDS(paste0(work_path, "/tree_analysis/sibling_cells/fc_celltype_backbone_tree_reshuff_global_all.rds"))

df = fc_all_full %>% select(celltype, full_log2_fc = fold_change) %>%
    left_join(fc_all_backbone %>% select(celltype, backbone_log2_fc = fold_change))

cor_test <- cor.test(df$full_log2_fc, df$backbone_log2_fc,
                     method = "spearman")
cat("\nSpearman correlation between full and backbone:", round(cor_test$estimate, 3),
    ", p =", format.pval(cor_test$p.value, digits = 3), "\n")

### Spearman correlation between full and backbone: 0.943 , p = <2e-16



##################################################################
### Repeat this analysis, but using a clade-restricted permutation

split = castor::split_tree_at_height(tree, height = 7)

tip_table <- data.frame(
  tip   = tree$tip.label,
  clade = split$clade2subtree[seq_along(tree$tip.label)]
)
### 3848 clones clades at E7.0, of mean size 170 and median size 32


fc_all <- compute_fold_change(res_all$subtree, res_all$sibling_pairs, res_all$sibling_celltypes, cell_meta, strata = tip_table)
fc_B1 <- compute_fold_change(res_B1$subtree, res_B1$sibling_pairs, res_B1$sibling_celltypes, cell_meta, strata = tip_table)
fc_B2 <- compute_fold_change(res_B2$subtree, res_B2$sibling_pairs, res_B2$sibling_celltypes, cell_meta, strata = tip_table)


saveRDS(fc_all, paste0(work_path, "/tree_analysis/sibling_cells/fc_celltype_backbone_tree_reshuff_clade_all.rds"))
saveRDS(fc_B1, paste0(work_path, "/tree_analysis/sibling_cells/fc_celltype_backbone_tree_reshuff_clade_B1.rds"))
saveRDS(fc_B2, paste0(work_path, "/tree_analysis/sibling_cells/fc_celltype_backbone_tree_reshuff_clade_B2.rds"))


sum(fc_all$n_obs)/sum(fc_all$n_null_mean)
### an 2.8-fold enrichment over a permuted null.

x = read.csv(paste0(work_path, "/tree_analysis/sibling_cells/fc_celltype_full_tree_reshuff_global.csv"))
celltype_common = unique(x$celltype)


table_out = rbind(fc_B1 %>% filter(celltype %in% celltype_common) %>% mutate(blastomere = "Blastomere_A"),
  fc_B2 %>% filter(celltype %in% celltype_common) %>% mutate(blastomere = "Blastomere_B"),
  fc_all %>% filter(celltype %in% celltype_common) %>% mutate(blastomere = "Both"))

write.csv(table_out, paste0(work_path, "/tree_analysis/sibling_cells/fc_celltype_backbone_tree_reshuff_clade.csv"), row.names=F, quote=F)


setwd('/net/shendure/vol2/projects/cxqiu/work/tapemouse/tree_analysis/sibling_cells')
x1 = read.csv('fc_celltype_full_tree_reshuff_global.csv')
x2 = read.csv('fc_celltype_full_tree_reshuff_clade.csv')
x3 = read.csv('fc_celltype_backbone_tree_reshuff_global.csv')
x4 = read.csv('fc_celltype_backbone_tree_reshuff_clade.csv')

x = x1 %>% left_join(x2, by = c("celltype", "blastomere")) %>%
left_join(x3, by = c("celltype", "blastomere")) %>%
left_join(x4, by = c("celltype", "blastomere"))

write.csv(x, paste0(work_path, "/tree_analysis/sibling_cells/fc_celltype.csv"), row.names=F, quote=F)



x = res_all$sibling_celltypes %>% filter(celltype_A == celltype_B) %>%
    group_by(celltype_A) %>% tally() %>% filter(n >= 5)
# n = 88 cell types with at least five same-type sibling pairs

fc_all_sub = fc_all %>% filter(celltype %in% x$celltype_A)

> summary(fc_all_sub$fold_change)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  1.890   4.595   6.245   7.484   9.027  30.769

# enrichments remained above 1.0 for all 88 cell types with at least five same-type sibling pairs


fc_all_reshuff_global = readRDS(paste0(work_path, "/tree_analysis/sibling_cells/fc_celltype_backbone_tree_reshuff_global_all.rds"))
fc_all_reshuff_clade = readRDS(paste0(work_path, "/tree_analysis/sibling_cells/fc_celltype_backbone_tree_reshuff_clade_all.rds"))

fc_all_reshuff_global[fc_all_reshuff_global$celltype == "Melanocyte cells",]
fc_all_reshuff_clade[fc_all_reshuff_clade$celltype == "Melanocyte cells",]
### 47 -> 31

fc_all_reshuff_global[fc_all_reshuff_global$celltype == "Dorsal telencephalon",]
fc_all_reshuff_clade[fc_all_reshuff_clade$celltype == "Dorsal telencephalon",]
### 50 -> 18

df = fc_all_reshuff_global %>% select(celltype, fc_global = fold_change) %>%
    inner_join(fc_all_reshuff_clade %>% select(celltype, fc_clade = fold_change), by = "celltype") %>%
    filter(celltype %in% x$celltype_A)

cor.test(df$fc_global, df$fc_clade, method = "spearman")
# Spearman rho = 0.45; p = 1e-5


