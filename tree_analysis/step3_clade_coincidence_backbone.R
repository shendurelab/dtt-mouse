
########################
### Clade Coincidence

source("~/work/scripts/utils.R")
library(dplyr)
library(tidyr)
library(ape)
library(ggplot2)
library(ggrepel)
library(phangorn)

work_path <- "/net/shendure/vol2/projects/cxqiu/work/tapemouse"

# ---- Load data ----
cell_meta <- read.table(paste0(work_path, "/tree_analysis/cell_metadata.v8.txt"),
                        header = TRUE, sep = "\t")

tree <- read.tree(paste0(work_path, "/tree_analysis/merged_minB2h_lineage_constrained.nwk"))

B1 <- read.table(paste0(work_path, "/tree_analysis/e3v8.B1_tape_consensus.tsv.gz"), header = TRUE)
B2 <- read.table(paste0(work_path, "/tree_analysis/e3v8.B2_tape_consensus.tsv.gz"), header = TRUE)

tree_tips_B1 <- tree$tip.label[tree$tip.label %in% B1$cell_id]
tree_tips_B2 <- tree$tip.label[tree$tip.label %in% B2$cell_id]


# ---- Internal nodes ----
internal_nodes <- function(tree, rescale = 1.5){
  depths <- node.depth.edgelength(tree)

  n_tips <- Ntip(tree)
  n_int  <- Nnode(tree)

  internal <- data.frame(
    node = paste0("node_", (n_tips + 1):(n_tips + n_int)),
    node_depth  = depths[(n_tips + 1):(n_tips + n_int)] + rescale
  )

  return(internal)
}


# ---- Function: cut tree into K clades (memory-efficient) ----
maximal_clades <- function(tree, K = 60, min_size = 3) {
  n_tips  <- Ntip(tree)
  n_nodes <- n_tips + Nnode(tree)
  
  desc <- Descendants(tree, 1:n_nodes, type = "tips")
  size <- lengths(desc)
  
  parent_of <- integer(n_nodes)
  parent_of[tree$edge[, 2]] <- tree$edge[, 1]
  
  clade_roots <- c()
  for (u in (n_tips + 1):n_nodes) {
    if (size[u] < min_size || size[u] > K) next
    if (parent_of[u] != 0 && size[parent_of[u]] <= K) next
    clade_roots <- c(clade_roots, u)
  }
  
  cat("K =", K, ", n maximal clades =", length(clade_roots), "\n")
  
  clade_ids   <- rep(NA_integer_, n_tips)
  clade_nodes <- rep(NA_integer_, n_tips)                       # new
  for (i in seq_along(clade_roots)) {
    clade_ids[desc[[clade_roots[i]]]]   <- i
    clade_nodes[desc[[clade_roots[i]]]] <- clade_roots[i]        # new
  }
  
  in_clade <- !is.na(clade_ids)
  cat("Cells covered:", sum(in_clade), "/", n_tips, "\n")
  
  data.frame(cell_id     = tree$tip.label[in_clade],
             clade       = clade_ids[in_clade],
             node  = clade_nodes[in_clade])                # new
}


# ---- Function: compute observed co-occurrence matrix (clade-based) ----
compute_cooccur_matrix <- function(clade_df, cell_meta, all_celltypes) {
  df <- clade_df %>%
    left_join(cell_meta, by = "cell_id") %>%
    filter(!is.na(celltype), celltype %in% all_celltypes)
  
  df$celltype <- factor(df$celltype, levels = all_celltypes)
  
  # cell type × clade count matrix
  counts <- table(df$celltype, df$clade)      # rows = cell types, cols = clades
  
  # Presence indicator: 1 if a type is in the clade at all
  present <- counts > 0
  
  # Off-diagonal: how many clades contain BOTH type i and type j
  # This is presence[i, ] %*% presence[j, ]^T summed across clades
  obs <- present %*% t(present)
  
  # Diagonal: how many clades contain type i with ≥ 2 cells
  diag(obs) <- rowSums(counts >= 2)
  
  # Return as plain matrix with names
  dimnames(obs) <- list(all_celltypes, all_celltypes)
  as.matrix(obs)
}


# ---- Function: compute expected co-occurrence matrix under random assignment ----
compute_expected_matrix <- function(clade_df, cell_meta, all_celltypes,
                                     n_perm = 50, seed = 1) {
  set.seed(seed)
  df <- clade_df %>%
    left_join(cell_meta, by = "cell_id") %>%
    filter(!is.na(celltype), celltype %in% all_celltypes)
  df$celltype <- factor(df$celltype, levels = all_celltypes)
  
  exp_mat <- matrix(0, nrow = length(all_celltypes), ncol = length(all_celltypes),
                    dimnames = list(all_celltypes, all_celltypes))
  
  for (perm in 1:n_perm) {
    df_shuf <- df
    df_shuf$celltype <- sample(df$celltype)
    
    counts <- table(df_shuf$celltype, df_shuf$clade)
    present <- counts > 0
    
    m <- present %*% t(present)
    diag(m) <- rowSums(counts >= 2)
    
    exp_mat <- exp_mat + m
  }
  
  exp_mat / n_perm
}



# ---- Full pipeline for one blastomere at one K ----
run_clade_analysis <- function(tree, cell_meta, clade_tips, blastomere_name, K,
                                all_celltypes, n_perm = 100) {
  cat(blastomere_name, ", K =", K, "...\n")
  subtree <- keep.tip(tree, clade_tips)
  cat("  Subtree tips:", Ntip(subtree), "\n")
  
  clade_df <- maximal_clades(subtree, K)
  cat("  Cell types (shared):", length(all_celltypes), "\n")

  if(blastomere_name == "All"){
    internal <- internal_nodes(subtree, rescale = 0)
  } else {
    internal <- internal_nodes(subtree, rescale = 1.5)
  }
  
  clade_df$node = paste0("node_", clade_df$node)
  clade_df_x = clade_df %>% left_join(internal, by = "node")
  clade_df$node_depth = clade_df_x$node_depth
  
  obs <- compute_cooccur_matrix(clade_df, cell_meta, all_celltypes)
  exp_mat <- compute_expected_matrix(clade_df, cell_meta, all_celltypes, n_perm)

  enr <- log2((obs + 1) / (exp_mat + 1))
  
  list(obs = obs, exp = exp_mat, log2_enr = enr,
       K = K, blastomere = blastomere_name, clade_df = clade_df)
}


# ---- Run for each blastomere at K = 15 and K = 60 ----

res_B1_K15 = readRDS(paste0(work_path, "/tree_analysis/clade_coincodence/res_B1_K15.rds"))$log2_enr
shared_celltypes = rownames(res_B1_K15)


for (K in c(15, 200)) {
  res <- run_clade_analysis(tree, cell_meta, tree_tips_B1, "Blastomere A", K, shared_celltypes)
  saveRDS(res, paste0(work_path, "/tree_analysis/clade_coincodence/res_B1_K", K, "_backbone.rds"))

  res <- run_clade_analysis(tree, cell_meta, tree_tips_B2, "Blastomere B", K, shared_celltypes)
  saveRDS(res, paste0(work_path, "/tree_analysis/clade_coincodence/res_B2_K", K, "_backbone.rds"))
}



### off-diagonal correlation

res_B1_K15 = readRDS(paste0(work_path, "/tree_analysis/clade_coincodence/res_B1_K15.rds"))$log2_enr
res_B2_K15 = readRDS(paste0(work_path, "/tree_analysis/clade_coincodence/res_B2_K15.rds"))$log2_enr

res_B1_K200 = readRDS(paste0(work_path, "/tree_analysis/clade_coincodence/res_B1_K200.rds"))$log2_enr
res_B2_K200 = readRDS(paste0(work_path, "/tree_analysis/clade_coincodence/res_B2_K200.rds"))$log2_enr

res_B1_backbone_K15 = readRDS(paste0(work_path, "/tree_analysis/clade_coincodence/res_B1_K15_backbone.rds"))$log2_enr
res_B2_backbone_K15 = readRDS(paste0(work_path, "/tree_analysis/clade_coincodence/res_B2_K15_backbone.rds"))$log2_enr

res_B1_backbone_K200 = readRDS(paste0(work_path, "/tree_analysis/clade_coincodence/res_B1_K200_backbone.rds"))$log2_enr
res_B2_backbone_K200 = readRDS(paste0(work_path, "/tree_analysis/clade_coincodence/res_B2_K200_backbone.rds"))$log2_enr




cor.test(res_B1_K15[upper.tri(res_B1_K15)], res_B1_backbone_K15[upper.tri(res_B1_backbone_K15)], method = "spearman")
## 0.88

cor.test(res_B2_K15[upper.tri(res_B2_K15)], res_B2_backbone_K15[upper.tri(res_B2_backbone_K15)], method = "spearman")
## 0.83

cor.test(res_B1_K200[upper.tri(res_B1_K200)], res_B1_backbone_K200[upper.tri(res_B1_backbone_K200)], method = "spearman")
## 0.92

cor.test(res_B2_K200[upper.tri(res_B2_K200)], res_B2_backbone_K200[upper.tri(res_B2_backbone_K200)], method = "spearman")
## 0.90





