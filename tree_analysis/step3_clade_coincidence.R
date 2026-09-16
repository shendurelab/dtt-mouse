
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

tree <- read.tree(paste0(work_path, "/tree_analysis/merged_full_placed.nwk"))

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

# Subset cell types represented by at least 100 cells in each blastomere subtree
cell_meta_B1 = cell_meta[cell_meta$cell_id %in% tree_tips_B1,] %>%
    group_by(celltype) %>% tally() %>%
    filter(n >= 100)
cell_meta_B2 = cell_meta[cell_meta$cell_id %in% tree_tips_B2,] %>%
    group_by(celltype) %>% tally() %>%
    filter(n >= 100)

shared_celltypes = intersect(cell_meta_B1$celltype, cell_meta_B2$celltype)
# n = 88 cell types


for (K in c(15, 200)) {
  res <- run_clade_analysis(tree, cell_meta, tree_tips_B1, "Blastomere A", K, shared_celltypes)
  saveRDS(res, paste0(work_path, "/tree_analysis/clade_coincodence/res_B1_K", K, ".rds"))

  res <- run_clade_analysis(tree, cell_meta, tree_tips_B2, "Blastomere B", K, shared_celltypes)
  saveRDS(res, paste0(work_path, "/tree_analysis/clade_coincodence/res_B2_K", K, ".rds"))

  res <- run_clade_analysis(tree, cell_meta, c(tree_tips_B1, tree_tips_B2), "All", K, shared_celltypes)
  saveRDS(res, paste0(work_path, "/tree_analysis/clade_coincodence/res_all_K", K, ".rds"))
}




# ---- Report the number of clades and mean/median of internal nodes ----

for (K in c(15, 200)) {
  res <- readRDS(paste0(work_path, "/tree_analysis/clade_coincodence/res_B1_K", K, ".rds"))
  k_node = unique(res$clade_df[,c("node", "node_depth")])
  print(paste0("B1:K=", K, ", # = ", nrow(k_node), 
    ", ", round(mean(k_node$node_depth),2), " +/- ", round(sd(k_node$node_depth),2), "; median = ", round(median(k_node$node_depth),2) ))

  res <- readRDS(paste0(work_path, "/tree_analysis/clade_coincodence/res_B2_K", K, ".rds"))
  k_node = unique(res$clade_df[,c("node", "node_depth")])
  print(paste0("B2:K=", K, ", # = ", nrow(k_node), 
    ", ", round(mean(k_node$node_depth),2), " +/- ", round(sd(k_node$node_depth),2), "; median = ", round(median(k_node$node_depth),2) ))

  res <- readRDS(paste0(work_path, "/tree_analysis/clade_coincodence/res_all_K", K, ".rds"))
  k_node = unique(res$clade_df[,c("node", "node_depth")])
  print(paste0("All:K=", K, ", # = ", nrow(k_node), 
    ", ", round(mean(k_node$node_depth),2), " +/- ", round(sd(k_node$node_depth),2), "; median = ", round(median(k_node$node_depth),2) ))
}

[1] "B1:K=15, # = 82749, 11.16 +/- 1.58; median = 11.27"
[1] "B2:K=15, # = 57969, 11.08 +/- 1.67; median = 11.16"
[1] "All:K=15, # = 140718, 11.13 +/- 1.62; median = 11.23"
[1] "B1:K=200, # = 10047, 9.2 +/- 1.6; median = 8.85"
[1] "B2:K=200, # = 6769, 9.11 +/- 1.61; median = 8.79"
[1] "All:K=200, # = 16816, 9.16 +/- 1.6; median = 8.83"




########################
# ---- Plot heatmap ----

library(reshape2)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

work_path <- "/net/shendure/vol2/projects/cxqiu/work/tapemouse"

# ---- Load data ----
res_B1_K15  <- readRDS(paste0(work_path, "/tree_analysis/clade_coincodence/res_B1_K15.rds"))
res_B2_K15  <- readRDS(paste0(work_path, "/tree_analysis/clade_coincodence/res_B2_K15.rds"))
res_all_K15 <- readRDS(paste0(work_path, "/tree_analysis/clade_coincodence/res_all_K15.rds"))
res_B1_K200 <- readRDS(paste0(work_path, "/tree_analysis/clade_coincodence/res_B1_K200.rds"))
res_B2_K200 <- readRDS(paste0(work_path, "/tree_analysis/clade_coincodence/res_B2_K200.rds"))
res_all_K200<- readRDS(paste0(work_path, "/tree_analysis/clade_coincodence/res_all_K200.rds"))

library(ggplot2)
library(reshape2)
library(scales)
library(patchwork)

# ---- 1. Cell types + category info ----

x = res_all_K15$log2_enr
x2 <- x
diag(x2) <- 0
hc  <- hclust(dist(x2), method = "ward.D2")
ord <- rownames(x)[hc$order]
ord

write.table(ord, paste0(work_path, "/tree_analysis/clade_coincodence/cell_type_categories_old.tsv"), row.names=F, col.names=F, sep="\t", quote=F)

celltype_list <- read.table(paste0(work_path, "/tree_analysis/clade_coincodence/cell_type_categories_new.tsv"),
                            header = TRUE, sep = "\t")
ordered_ct <- celltype_list$cell_type
n_ct       <- length(ordered_ct)

# runs of categories along ordered_ct
cat_runs      <- rle(celltype_list$category)
cat_ends      <- cumsum(cat_runs$lengths)          # last position of each block
cat_starts    <- c(1, head(cat_ends, -1) + 1)
cat_mids      <- (cat_starts + cat_ends) / 2       # centers, for x-axis labels
cat_names     <- cat_runs$values
inner_sep_pos <- head(cat_ends, -1) + 0.5          # between-block positions

# ---- 2. Shared limits across K=15 and K=200 ----
max_offdiag <- function(mat) { m <- mat; diag(m) <- NA; max(m, na.rm = TRUE) }

lims_shared <- c(
  min(c(res_all_K15$log2_enr, res_all_K200$log2_enr), na.rm = TRUE),
  max(c(max_offdiag(res_all_K15$log2_enr),
        max_offdiag(res_all_K200$log2_enr)))
)

lims_shared[1] = -5
lims_shared[2] = 3.5

# ---- 3. Diverging palette: cyan -> dark blue -> black -> orange -> yellow ----
make_palette <- function(lims, n_total = 200) {
  zero_pos <- (0 - lims[1]) / (lims[2] - lims[1])
  zero_pos <- max(0.05, min(0.95, zero_pos))
  n_neg <- round(n_total * zero_pos)
  n_pos <- n_total - n_neg
  neg_ramp <- colorRampPalette(c("#5FCDF0", "#1E5E8E", "#0A2540", "#000000"))(n_neg)
  pos_ramp <- colorRampPalette(c("#000000", "#3a0a00", "#8a1c00", "#d15200",
                                 "#ff8f1f", "#ffc140", "#ffec6b"))(n_pos)
  c(neg_ramp, pos_ramp)
}
Colors_shared <- make_palette(lims_shared)

# ---- 4. Plot function ----
plot_enrichment <- function(res, ordered_ct, colors, limits) {
  mat <- res$log2_enr[ordered_ct, ordered_ct]  # enforce order
  df  <- melt(mat, varnames = c("CT_i", "CT_j"), value.name = "log2_enr")
  df$CT_i <- factor(df$CT_i, levels = rev(ordered_ct))   # first cell type on top
  df$CT_j <- factor(df$CT_j, levels = ordered_ct)

  # horizontal separator y-positions must be flipped because y-levels are reversed
  h_sep <- n_ct - inner_sep_pos + 1

  ggplot(df, aes(x = CT_j, y = CT_i, fill = log2_enr)) +
    geom_tile() +
    geom_vline(xintercept = inner_sep_pos, color = "white", linewidth = 0.4) +
    geom_hline(yintercept = h_sep,         color = "white", linewidth = 0.4) +
    scale_fill_gradientn(colors = colors,
                         name   = expression(log[2]~"enrichment of clade co-occurrence"),
                         limits = limits, oob = scales::squish) +
    scale_x_discrete(breaks = ordered_ct[round(cat_mids)],
                     labels = cat_names,
                     expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_fixed() +
    labs(x = NULL, y = NULL,
         title = paste0(res$blastomere, "   K = ", res$K)) +
    theme_minimal(base_size = 7) +
    theme(axis.text.x  = element_text(angle = 30, hjust = 1, size = 8),
          axis.text.y  = element_text(size = 5),
          panel.grid   = element_blank(),
          plot.title   = element_text(hjust = 0.5, size = 10))
}

# ---- 5. Plots + combine (shared colorbar) ----
p_K15  <- plot_enrichment(res_all_K15,  ordered_ct, Colors_shared, lims_shared)
p_K200 <- plot_enrichment(res_all_K200, ordered_ct, Colors_shared, lims_shared)

combined <- (p_K15 / p_K200) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave("~/share/clade_enrichment_combined_all.pdf", combined, width = 7, height = 15)


p_K15_B1  <- plot_enrichment(res_B1_K15,  ordered_ct, Colors_shared, lims_shared)
p_K200_B1 <- plot_enrichment(res_B1_K200, ordered_ct, Colors_shared, lims_shared)

p_K15_B2  <- plot_enrichment(res_B2_K15,  ordered_ct, Colors_shared, lims_shared)
p_K200_B2 <- plot_enrichment(res_B2_K200, ordered_ct, Colors_shared, lims_shared)

combined <- (p_K15_B1 + p_K200_B1) / (p_K15_B2 + p_K200_B2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave("~/share/clade_enrichment_combined_blastomere.pdf", combined, width = 15, height = 15)


### off-diagonal correlation

x = res_B1_K15$log2_enr
y = res_B2_K15$log2_enr

x = x[upper.tri(x)]
y = y[upper.tri(y)]

cor.test(x, y, method = "spearman")
### 0.7927594
### < 2.2e-16












