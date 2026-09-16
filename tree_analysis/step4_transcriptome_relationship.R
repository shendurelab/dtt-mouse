
##################################################
### Relationship between transcriptome and lineags

source("~/work/scripts/utils.R")
library(dplyr)
library(tidyr)
library(ape)
library(ggplot2)
library(ggrepel)

work_path <- "/net/shendure/vol2/projects/cxqiu/work/tapemouse"

# ---- Load data ----
cell_meta <- read.table(paste0(work_path, "/tree_analysis/cell_metadata.v8.txt"),
                        header = TRUE, sep = "\t")

tree <- read.tree(paste0(work_path, "/tree_analysis/merged_full_placed.nwk"))

# ---- Function 1: extract sibling pairs and annotate ----
analyze_clade <- function(subtree, cell_meta) {

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

  sibling_celltypes <- sibling_pairs %>%
    left_join(cell_meta %>% select(cell_id, celltype), by = c("cell_A" = "cell_id")) %>%
    rename(celltype_A = celltype) %>%
    left_join(cell_meta %>% select(cell_id, celltype), by = c("cell_B" = "cell_id")) %>%
    rename(celltype_B = celltype)
  
  return(sibling_celltypes)
}

sibling = analyze_clade(tree, cell_meta)

######################################
### Step-1: create table (randomly selected 1000 pairs of sibling from the same cell type and matched non-sibling pairs)

set.seed(1)
MAX_PAIRS <- 10000

sib_same <- sibling %>%
    filter(celltype_A == celltype_B) %>%
    transmute(cell_A, cell_B, celltype = celltype_A)

sib_sampled <- sib_same %>%
    group_by(celltype) %>%
    slice_sample(n = MAX_PAIRS) %>%
    ungroup() %>%
    mutate(sibling = TRUE)

sib_key <- sib_same %>%
    transmute(key = paste(pmin(cell_A, cell_B), pmax(cell_A, cell_B), sep = "|")) %>%
    pull(key) %>%
    unique()

ct_to_cells <- split(cell_meta$cell_id, cell_meta$celltype)

nonsib_sampled <- sib_sampled %>%
    count(celltype, name = "n_pairs") %>%
    rowwise() %>%
    do({
        ct <- .$celltype
        n  <- .$n_pairs
        pool <- ct_to_cells[[ct]]
        if (length(pool) < 2) return(tibble())
        need <- n
        picked <- tibble()
        tries <- 0
        while (nrow(picked) < n && tries < 20) {
            m  <- max(need * 2, 50)
            a  <- sample(pool, m, replace = TRUE)
            b  <- sample(pool, m, replace = TRUE)
            df <- tibble(cell_A = a, cell_B = b) %>%
                filter(cell_A != cell_B) %>%
                mutate(key = paste(pmin(cell_A, cell_B), pmax(cell_A, cell_B), sep = "|")) %>%
                filter(!key %in% sib_key) %>%
                distinct(key, .keep_all = TRUE) %>%
                select(-key)
            picked <- bind_rows(picked, df) %>% distinct()
            need   <- n - nrow(picked)
            tries  <- tries + 1
        }
        picked %>%
            slice_head(n = n) %>%
            mutate(celltype = ct)
    }) %>%
    ungroup() %>%
    mutate(sibling = FALSE)

pairs_tbl <- bind_rows(sib_sampled, nonsib_sampled) %>%
    select(cell_A, cell_B, celltype, sibling)

pairs_tbl %>% count(celltype, sibling)

write.table(pairs_tbl, paste0(work_path, "/tree_analysis/transcriptome_relationship/sibling_pairs.txt"), row.names=F, col.names=T, sep="\t", quote=F)


################
### calculate distance in PCA space

import numpy as np
import pandas as pd

work_path = "/net/shendure/vol2/projects/cxqiu/work/tapemouse"

pca_coor = pd.read_csv(f"{work_path}/transcriptome_analysis/adata_integration.pca.csv", index_col = 0)
pd_meta = pd.read_csv(f"{work_path}/transcriptome_analysis/adata_integration.obs.csv", index_col = 0)

sib = pd.read_csv(f"{work_path}/tree_analysis/transcriptome_relationship/sibling_pairs.txt", sep="\t")

# --- Map cell_id -> row index in pca_coor (both share row order with pd_meta) ---
idx = pd.Series(np.arange(len(pd_meta)), index=pd_meta.index)

# --- Look up row indices for each pair; drop pairs with missing cells ---
iA = idx.reindex(sib['cell_A']).values
iB = idx.reindex(sib['cell_B']).values
ok = ~(np.isnan(iA) | np.isnan(iB))

sib = sib.loc[ok].reset_index(drop=True)
iA  = iA[ok].astype(np.int64)
iB  = iB[ok].astype(np.int64)

# --- Ensure pca_coor is a numpy array (fast row indexing) ---
X = pca_coor.values if isinstance(pca_coor, pd.DataFrame) else pca_coor

# --- Vectorized Euclidean distance ---
diff = X[iA] - X[iB]
sib['euclid_pc50'] = np.sqrt(np.einsum('ij,ij->i', diff, diff))

print(sib.head())
sib.to_csv(f"{work_path}/tree_analysis/transcriptome_relationship/sibling_pairs_dist.txt",
           sep="\t", index=False)



################
### making plot

res = read.table(paste0(work_path, "/tree_analysis/transcriptome_relationship/sibling_pairs_dist.txt"), sep="\t", header=T)
res <- res %>%
    mutate(sibling = if_else(sibling == "True", "sibling", "not_sibling"),
           sibling = factor(sibling, levels = c("sibling", "not_sibling")))

p = ggplot(res, aes(x = sibling, y = euclid_pc50, fill = sibling)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
    scale_fill_manual(values = c("sibling" = "#d1495b", "not_sibling" = "#8d99ae")) +
    coord_cartesian(ylim = c(0, 25)) +           # zoom, doesn't drop data
    labs(x = NULL, y = "Euclidean distance") +
    theme_classic(base_size = 11) +
    theme(legend.position = "none")

fit = wilcox.test(res$euclid_pc50[res$sibling == "sibling"], res$euclid_pc50[res$sibling != "sibling"])
### p < 2.2e-16

ggsave("~/share/sibling_pairs_dist.pdf", p, height=5, width=5)






##################################################
### Relationship between transcriptome and lineags

source("~/work/scripts/utils.R")
library(dplyr)
library(tidyr)
library(ape)
library(ggplot2)
library(ggrepel)
library(castor)

work_path <- "/net/shendure/vol2/projects/cxqiu/work/tapemouse"

# ---- Load data ----
cell_meta <- read.table(paste0(work_path, "/tree_analysis/cell_metadata.v8.txt"),
                        header = TRUE, sep = "\t")
tree <- read.tree(paste0(work_path, "/tree_analysis/merged_full_placed.nwk"))
set.seed(1)
PAIRS_PER_CT_CAP <- 100000
PAIRS_PER_BIN    <- 10000

# --- 1. Restrict cell_meta to cells on the tree ---
tips_on_tree   <- tree$tip.label
cell_meta_tree <- cell_meta %>% filter(cell_id %in% tips_on_tree)
tip_idx        <- setNames(seq_along(tree$tip.label), tree$tip.label)
ct_to_cells    <- split(cell_meta_tree$cell_id, cell_meta_tree$celltype)

# --- Precompute distance from every tip to the root ---
root_dist_all  <- castor::get_all_distances_to_root(tree)   # length: Ntips + Nnodes
tip_root_dist  <- root_dist_all[seq_along(tree$tip.label)]  # tip depths only

# --- 2. Sample same-celltype pairs (candidate pool) -----------------------
candidates <- lapply(names(ct_to_cells), function(ct) {
    pool <- ct_to_cells[[ct]]
    if (length(pool) < 2) return(NULL)
    n_possible <- choose(length(pool), 2)
    n_draw     <- min(n_possible, PAIRS_PER_CT_CAP)
    m <- ceiling(n_draw * 1.2)
    a <- sample(pool, m, replace = TRUE)
    b <- sample(pool, m, replace = TRUE)
    tibble(cell_A = a, cell_B = b, celltype = ct) %>%
        filter(cell_A != cell_B) %>%
        distinct(cell_A, cell_B, .keep_all = TRUE) %>%
        slice_head(n = n_draw)
}) %>% bind_rows()

# --- 3. Compute pairwise tree distance and derive MRCA depth --------------
iA <- tip_idx[candidates$cell_A]
iB <- tip_idx[candidates$cell_B]

d_AB     <- get_pairwise_distances(tree, iA, iB)
depth_A  <- tip_root_dist[iA]
depth_B  <- tip_root_dist[iB]

candidates$lineage_dist <- d_AB
candidates$mrca_depth   <- (depth_A + depth_B - d_AB) / 2   # depth of MRCA from root

# If your tree is rooted at E0 (or a known root age), mrca_depth = E-day of MRCA
# Otherwise treat it as an ordinal depth measure

# --- 4. Bin MRCA depth at width 1 -----------------------------------------
candidates <- candidates %>%
    mutate(dist_bin = floor(mrca_depth))

# --- 5. Down-sample per bin -----------------------------------------------
sampled <- candidates %>%
    group_by(dist_bin) %>%
    slice_sample(n = PAIRS_PER_BIN) %>%
    ungroup() %>%
    select(cell_A, cell_B, lineage_dist, mrca_depth, celltype, dist_bin)

# --- 6. Inspect -----------------------------------------------------------
sampled %>% count(dist_bin) %>% arrange(dist_bin)
head(sampled)

write.table(sampled, paste0(work_path, "/tree_analysis/transcriptome_relationship/distance_pairs.txt"), row.names=F, col.names=T, sep="\t", quote=F)



################
### calculate distance in PCA space

import numpy as np
import pandas as pd

work_path = "/net/shendure/vol2/projects/cxqiu/work/tapemouse"

pca_coor = pd.read_csv(f"{work_path}/transcriptome_analysis/adata_integration.pca.csv", index_col = 0)
pd_meta = pd.read_csv(f"{work_path}/transcriptome_analysis/adata_integration.obs.csv", index_col = 0)

sib = pd.read_csv(f"{work_path}/tree_analysis/transcriptome_relationship/distance_pairs.txt", sep="\t")

# --- Map cell_id -> row index in pca_coor (both share row order with pd_meta) ---
idx = pd.Series(np.arange(len(pd_meta)), index=pd_meta.index)

# --- Look up row indices for each pair; drop pairs with missing cells ---
iA = idx.reindex(sib['cell_A']).values
iB = idx.reindex(sib['cell_B']).values
ok = ~(np.isnan(iA) | np.isnan(iB))

sib = sib.loc[ok].reset_index(drop=True)
iA  = iA[ok].astype(np.int64)
iB  = iB[ok].astype(np.int64)

# --- Ensure pca_coor is a numpy array (fast row indexing) ---
X = pca_coor.values if isinstance(pca_coor, pd.DataFrame) else pca_coor

# --- Vectorized Euclidean distance ---
diff = X[iA] - X[iB]
sib['euclid_pc50'] = np.sqrt(np.einsum('ij,ij->i', diff, diff))

print(sib.head())
sib.to_csv(f"{work_path}/tree_analysis/transcriptome_relationship/distance_pairs_dist.txt",
           sep="\t", index=False)



################
### making plot

res = read.table(paste0(work_path, "/tree_analysis/transcriptome_relationship/distance_pairs_dist.txt"), sep="\t", header=T)
res$MRCA_day = paste0("E", res$dist_bin)
res$MRCA_day = factor(res$MRCA_day, levels = paste0("E", 0:13))

# --- Plot ---
p <- ggplot(res, aes(x = MRCA_day, y = euclid_pc50, fill = MRCA_day)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
    scale_fill_viridis_d(option = "D") +
    coord_flip(ylim = c(0, 25)) +
    labs(x = "MRCA", y = "Euclidean distance") +
    theme_classic(base_size = 11) +
    theme(legend.position = "none")

ggsave("~/share/distance_pairs_dist.pdf", p, height=5, width=4)

cor.test(res$dist_bin, res$euclid_pc50, method = "spearman")
### r = -0.1, p = < 1e-250






