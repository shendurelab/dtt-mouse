

#################################################################################################
### Infer the transcriptional states for 20% of hold-out internal nodes by its cousin nodes


#################################################################
### Step-1: IMPUTE the PCA coordinates of hold-out internal nodes

source("~/work/scripts/utils.R")
library(dplyr)
library(tidyr)
library(ape)
library(ggplot2)
library(ggrepel)
library(phangorn)
library(castor)

work_path <- "/net/shendure/vol2/projects/cxqiu/work/tapemouse"

tree <- read.tree(paste0(work_path, "/tree_analysis/merged_full_placed.nwk"))
internal <- read.table(paste0(work_path, "/tree_analysis/impute_nodes/internal.tsv"), header = TRUE)
internal$ape_id <- as.integer(sub("node_", "", internal$node))
pca_coor <- read.table(paste0(work_path, "/tree_analysis/impute_nodes/pca_impute_full.tsv"))
pca_coor <- pca_coor[internal$node, ]

pd_impute <- read.table(paste0(work_path, "/tree_analysis/impute_nodes/pd_nodes_infer_200.txt"), header = TRUE, sep="\t")
node_impute <- pd_impute %>% filter(celltype != "missing", cell_id == "internal_node") %>% pull(node_id)

day_list = c("E13.25","E13.0","E12.75","E12.5","E12.25","E12.0","E11.75","E11.5","E11.25","E11.0","E10.75","E10.5","E10.25","E10.0","E9.75","E9.5","E9.25","E9.0","E8.75","E8.5")


# ---- Setup: parent + children lookups (needed to get ancestors/descendants) ----

n_tips <- Ntip(tree)
n_int  <- Nnode(tree)
n_total <- n_tips + n_int

parent_of <- integer(n_total)
parent_of[tree$edge[, 2]] <- tree$edge[, 1]

children_of <- vector("list", n_total)
for (i in seq_len(nrow(tree$edge))) {
  pa <- tree$edge[i, 1]; ch <- tree$edge[i, 2]
  children_of[[pa]] <- c(children_of[[pa]], ch)
}

get_ancestors <- function(node) {
  ancestors <- integer(0)
  p <- parent_of[node]
  while (p != 0) { ancestors <- c(ancestors, p); p <- parent_of[p] }
  ancestors
}

get_descendants <- function(node) {
  descendants <- integer(0)
  stack <- node
  while (length(stack) > 0) {
    current <- stack[1]; stack <- stack[-1]
    ch <- children_of[[current]]
    if (length(ch) > 0) { descendants <- c(descendants, ch); stack <- c(stack, ch) }
  }
  descendants
}





# ---- Hold-out setup ----

args = commandArgs(trailingOnly=TRUE)
kk = as.numeric(args[1])
seed = kk-1

set.seed(seed)
hold_out_ape_ids <- internal %>% filter(day %in% day_list, node %in% node_impute) %>%
    group_by(day) %>% slice_sample(n = 1000) %>% pull(ape_id)

# ---- Group internal nodes by day (for fast lookup) ----

nodes_by_day <- split(internal$ape_id, internal$day)

# For each query, we'll need to know its day
day_of <- setNames(internal$day, internal$ape_id)

# ---- Find K nearest cousins for one query ----

K <- 10

nearest_cousins <- function(q) {
  # Candidates = same-day internal nodes, minus q + ancestors + descendants
  candidates <- nodes_by_day[[ day_of[as.character(q)] ]]
  excluded   <- c(q, get_ancestors(q), get_descendants(q))
  candidates <- setdiff(candidates, excluded)
  
  d <- castor::get_pairwise_distances(tree,
                                      A = rep(q, length(candidates)),
                                      B = candidates)
  candidates[order(d)[1:K]]
}

# ---- Test on 10 queries ----

t0 <- Sys.time()
test_result <- t(sapply(hold_out_ape_ids, nearest_cousins))
t1 <- Sys.time()
print(t1 - t0)

rownames(test_result) = hold_out_ape_ids
colnames(test_result) = paste0("NN_", 1:ncol(test_result))

write.table(test_result, paste0(work_path, "/tree_analysis/impute_nodes/impute_by_cousin/node_NN_seed_", seed, ".txt"), row.names=T, col.names=T, sep="\t", quote=F)



#################################################################
### Step-2: IMPUTE the PCA coordinates of hold-out internal nodes

import numpy as np
import pandas as pd
from collections import defaultdict, deque

work_path = "/net/shendure/vol2/projects/cxqiu/work/tapemouse"

pca_coor = pd.read_csv(f"{work_path}/tree_analysis/impute_nodes/pca_impute_full.tsv", sep="\t", index_col=0)

for seed in range(10):
    print(seed)
    NN = pd.read_csv(f"{work_path}/tree_analysis/impute_nodes/impute_by_cousin/node_NN_seed_{seed}.txt", sep="\t", index_col=0)
    ref = pca_coor.index
    center = "node_" + NN.index.astype(str)
    missing_center = center[~center.isin(ref)]
    print(f"center nodes: {len(missing_center)} / {len(center)} missing")
    for col in NN.columns:
        labels = "node_" + NN[col].astype(str)
        missing = labels[~labels.isin(ref)]
        print(f"{col}: {len(missing)} / {len(labels)} missing")
    M = pca_coor.to_numpy()
    center = "node_" + NN.index.astype(str)
    row_pos = pca_coor.index.get_indexer(center)
    assert (row_pos != -1).all(), "some NN center nodes not in pca_coor"
    nbr_pos = pca_coor.index.get_indexer(("node_" + NN.astype(str)).to_numpy().ravel())
    assert (nbr_pos != -1).all()
    mean = M[nbr_pos].reshape(*NN.shape, -1).mean(axis=1)    # (20000, 50)
    pca_coor_new = pca_coor.copy()
    pca_coor_new.iloc[row_pos, :] = mean
    pca_coor_new.to_csv(f"{work_path}/tree_analysis/impute_nodes/impute_by_cousin/pca_impute_seed_{seed}.tsv", sep="\t")




###################################################################
#### Step-3: Identify mutual nearest neighbors to assign cell types

import numpy as np
import pandas as pd
import sys

work_path = "/net/shendure/vol2/projects/cxqiu/work/tapemouse"

pca_coor = pd.read_csv(f"{work_path}/transcriptome_analysis/adata_integration_early.pca.csv", index_col = 0)
pd_meta = pd.read_csv(f"{work_path}/tree_analysis/impute_nodes/pd.txt", sep = "\t", index_col = 0)

pca_coor.index   = pd_meta.index
pca_coor.columns = [f"PC_{i}" for i in range(1, pca_coor.shape[1] + 1)]

seed = int(sys.argv[2] - 1)

internal = pd.read_csv(f"{work_path}/tree_analysis/impute_nodes/internal.tsv", sep="\t")
internal.index = internal["node"]
pca_impute = pd.read_csv(f"{work_path}/tree_analysis/impute_nodes/impute_by_cousin/pca_impute_seed_{seed}.tsv", sep="\t", index_col=0)

day_list = ["E13.25","E13.0","E12.75","E12.5","E12.25","E12.0","E11.75","E11.5","E11.25","E11.0","E10.75","E10.5","E10.25","E10.0","E9.75","E9.5","E9.25","E9.0","E8.75","E8.5"]

day = day_list[int(sys.argv[1]) - 1]
print(day)

pd_node  = internal[internal["day"] == day]
pca_node = pca_impute.loc[pd_node["node"].values]

keep = ~pca_node.isna().any(axis=1)
pca_node = pca_node.loc[keep]
pd_node  = pd_node.loc[keep.values]

pd_cell  = pd_meta[pd_meta["day"] == day]
pca_cell = pca_coor.loc[pd_cell["cell_id"].values]

from annoy import AnnoyIndex
from scipy.sparse import csr_matrix

k_list = [200]

n_dim = pca_node.shape[1]
n_trees = 150

# ---- Build Annoy index on internal nodes; query with cells ----
ann_node = AnnoyIndex(n_dim, "euclidean")
for i, vec in enumerate(pca_node.values):
    ann_node.add_item(i, vec)

ann_node.build(n_trees)

# ---- Build Annoy index on cells; query with internal nodes ----
ann_cell = AnnoyIndex(n_dim, "euclidean")
for i, vec in enumerate(pca_cell.values):
    ann_cell.add_item(i, vec)

ann_cell.build(n_trees)

for k in k_list:
    
    idx_cell_to_node = np.array([
    ann_node.get_nns_by_vector(vec, k) for vec in pca_cell.values
    ])   # (n_cells, k)

    idx_node_to_cell = np.array([
        ann_cell.get_nns_by_vector(vec, k) for vec in pca_node.values
    ])   # (n_nodes, k)

    # ---- Build boolean membership matrices ----
    n_cells = pca_cell.shape[0]
    n_nodes = pca_node.shape[0]

    rows1 = np.repeat(np.arange(n_cells), k)
    cols1 = idx_cell_to_node.ravel()
    m1 = csr_matrix((np.ones_like(rows1, dtype=np.int8), (rows1, cols1)),
                    shape=(n_cells, n_nodes))

    rows2 = idx_node_to_cell.ravel()
    cols2 = np.repeat(np.arange(n_nodes), k)
    m2 = csr_matrix((np.ones_like(rows2, dtype=np.int8), (rows2, cols2)),
                    shape=(n_cells, n_nodes))

    # ---- Mutual nearest neighbors ----
    mutual = m1.multiply(m2).tocoo()

    cell_ids = pca_cell.index.to_numpy()
    node_ids = pca_node.index.to_numpy()

    hits_df = pd.DataFrame({
        "cell_id": cell_ids[mutual.row],
        "node_id": node_ids[mutual.col],
        "day":     day
    })

    print("unique cells:", hits_df["cell_id"].nunique())
    print("unique nodes:", hits_df["node_id"].nunique())

    hits_df.to_csv(f"{work_path}/tree_analysis/impute_nodes/impute_by_cousin/assign/hits_df_{day}_{k}_seed_{seed}.txt", sep="\t", index=False)


### HOW to submit?

for seed in $(seq 0 9); do
cat > run_NN_"$seed".sh <<EOF
#!/bin/bash
num=\${SGE_TASK_ID}
/net/gs/vol1/home/cxqiu/miniconda/miniconda/bin/python /net/gs/vol1/home/cxqiu/bin/run_NN.py \${num} $seed
EOF
chmod +x run_NN_"$seed".sh
done

for seed in $(seq 0 9); do
qsub -t 1-20 -l mfree=50G,hostname="$gpu_node",gpgpu=1 run_NN_"$seed".sh
done


############################################################################
### Step-4: Assign cel type labels for internal nodes based on its neighbors


source("~/work/scripts/utils.R")
library(dplyr)
library(tidyr)
library(ape)
library(ggplot2)
library(ggrepel)
library(phangorn)

work_path <- "/net/shendure/vol2/projects/cxqiu/work/tapemouse"

internal = read.table(paste0(work_path, "/tree_analysis/impute_nodes/internal.tsv"), header=T, sep="\t")

pd = read.table(paste0(work_path, "/tree_analysis/impute_nodes/pd.txt"), header=T, row.names=1, sep="\t")

pd_E13.5 = read.table(paste0(work_path, "/tree_analysis/impute_nodes/pd_E13.5.tsv"), header=T, row.names=1, sep="\t")

day_list = c("E13.25","E13.0","E12.75","E12.5","E12.25","E12.0","E11.75","E11.5","E11.25","E11.0","E10.75","E10.5","E10.25","E10.0","E9.75","E9.5","E9.25","E9.0","E8.75","E8.5")

args = commandArgs(trailingOnly=TRUE)
seed = as.numeric(args[1])

pd_internal_list = list()
for(day in day_list){
    print(day)
    pd_internal_list[[day]] = read.table(paste0(work_path, "/tree_analysis/impute_nodes/impute_by_cousin/assign/hits_df_", day, "_200_seed_", seed, ".txt"), header=T, sep='\t')
}
pd_internal = do.call(rbind, pd_internal_list)

pd_internal_x = pd_internal %>% left_join(pd[,c("cell_id", "celltype")], by = "cell_id")

pd_internal_y = pd_internal_x %>% group_by(node_id, celltype) %>% tally() %>% ungroup() %>%
    group_by(node_id) %>% slice_max(order_by = n, n = 1, with_ties = FALSE)
internal_x = internal %>% select(node_id = node, day) %>% left_join(pd_internal_y[,c("node_id", "celltype")], by = "node_id")
internal_x$celltype[is.na(internal_x$celltype)] = "missing"
write.table(internal_x, paste0(work_path, "/tree_analysis/impute_nodes/impute_by_cousin/pd_nodes_infer_200_seed_", seed, ".txt"), row.names=F, sep="\t", quote=F)


###################################
### Step-5: Summarizing the results

source("~/work/scripts/utils.R")
work_path <- "/net/shendure/vol2/projects/cxqiu/work/tapemouse"

day_list = c("E13.25","E13.0","E12.75","E12.5","E12.25","E12.0","E11.75","E11.5","E11.25","E11.0","E10.75","E10.5","E10.25","E10.0","E9.75","E9.5","E9.25","E9.0","E8.75","E8.5")

pd_orig = read.table(paste0(work_path, "/tree_analysis/impute_nodes/pd_nodes_infer_200.txt"), sep="\t", header=T)
pd_orig = pd_orig[pd_orig$cell_id == "internal_node",]
pd_orig = pd_orig[,c("node_id", "celltype")]
colnames(pd_orig) = c("node_id", "celltype_orig")

df_list = list()
for(seed in c(0:9)){
    print(seed)
    node_NN = read.table(paste0(work_path, "/tree_analysis/impute_nodes/impute_by_cousin/node_NN_seed_", seed, ".txt"), sep="\t", header=T)
    pd_mask = read.table(paste0(work_path, "/tree_analysis/impute_nodes/impute_by_cousin/pd_nodes_infer_200_seed_", seed, ".txt"), sep="\t", header=T)
    pd_mask = pd_mask[pd_mask$node_id %in% paste0("node_", rownames(node_NN)),]
    pd_mask = pd_mask %>% left_join(pd_orig, by = "node_id")
    pd_mask$seed = seed
    df_list[[seed + 1]] = pd_mask
}
df = do.call(rbind, df_list)

saveRDS(df, paste0(work_path, "/tree_analysis/impute_nodes/impute_by_cousin/internal_nodes_assigned_accuracy.rds"))


sum(df$celltype == df$celltype_orig)
49103/200000 = 25%


### stratified by embryonic day

df_plot <- df %>%
  group_by(day, seed) %>%
  summarise(
    total_n = n(),
    kept_n  = sum(celltype_orig == celltype),
    pct     = 100 * kept_n / total_n,
    .groups = "drop"
  ) %>%
  mutate(day = factor(day, levels = rev(day_list)))

df_plot %>% group_by(day) %>% summarize(median_pct = median(pct))

   day    median_pct
   <fct>       <dbl>
 1 E8.5         17.0
 2 E8.75        16.8
 3 E9.0         17.4
 4 E9.25        20.6
 5 E9.5         20.8
 6 E9.75        22.0
 7 E10.0        21.7
 8 E10.25       22.9
 9 E10.5        22.4
10 E10.75       25.2
11 E11.0        24.8
12 E11.25       25.4
13 E11.5        27.4
14 E11.75       26.6
15 E12.0        25.6
16 E12.25       28.2
17 E12.5        25.2
18 E12.75       33
19 E13.0        38.6
20 E13.25       28.8

df_plot$day = factor(df_plot$day, levels = rev(day_list))

p = ggplot(df_plot, aes(day, pct, fill = day)) + 
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.5) +
    theme_classic(base_size = 10) +
    scale_fill_manual(values=day_color_plate) +
    theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("~/share/boxplot_assign_accuracy_by_cousin.pdf", p, width = 6, height = 4)





