

#################################################################################################
### Infer the transcriptional states for individiul internal nodes (randomly masking 20% of tips)


########################################################
### Step-2: IMPUTE the PCA coordinates of internal nodes

import numpy as np
import pandas as pd
from collections import defaultdict, deque

work_path = "/net/shendure/vol2/projects/cxqiu/work/tapemouse"

pca_coor = pd.read_csv(f"{work_path}/transcriptome_analysis/adata_integration_early.pca.csv", index_col = 0)
pd_meta = pd.read_csv(f"{work_path}/tree_analysis/impute_nodes/pd.txt", sep = "\t", index_col = 0)
pca_coor.index   = pd_meta.index
pca_coor.columns = [f"PC_{i}" for i in range(1, pca_coor.shape[1] + 1)]

pd_E135 = pd.read_csv(f"{work_path}/tree_analysis/impute_nodes/pd_E13.5.tsv", sep = "\t", index_col = 0)
pca_impute = pca_coor.loc[pd_E135['cell_id']]
pca_impute.index = pd_E135['node_id'].values
### pca_impute is pca coor of E13.5 cells and row ids are node_id

internal = pd.read_csv(f"{work_path}/tree_analysis/impute_nodes/internal.tsv", sep="\t")
edge             = pd.read_csv(f"{work_path}/tree_analysis/impute_nodes/edge.tsv", sep="\t")
children_map = edge.groupby("parent")["child"].apply(list).to_dict()

n_pc         = pca_impute.shape[1]
tip_ids      = pca_impute.index.tolist()
internal_ids = internal["node"].tolist()



# ---- Topological sort: guarantee children-before-parent ----
internal_set = set(internal_ids)
parents_of   = defaultdict(list)   # child -> list of parents
child_count  = defaultdict(int)    # parent -> number of children still to process

for parent, child in zip(edge["parent"], edge["child"]):
    if parent in internal_set:
        child_count[parent] += 1
        parents_of[child].append(parent)

# Seed queue with tips and any internal with no children in internal_set
queue = deque(tip_ids)
topo_order = []
for n in internal_ids:
    if child_count[n] == 0:
        topo_order.append(n)
        queue.append(n)

while queue:
    node = queue.popleft()
    for p in parents_of.get(node, []):
        child_count[p] -= 1
        if child_count[p] == 0:
            topo_order.append(p)
            queue.append(p)

assert len(topo_order) == len(internal_ids), \
    f"topo sort produced {len(topo_order)} nodes, expected {len(internal_ids)}"

internal_ids = topo_order
# ---- Now internal_ids is ordered by topology (from bottom to top, children is always visited before its parent)



# ---- Randomly permuting the tips ----
for seed in range(10):
    all_ids = tip_ids + internal_ids
    idx_of  = {nid: i for i, nid in enumerate(all_ids)}
    all_pcs = np.full((len(all_ids), n_pc), np.nan, dtype=np.float32)
    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(tip_ids))
    all_pcs[:len(tip_ids)] = pca_impute.to_numpy()[perm]
    for cnt, nd in enumerate(internal_ids, 1):
        if cnt % 10000 == 0:
            print(f"Seed={seed}: {cnt}/{len(internal_ids)}")
        children = children_map.get(nd)
        child_idx = [idx_of[c] for c in children]
        all_pcs[idx_of[nd]] = all_pcs[child_idx].mean(axis=0)
    pca_impute_full = pd.DataFrame(all_pcs, index=all_ids, columns=pca_impute.columns)
    pca_impute_full.to_csv(f"{work_path}/tree_analysis/impute_nodes/permute_tips/pca_impute_seed_{seed}.tsv", sep="\t")



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

seed = int(sys.argv[2])

internal = pd.read_csv(f"{work_path}/tree_analysis/impute_nodes/internal.tsv", sep="\t")
internal.index = internal["node"]
pca_impute = pd.read_csv(f"{work_path}/tree_analysis/impute_nodes/permute_tips/pca_impute_seed_{seed}.tsv", sep="\t", index_col=0)

day_list = ["E13.25","E13.0","E12.75","E12.5","E12.25","E12.0","E11.75","E11.5","E11.25","E11.0","E10.75","E10.5","E10.25","E10.0","E9.75","E9.5","E9.25","E9.0","E8.75","E8.5"]

day = day_list[int(sys.argv[1]) - 1]
print(day)

pd_node  = internal[internal["day"] == day]
pca_node = pca_impute.loc[pd_node["node"].values]

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

    hits_df.to_csv(f"{work_path}/tree_analysis/impute_nodes/permute_tips/assign/hits_df_{day}_{k}_seed_{seed}.txt", sep="\t", index=False)





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
    pd_internal_list[[day]] = read.table(paste0(work_path, "/tree_analysis/impute_nodes/permute_tips/assign/hits_df_", day, "_200_seed_", seed, ".txt"), header=T, sep='\t')
}
pd_internal = do.call(rbind, pd_internal_list)

pd_internal_x = pd_internal %>% left_join(pd[,c("cell_id", "celltype")], by = "cell_id")
write.table(pd_internal_x, paste0(work_path, "/tree_analysis/impute_nodes/permute_tips/internal_assign_cells_200_seed_", seed, ".txt"), row.names=F, col.names=F, sep="\t", quote=F)

node_cell_num = pd_internal_x %>% group_by(node_id) %>% tally()
print(paste0('Seed = ', seed, ', ', 100*round(nrow(node_cell_num)/nrow(internal),2), "% of nodes assigned MNNs, ", round(mean(node_cell_num$n), 2), " cells per node"))

pd_internal_y = pd_internal_x %>% group_by(node_id, celltype) %>% tally() %>% ungroup() %>%
    group_by(node_id) %>% slice_max(order_by = n, n = 1, with_ties = FALSE)
internal_x = internal %>% select(node_id = node, day) %>% left_join(pd_internal_y[,c("node_id", "celltype")], by = "node_id")
internal_x$celltype[is.na(internal_x$celltype)] = "missing"
write.table(internal_x, paste0(work_path, "/tree_analysis/impute_nodes/permute_tips/pd_nodes_infer_200_seed_", seed, ".txt"), row.names=F, sep="\t", quote=F)



[1] "Seed = 0, 73% of nodes assigned MNNs, 46.41 cells per node"
[1] "Seed = 1, 73% of nodes assigned MNNs, 46.39 cells per node"
[1] "Seed = 2, 73% of nodes assigned MNNs, 46.52 cells per node"
[1] "Seed = 3, 73% of nodes assigned MNNs, 46.51 cells per node"
[1] "Seed = 4, 73% of nodes assigned MNNs, 46.45 cells per node"
[1] "Seed = 5, 73% of nodes assigned MNNs, 46.44 cells per node"
[1] "Seed = 6, 73% of nodes assigned MNNs, 46.45 cells per node"
[1] "Seed = 7, 73% of nodes assigned MNNs, 46.58 cells per node"
[1] "Seed = 8, 73% of nodes assigned MNNs, 46.53 cells per node"
[1] "Seed = 9, 73% of nodes assigned MNNs, 46.44 cells per node"


pd_orig = read.table(paste0(work_path, "/tree_analysis/impute_nodes/pd_nodes_infer_200.txt"), sep="\t", header=T)
pd_orig = pd_orig %>% filter(cell_id == "internal_node") %>% rename(day_orig = day, celltype_orig = celltype) %>% select(-cell_id)

df_list = list()
for(seed in c(0:9)){
    print(seed)
    pd_mask = read.table(paste0(work_path, "/tree_analysis/impute_nodes/permute_tips/pd_nodes_infer_200_seed_", seed, ".txt"), sep="\t", header=T)
    df_list[[seed + 1]] = pd_mask %>% left_join(pd_orig, by = "node_id") %>% filter(day %in% day_list) %>% mutate(seed = seed)
}
df = do.call(rbind, df_list)

saveRDS(df, paste0(work_path, "/tree_analysis/impute_nodes/permute_tips/internal_nodes_assigned_accuracy.rds"))


df_robust = df %>% filter(celltype_orig != "missing" & celltype_orig == celltype) %>%
    group_by(node_id, day_orig, celltype_orig) %>% tally() %>% filter(n >= 5)
# n = 559
# 559/926852 = 0.0006031168 of internal nodes are identical between permutation and original (>= 5 random permutation)






