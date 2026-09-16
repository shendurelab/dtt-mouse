



##################################################################
### Infer the transcriptional states for individiul internal nodes

source("~/work/scripts/utils.R")
library(dplyr)
library(tidyr)
library(ape)
library(ggplot2)
library(ggrepel)
library(phangorn)

work_path <- "/net/shendure/vol2/projects/cxqiu/work/tapemouse"


###################################################################################
### Step-1: Preparing the data (tree, E13.5 PCA)
### In the integration, assiging day and cell types to each JAX and TapeMouse cells
### pd has the same cell order in the pca coor

pd <- read.csv(paste0(work_path, "/transcriptome_analysis/adata_integration_early.obs.csv"), row.names=1)
pd$cell_id <- rownames(pd)

pd_1 = pd[pd$dataset == "jax",]
pd_jax <- readRDS("/net/shendure/vol2/projects/cxqiu/JAX_rna_mm39/pd.rds")
pd_1_x = pd_1 %>% left_join(pd_jax, by = "cell_id")
pd_1$day = pd_1_x$day
pd_1$major_trajectory = pd_1_x$major_trajectory
pd_1$celltype = pd_1_x$celltype


pd_2 = pd[pd$dataset == "tapemouse",]
cell_meta <- read.table(paste0(work_path, "/tree_analysis/cell_metadata.v8.txt"),
                        header = TRUE, sep = "\t")
pd_2_x = pd_2 %>% left_join(cell_meta, by = "cell_id")
pd_2$day = "E13.5"
pd_2$major_trajectory = pd_2_x$major_trajectory
pd_2$celltype = pd_2_x$celltype

pd_x = rbind(pd_2, pd_1)
sum(pd$cell_id == pd_x$cell_id)

write.table(pd_x[,c("cell_id", "dataset", "day", "major_trajectory", "celltype")], paste0(work_path, "/tree_analysis/impute_nodes/pd.txt"), row.names=T, col.names=T, sep="\t", quote=F)


##########################
# ---- Load tree data ----

tree <- read.tree(paste0(work_path, "/tree_analysis/merged_full_placed.nwk"))

depths <- node.depth.edgelength(tree)

edge <- tree$edge
edge[,1] <- paste0("node_", edge[,1])
edge[,2] <- paste0("node_", edge[,2])

n_tips <- Ntip(tree)
n_int  <- Nnode(tree)

internal <- data.frame(
  node = paste0("node_", (n_tips + 1):(n_tips + n_int)),
  node_depth  = depths[(n_tips + 1):(n_tips + n_int)]
)

# Define bin edges: 0 as its own bin, then 1.5, 1.75, 2.0, ..., 13.25
bin_edges <- c(0, seq(1.5, 13.25, by = 0.25))

# Assign each node to a numeric bin (labeled by the left edge)
internal$bin <- bin_edges[findInterval(internal$node_depth, bin_edges)]

# Format bin as "E8.0", "E8.25", "E8.5", etc.
internal$day <- paste0("E", ifelse(internal$bin == floor(internal$bin),
                                   sprintf("%.1f", internal$bin),
                                   as.character(internal$bin)))

write.table(internal[,c("node", "day")], paste0(work_path, "/tree_analysis/impute_nodes/internal.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(edge, paste0(work_path, "/tree_analysis/impute_nodes/edge.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = c("parent", "child"))


# ---- Impute internal nodes' transcriptome
pd = read.table(paste0(work_path, "/tree_analysis/impute_nodes/pd.txt"), header=T, sep="\t", row.names=1)
pd_E13.5 = pd[pd$cell_id %in% tree$tip.label,]

tips_df <- data.frame(
  node_id = paste0("node_", seq_along(tree$tip.label)),
  cell_id = tree$tip.label,
  stringsAsFactors = FALSE
)

pd_E13.5 = pd_E13.5 %>% left_join(tips_df, by = "cell_id") %>% as.data.frame()

write.table(pd_E13.5, paste0(work_path, "/tree_analysis/impute_nodes/pd_E13.5.tsv"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)






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



# ---- Build combined matrix ----
all_ids = tip_ids + internal_ids
idx_of  = {nid: i for i, nid in enumerate(all_ids)}

all_pcs = np.full((len(all_ids), n_pc), np.nan, dtype=np.float32)
all_pcs[:len(tip_ids)] = pca_impute.to_numpy()

# ---- Fill each internal node = mean of its direct children ----
for cnt, nd in enumerate(internal_ids, 1):
    if cnt % 10000 == 0:
        print(f"{cnt}/{len(internal_ids)}")
    children = children_map.get(nd)
    if children is None:
        continue
    child_idx = [idx_of[c] for c in children]
    all_pcs[idx_of[nd]] = all_pcs[child_idx].mean(axis=0)

# Sanity check
print(np.isnan(all_pcs).any(axis=1).sum(), "rows still contain NaN")

pca_impute_full = pd.DataFrame(all_pcs, index=all_ids, columns=pca_impute.columns)
pca_impute_full.to_csv(f"{work_path}/tree_analysis/impute_nodes/pca_impute_full.tsv", sep="\t")
### It has pca coors for all the internal nodes and tips



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

internal = pd.read_csv(f"{work_path}/tree_analysis/impute_nodes/internal.tsv", sep="\t")
internal.index = internal["node"]
pca_impute = pd.read_csv(f"{work_path}/tree_analysis/impute_nodes/pca_impute_full.tsv", sep="\t", index_col=0)

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

    hits_df.to_csv(f"{work_path}/tree_analysis/impute_nodes/assign/hits_df_{day}_{k}.txt", sep="\t", index=False)





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

for(k in c(200)){
    pd_internal_list = list()
    for(day in day_list){
        pd_internal_list[[day]] = read.table(paste0(work_path, "/tree_analysis/impute_nodes/assign/hits_df_", day, "_", k, ".txt"), header=T, sep='\t')
    }
    pd_internal = do.call(rbind, pd_internal_list)
    
    pd_internal_x = pd_internal %>% left_join(pd[,c("cell_id", "celltype")], by = "cell_id")
    write.table(pd_internal_x, paste0(work_path, "/tree_analysis/impute_nodes/internal_assign_cells_", k, ".txt"), row.names=F, col.names=F, sep="\t", quote=F)
    
    node_cell_num = pd_internal_x %>% group_by(node_id) %>% tally()
    print(paste0('K = ', k, ', ', 100*round(nrow(node_cell_num)/nrow(internal),2), "% of nodes assigned MNNs, ", round(mean(node_cell_num$n), 2), " cells per node"))

    pd_internal_y = pd_internal_x %>% group_by(node_id, celltype) %>% tally() %>% ungroup() %>%
        group_by(node_id) %>% slice_max(order_by = n, n = 1, with_ties = FALSE)
    internal_x = internal %>% select(node_id = node, day) %>% left_join(pd_internal_y[,c("node_id", "celltype")], by = "node_id")
    internal_x$celltype[is.na(internal_x$celltype)] = "missing"
    internal_x$cell_id = "internal_node"
    pd_out = rbind(internal_x, pd_E13.5[,c("node_id", "day", "celltype", "cell_id")])
    write.table(pd_out, paste0(work_path, "/tree_analysis/impute_nodes/pd_nodes_infer_", k, ".txt"), row.names=F, sep="\t", quote=F)
}


[1] "K = 200, 84% of nodes assigned MNNs, 66.6 cells per node"

df = pd_internal_x %>% group_by(node_id) %>% tally()

p = ggplot(df, aes(x = n)) +
    geom_histogram(bins = 30, fill = "#8d99ae", color = "white") +
    labs(x = "n", y = "count") +
    theme_classic()
ggsave("~/share/hist_assignment.pdf", p, width = 3.5, height=5)


x = pd_out %>% filter(cell_id == "internal_node", celltype != "missing") %>%
    group_by(day) %>% tally() %>% 
    left_join(pd_out %>% filter(cell_id == "internal_node") %>%
    group_by(day) %>% tally() %>% rename(total_n = n), by = "day") %>%
    mutate(frac = 100*n/total_n) %>% arrange(frac)

print(sum(x$n)/sum(x$total_n))
#87.9%; 926852/1054879

day         n total_n  frac
E8.5     6905   14930  46.2
E8.75    7692   16753  45.9
E9.0     9503   18647  51.0
E9.25   12560   20927  60.0
E9.5    16778   23731  70.7
E9.75   19826   25548  77.6
E10.0   22851   29781  76.7
E10.25  26607   30871  86.2
E10.5   28046   33525  83.7
E10.75  34294   37133  92.4
E11.0   37344   40018  93.3
E11.25  38885   41192  94.4
E11.5   43363   45318  95.7
E11.75  47803   49578  96.4
E12.0   47551   49040  97.0
E12.25  50363   51505  97.8
E12.5   67436   70198  96.1
E12.75  59633   61838  96.4
E13.0   49526   50988  97.1
E13.25 299886  343358  87.3

E11.5   95.7
E11.75  96.4
E12.0   97.0
E12.25  97.8
E12.5   96.1
E12.75  96.4
E13.0   97.1

###############################################################################
### Step-5: Identifying the trajectories giving rise to each cell type at E13.5


import numpy as np
import pandas as pd
from tqdm import tqdm

work_path = "/net/shendure/vol2/projects/cxqiu/work/tapemouse"

k = 200
pd_node = pd.read_csv(f"{work_path}/tree_analysis/impute_nodes/pd_nodes_infer_{k}.txt", sep="\t")

day_list = ["E13.5","E13.25","E13.0","E12.75","E12.5","E12.25","E12.0",
            "E11.75","E11.5","E11.25","E11.0","E10.75","E10.5",
            "E10.25","E10.0","E9.75","E9.5","E9.25","E9.0",
            "E8.75","E8.5"]
day_set = set(day_list)

edge = pd.read_csv(f"{work_path}/tree_analysis/impute_nodes/edge.tsv", sep="\t")
parent_map = dict(zip(edge['child'], edge['parent']))
day_map = dict(zip(pd_node['node_id'], pd_node['day']))
celltype_map = dict(zip(pd_node['node_id'], pd_node['celltype']))

leaves = pd_node.loc[pd_node['cell_id'] != 'internal_node', 'node_id'].tolist()

records = []
for leaf_id in tqdm(leaves, desc="Tracing ancestors"):
    row = {'leaf_node_id': leaf_id}
    trajectory = []
    d = day_map.get(leaf_id)
    ct = celltype_map.get(leaf_id)
    row[d] = ct
    trajectory.append(ct)
    current = leaf_id
    while current in parent_map:
        current = parent_map[current]
        d = day_map.get(current)
        if d in day_set:
            ct = celltype_map.get(current)
            if pd.isna(ct) or ct == 'missing':
                continue
            row[d] = ct
            if not trajectory or trajectory[-1] != ct:
                trajectory.append(ct)
    row['trajectory'] = '->'.join(reversed(trajectory))
    records.append(row)

trajectory_df = pd.DataFrame(records, columns=['leaf_node_id'] + day_list + ['trajectory'])


top_n = 5

counts = (trajectory_df
          .groupby('E13.5')['trajectory']
          .value_counts()
          .rename('n'))

pct = counts.groupby(level='E13.5').transform(lambda x: 100 * x / x.sum()).map(lambda v: f"{v:.2f}%")

result = pd.concat([counts, pct.rename('pct')], axis=1).reset_index()

top = (result[result['E13.5'] != result['trajectory']]
       .sort_values(['E13.5', 'n'], ascending=[True, False])
       .groupby('E13.5')
       .head(top_n)
       .reset_index(drop=True))

top.to_csv(f"{work_path}/tree_analysis/impute_nodes/top_trajectories_{k}.csv",
           index=False)


###############################################################################
### Step-6: Perform UMAP on the inferred PCA coordinates

import numpy as np
import pandas as pd
import umap
import pickle

work_path = "/net/shendure/vol2/projects/cxqiu/work/tapemouse"

pca_original = pd.read_csv(f"{work_path}/transcriptome_analysis/adata_integration_early.pca.csv", index_col=0)
pd_original = pd.read_csv(f"{work_path}/transcriptome_analysis/adata_integration_early.obs.csv", index_col=0)
pd_original_add_more = pd.read_csv(f"{work_path}/tree_analysis/impute_nodes/pd.txt", sep="\t", index_col=0)
pd_original.index.equals(pd_original_add_more.index)
pd_original['day'] = pd_original_add_more['day']
pca_original.index = pd_original.index

node_assign = pd.read_csv(f"{work_path}/tree_analysis/impute_nodes/internal_assign_cells_200.txt", sep = "\t")
node_assign.columns = ['cell_id', 'node', 'day', 'celltype']
merged = node_assign.merge(pca_original, left_on='cell_id', right_index=True, how='inner')
pc_cols = pca_original.columns.tolist()
node_pca = merged.groupby('node')[pc_cols].mean()
### n = 926852 internal nodes

node_pca.to_csv(f"{work_path}/tree_analysis/impute_nodes/node_pca.csv",
           index=True)

# --- Filter to jax cells ---
mask_jax = pd_original['dataset'] == 'jax'
pd_sub = pd_original.loc[mask_jax].copy()
pca_sub = pca_original.loc[mask_jax].copy()   # or: pca_original.iloc[mask_jax.values]

# --- Cap each day at 50,000 cells ---
rng = np.random.default_rng(0)
idx_keep = (
    pd_sub.groupby('day', group_keys=False)
          .apply(lambda x: x.sample(n=min(len(x), 50000), random_state=0))
          .index
)

pd_sub  = pd_sub.loc[idx_keep]
pca_sub = pca_sub.loc[idx_keep]

# --- Sanity check ---
print(pd_sub.shape, pca_sub.shape)
print(pd_sub['day'].value_counts().sort_index())
# n = 1046253 cells

reducer = umap.UMAP(
    n_components   = 2,
    n_neighbors    = 30,
    min_dist       = 0.3,
    metric         = 'euclidean',
    random_state   = 0,
    verbose        = True,
)

embedding = reducer.fit_transform(pca_sub.values)

pd_sub['cell_id'] = pd_sub.index
umap_df = pd_sub[['cell_id', 'day']].copy()
umap_df['UMAP_1'] = embedding[:, 0]
umap_df['UMAP_2'] = embedding[:, 1]

umap_df.to_csv(f"{work_path}/tree_analysis/impute_nodes/umap/UMAP_backbone.csv",
           index=True)

with open(f"{work_path}/tree_analysis/impute_nodes/umap/umap_model.pkl", "wb") as f:
    pickle.dump(reducer, f)

new_embedding = reducer.transform(node_pca.values)

new_embedding_df = pd.DataFrame(
    new_embedding,
    index=node_pca.index,
    columns=['UMAP_1', 'UMAP_2']
)

new_embedding_df.to_csv(f"{work_path}/tree_analysis/impute_nodes/umap/UMAP_new.csv",
           index=True)



### PLOT

source("~/work/scripts/utils.R")
work_path <- "/net/shendure/vol2/projects/cxqiu/work/tapemouse"

df = read.table(paste0(work_path, "/tree_analysis/impute_nodes/pd_nodes_infer_200.txt"), sep="\t", header=T)

umap_coor = read.csv(paste0(work_path, "/tree_analysis/impute_nodes/umap/UMAP_new.csv"))
colnames(umap_coor) = c("node_id", "UMAP_1", "UMAP_2")

df = df %>% left_join(umap_coor, by = "node_id") %>% filter(!is.na(UMAP_1))
### 926852 internal nodes assigned

df$day = factor(df$day, level = names(day_color_plate))
table(df$day)

major_trajectory_table = read.table(paste0(work_path, "/tree_analysis/major_trajectory_celltype_table.txt"), sep="\t", header=T)
df = df %>% left_join(major_trajectory_table, by = "celltype")
### 926852 cells

pd_jax = readRDS("/net/shendure/vol2/projects/cxqiu/JAX_rna_mm39/pd.rds")
umap_coor = read.csv(paste0(work_path, "/tree_analysis/impute_nodes/umap/UMAP_backbone.csv"))
df_backbone = umap_coor %>% left_join(pd_jax[,c("celltype", "major_trajectory", "cell_id")], by = "cell_id")
### n = 1046253 cells




p = ggplot() +
    geom_point(data = df_backbone, aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
    geom_point(data = df, aes(x = UMAP_1, y = UMAP_2, color = major_trajectory), size=0.15) +
    theme_void() +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=major_trajectory_color_plate)
ggsave("~/share/Fig7_umap_major_trajectory_new.png", p, dpi = 300, height = 5, width = 5)

p = ggplot() +
    geom_point(data = df, aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
    geom_point(data = df_backbone, aes(x = UMAP_1, y = UMAP_2, color = major_trajectory), size=0.15) +
    theme_void() +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=major_trajectory_color_plate)
ggsave("~/share/Fig7_umap_major_trajectory_jax.png", p, dpi = 300, height = 5, width = 5)



df_sub = df %>% group_by(day) %>% slice_sample(n = 50000)
p = ggplot() +
    geom_point(data = df_backbone, aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
    geom_point(data = df_sub[sample(1:nrow(df_sub)),], aes(x = UMAP_1, y = UMAP_2, color = day), size=0.15) +
    theme_void() +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=day_color_plate)
ggsave("~/share/Fig7_umap_day_new.png", p, dpi = 300, height = 5, width = 5)

p = ggplot() +
    geom_point(data = df, aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
    geom_point(data = df_backbone[sample(1:nrow(df_backbone)),], aes(x = UMAP_1, y = UMAP_2, color = day), size=0.15) +
    theme_void() +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=day_color_plate)
ggsave("~/share/Fig7_umap_day_jax.png", p, dpi = 300, height = 5, width = 5)


### Three important profiles:
### internal_assign_cells_200.txt: MNN pairs between internal nodes and JAX cells
### pd_nodes_infer_200.txt: internal nodes and tips' cell types
### pca_impute_full.tsv: internal nodes and tips' PCs, directly imputed by descendants
### node_pca.csv: internal nodes' PCs, after updating by its MNNs



