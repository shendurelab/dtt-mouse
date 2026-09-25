#!/usr/bin/env python3
"""STEP 1 -- clade-level co-occurrence on the FINAL tree, at K = 15 and K = 200.

Faithful port of draft_hierarchy/clade_cooccurrence/code/hv6_01_cooccur.py (itself a port of
step4_clade_coincidence.R) onto the final full placed tree: 1,340,794 tips, dated so 658
lineages exist at E6.5. Method unchanged; only the tree and the K values differ (K=200 replaces
the earlier K=60).

MAXIMAL TERMINAL-PROXIMATE CLADES. For each blastomere separately (A=B1, B=B2), subtree sizes
are counted over THAT blastomere's leaves only. An internal node is a maximal clade if
    3 <= size <= K   and   parent size > K
i.e. it cannot be grown further without breaching the cap.

CO-OCCURRENCE. obs[i,j] = number of maximal clades containing at least one cell of type i and
at least one of type j (presence per clade, multiplicities collapsed).

NULL. Tip-label permutation, evaluated EXACTLY in closed form (the analytic n_perm -> inf
limit, so no Monte-Carlo error):
    exp[a,b] = ncl - F(n_a) - F(n_b) + F(n_a + n_b),   F(m) = sum_c C(N-m, s_c) / C(N, s_c)
which is the expected number of clades containing both types when labels are shuffled over the
N covered tips and clade sizes s_c are held fixed. enr = log2((obs+1)/(exp+1)).

UNIVERSE. Cell types with >= MIN_CELLS (100) cells in EACH blastomere subtree.

Writes out/cooccur_{B1,B2}_K{15,200}_final.csv and out/c1_clades.json.
Read-only on inputs. Deterministic (no RNG anywhere).
"""
import gzip
import numpy as np, gzip, csv, collections, json, os, sys, time
from scipy.special import gammaln

import os as _os
_REPO = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))
def _support(name, env=None, hint=None):
    """Resolve an input under support_data/, overridable by env var."""
    if env:
        v = _os.environ.get(env)
        if v:
            return v
    p = _os.path.join(_REPO, "support_data", name)
    if not _os.path.exists(p) and hint:
        raise FileNotFoundError(f"{p} not found. {hint}")
    return p


HERE  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TREE  = _support("mergedtree_dttpq_v8.npz", "DTT_MERGED_TREE_NPZ", 'Generate it with: python3 tools/make_merged_tree_npz.py (derived from support_data/merged_full_placed.nwk).')
META  = _support("cell_metadata.v8.txt.gz", "DTT_CELL_METADATA")
ROUTE = _support("e3v8.routing_labels.tsv.gz", "DTT_ROUTING_LABELS")
MIN_CELLS = 100
KS = (15, 200)
os.makedirs(f"{HERE}/out", exist_ok=True)
log = lambda *a: print(*a, file=sys.stderr, flush=True)
t0 = time.time()

# ---- tree ---------------------------------------------------------------------------------
z = np.load(TREE, allow_pickle=True)
parent, node_time, is_leaf, names = z["parent"].astype(np.int64), z["time"], z["is_leaf"], z["names"]
N = len(parent)
leaf_node = np.flatnonzero(is_leaf)
leaf_ids = np.array([names[u] for u in leaf_node], object)
pp = parent.copy(); pp[0] = 0
internal = ~is_leaf
log(f"tree {os.path.basename(TREE)}: {N:,} nodes, {len(leaf_node):,} tips, "
    f"lineages@E6.5={int(((node_time[parent[parent>=0]]<6.5)&(node_time[parent>=0]>=6.5)).sum())}")

# Group nodes by depth so subtree sums accumulate level-by-level (vectorised per level).
# A genuine topological order is required: node_time cannot substitute for it, because
# zero-length branches give a parent and child identical times and an arbitrary tie-break
# would accumulate a child after its parent.
depth = np.zeros(N, np.int32)
todo = np.flatnonzero(parent >= 0)
while True:                      # depth[u] = depth[parent]+1, iterated to convergence
    nd = depth[parent[todo]] + 1
    if np.array_equal(nd, depth[todo]): break
    depth[todo] = nd
maxd = int(depth.max())
o = np.argsort(depth, kind="stable")
db = np.flatnonzero(np.r_[True, depth[o][1:] != depth[o][:-1], True])
levels = [o[db[i]:db[i + 1]] for i in range(len(db) - 1)]
log(f"levels: {len(levels)} (max depth {maxd}) [{time.time()-t0:.0f}s]")

def subtree_size(leafmask):
    """# of that blastomere's leaves below each node."""
    sz = np.zeros(N, np.int64)
    sz[leaf_node[leafmask]] = 1
    for d in range(len(levels) - 1, 0, -1):
        nd = levels[d]
        np.add.at(sz, pp[nd], sz[nd])
    return sz

def clade_of(sel):
    """map each node to its nearest selected ancestor-or-self (pointer jumping)."""
    anc = np.where(sel, np.arange(N), pp).astype(np.int64)
    for _ in range(80):
        nx = anc[anc]
        if np.array_equal(nx, anc): break
        anc = nx
    return anc

# ---- per-leaf cell type + blastomere -------------------------------------------------------
ct_of = {}
with gzip.open(META, "rt") as f:
    r = csv.reader(f, delimiter="\t"); h = next(r); ci, cei = h.index("cell_id"), h.index("celltype")
    for row in r:
        if len(row) > cei: ct_of[row[ci]] = row[cei]
bl_of = {}
with gzip.open(ROUTE, "rt") as f:
    f.readline()
    for line in f:
        p = line.rstrip("\n").split("\t"); bl_of[p[0]] = p[1]
leaf_ct = np.array([ct_of.get(x, "NA") for x in leaf_ids], object)
leaf_bl = np.array([1 if bl_of.get(x) == "B1" else (2 if bl_of.get(x) == "B2" else 0)
                    for x in leaf_ids], np.int8)
cntA = collections.Counter(leaf_ct[leaf_bl == 1]); cntB = collections.Counter(leaf_ct[leaf_bl == 2])
universe = sorted([t for t in cntA if t != "NA" and cntA[t] >= MIN_CELLS and cntB.get(t, 0) >= MIN_CELLS])
tix = {t: i for i, t in enumerate(universe)}; T = len(universe)
leaf_type = np.array([tix.get(t, -1) for t in leaf_ct], np.int64)
log(f"A={int((leaf_bl==1).sum()):,} B={int((leaf_bl==2).sum()):,} tips; "
    f"universe {T} cell types (>={MIN_CELLS} cells in each blastomere)")

# ---- co-occurrence -------------------------------------------------------------------------
def cooccur(bmask, K, size_b):
    size_par = size_b[pp]
    sel = internal & (size_b >= 3) & (size_b <= K) & (size_par > K)
    seln = np.flatnonzero(sel); ncl = len(seln)
    med_age = float(np.median(node_time[seln])) if ncl else float("nan")
    sizes_all = size_b[seln]
    covmask = bmask & (leaf_type >= 0)
    cov_node = leaf_node[covmask]; cov_type = leaf_type[covmask]
    anc = clade_of(sel); cr = anc[cov_node]; keep = sel[cr]
    crk = cr[keep]; tyk = cov_type[keep]
    ucl, cd = np.unique(crk, return_inverse=True); ncl_cov = len(ucl)
    Ncov = int(cd.size)
    pres = np.zeros((ncl_cov, T), bool); pres[cd, tyk] = True
    Pm = pres.astype(np.float64)
    obs = Pm.T @ Pm
    csize = np.bincount(cd, minlength=ncl_cov)
    su, sc = np.unique(csize, return_counts=True); sc = sc.astype(float); su = su.astype(np.int64)
    ncnt = np.bincount(tyk, minlength=T)
    LGN = gammaln(np.arange(Ncov + 2))
    Fmemo = {}
    def F(m):
        if m in Fmemo: return Fmemo[m]
        if m >= Ncov: Fmemo[m] = 0.0; return 0.0
        ok = (Ncov - m - su) >= 0; val = np.zeros(su.shape)
        s_ = su[ok]; val[ok] = np.exp(LGN[Ncov - m] - LGN[Ncov - m - s_] - LGN[Ncov] + LGN[Ncov - s_])
        r = float(np.dot(sc, val)); Fmemo[m] = r; return r
    Fn = np.array([F(int(v)) for v in ncnt])
    exp = np.empty((T, T))
    for a in range(T):
        for b in range(a, T):
            e = max(ncl_cov - Fn[a] - Fn[b] + F(int(ncnt[a] + ncnt[b])), 0.0)
            exp[a, b] = exp[b, a] = e
    enr = np.log2((obs + 1) / (exp + 1))
    return dict(obs=obs, exp=exp, enr=enr, ncl=ncl, ncl_cov=ncl_cov, med_age=med_age,
                Ncov=Ncov, mean_size=float(sizes_all.mean()), med_size=float(np.median(sizes_all)))

report = {}
print(f"\n{'blast':6}{'K':>5}{'#clades':>11}{'med anc age':>13}{'med size':>10}{'mean size':>11}{'covered tips':>14}")
for bl, blab in [(1, "B1"), (2, "B2")]:
    bmask = leaf_bl == bl
    size_b = subtree_size(bmask)
    for K in KS:
        R = cooccur(bmask, K, size_b)
        with open(f"{HERE}/out/cooccur_{blab}_K{K}_final.csv", "w", newline="") as f:
            w = csv.writer(f); w.writerow(["celltype_1", "celltype_2", "obs", "exp", "log2_enr"])
            for a in range(T):
                for b in range(T):
                    w.writerow([universe[a], universe[b], int(R["obs"][a, b]),
                                round(float(R["exp"][a, b]), 4), round(float(R["enr"][a, b]), 6)])
        report[f"{blab}_K{K}"] = {k: R[k] for k in ("ncl", "ncl_cov", "med_age", "Ncov", "mean_size", "med_size")}
        print(f"{blab:6}{K:>5}{R['ncl']:>11,}{('E%.2f'%R['med_age']):>13}"
              f"{R['med_size']:>10.0f}{R['mean_size']:>11.1f}{R['Ncov']:>14,}")
        log(f"  done {blab} K={K} [{time.time()-t0:.0f}s]")

json.dump({"universe": universe, "n_types": T, "min_cells": MIN_CELLS,
           "tree": os.path.basename(TREE), "n_tips": int(is_leaf.sum()), "report": report},
          open(f"{HERE}/out/c1_clades.json", "w"), indent=1)
print(f"\nwrote out/cooccur_{{B1,B2}}_K{{{','.join(map(str,KS))}}}_final.csv, out/c1_clades.json")
print("\n=== manuscript placeholders ===")
print(f"universe: {T} cell types (>= {MIN_CELLS} cells in each blastomere subtree)")
for K in KS:
    a, b = report[f"B1_K{K}"], report[f"B2_K{K}"]
    print(f"  K = {K}: {a['ncl']:,} (A) and {b['ncl']:,} (B) clades; "
          f"median ancestor age E{a['med_age']:.1f} (A) and E{b['med_age']:.1f} (B)")
