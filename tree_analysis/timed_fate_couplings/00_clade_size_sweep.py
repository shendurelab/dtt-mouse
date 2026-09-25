#!/usr/bin/env python3
"""STEP 4 -- sweep the clade-size cap K, to place cell-type couplings in developmental time.

Same maximal-terminal-proximate-clade definition and same exact closed-form tip-label-permutation
null as c1_cooccur.py; only K varies. Each K carries a median clade-ancestor age, so the sweep
converts a structural parameter into a developmental-time axis: small K -> recent, terminal-
proximate clades; large K -> ancient clades approaching the whole blastomere.

Writes out/c4_sweep.npz with, for every K and blastomere, the full 82x82 obs/exp/log2-enrichment
matrices plus clade count, median ancestor age and median clade size.
Deterministic (no RNG anywhere).
"""
import gzip
import numpy as np, gzip, csv, collections, json, os, sys, time
from scipy.special import gammaln

import os as _os
_REPO = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))
def _support(name, env):
    """Repo-relative input, overridable with an env var."""
    return _os.environ.get(env, _os.path.join(_REPO, "support_data", name))


HERE  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TREE  = "/Users/jay.shendure/Dropbox/claude/top_to_bottom_phylogeny/out/mergedtree_dttpq_v8.npz"
META  = _support("cell_metadata.v8.txt.gz", "DTT_CELL_METADATA")
ROUTE = _support("e3v8.routing_labels.tsv.gz", "DTT_ROUTING_LABELS")
MIN_CELLS = 100
KS = [3, 5, 8, 12, 15, 20, 30, 45, 65, 100, 150, 200, 300, 450, 700, 1000, 1500, 2500]
log = lambda *a: print(*a, file=sys.stderr, flush=True)
t0 = time.time()

z = np.load(TREE, allow_pickle=True)
parent, node_time, is_leaf, names = z["parent"].astype(np.int64), z["time"], z["is_leaf"], z["names"]
N = len(parent); leaf_node = np.flatnonzero(is_leaf)
leaf_ids = np.array([names[u] for u in leaf_node], object)
pp = parent.copy(); pp[0] = 0; internal = ~is_leaf
depth = np.zeros(N, np.int32); todo = np.flatnonzero(parent >= 0)
while True:
    nd = depth[parent[todo]] + 1
    if np.array_equal(nd, depth[todo]): break
    depth[todo] = nd
o = np.argsort(depth, kind="stable")
db = np.flatnonzero(np.r_[True, depth[o][1:] != depth[o][:-1], True])
levels = [o[db[i]:db[i + 1]] for i in range(len(db) - 1)]

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
log(f"{T} cell types; A={int((leaf_bl==1).sum()):,} B={int((leaf_bl==2).sum()):,}; "
    f"K = {KS}")

def subtree_size(leafmask):
    sz = np.zeros(N, np.int64); sz[leaf_node[leafmask]] = 1
    for d in range(len(levels) - 1, 0, -1):
        nd = levels[d]; np.add.at(sz, pp[nd], sz[nd])
    return sz

def clade_of(sel):
    anc = np.where(sel, np.arange(N), pp).astype(np.int64)
    for _ in range(80):
        nx = anc[anc]
        if np.array_equal(nx, anc): break
        anc = nx
    return anc

def cooccur(bmask, K, size_b):
    sel = internal & (size_b >= 3) & (size_b <= K) & (size_b[pp] > K)
    seln = np.flatnonzero(sel); ncl = len(seln)
    if ncl == 0: return None
    covmask = bmask & (leaf_type >= 0)
    cov_node = leaf_node[covmask]; cov_type = leaf_type[covmask]
    anc = clade_of(sel); cr = anc[cov_node]; keep = sel[cr]
    ucl, cd = np.unique(cr[keep], return_inverse=True); tyk = cov_type[keep]
    ncl_cov = len(ucl); Ncov = int(cd.size)
    pres = np.zeros((ncl_cov, T), bool); pres[cd, tyk] = True
    Pm = pres.astype(np.float64); obs = Pm.T @ Pm
    csize = np.bincount(cd, minlength=ncl_cov)
    su, sc = np.unique(csize, return_counts=True); sc = sc.astype(float); su = su.astype(np.int64)
    ncnt = np.bincount(tyk, minlength=T); LGN = gammaln(np.arange(Ncov + 2)); Fm = {}
    def F(m):
        if m in Fm: return Fm[m]
        if m >= Ncov: Fm[m] = 0.0; return 0.0
        ok = (Ncov - m - su) >= 0; val = np.zeros(su.shape); s_ = su[ok]
        val[ok] = np.exp(LGN[Ncov - m] - LGN[Ncov - m - s_] - LGN[Ncov] + LGN[Ncov - s_])
        Fm[m] = float(np.dot(sc, val)); return Fm[m]
    Fn = np.array([F(int(v)) for v in ncnt])
    exp = np.empty((T, T))
    for a in range(T):
        for b in range(a, T):
            exp[a, b] = exp[b, a] = max(ncl_cov - Fn[a] - Fn[b] + F(int(ncnt[a] + ncnt[b])), 0.0)
    return dict(obs=obs, exp=exp, enr=np.log2((obs + 1) / (exp + 1)), ncl=ncl,
                med_age=float(np.median(node_time[seln])), med_size=float(np.median(size_b[seln])),
                mean_size=float(size_b[seln].mean()), ncl_cov=ncl_cov)

out = {}
rows = []
print(f"{'K':>6}  {'blast':5}{'#clades':>10}{'med anc age':>13}{'med size':>10}")
for bl, blab in [(1, "B1"), (2, "B2")]:
    size_b = subtree_size(leaf_bl == bl)
    for K in KS:
        R = cooccur(leaf_bl == bl, K, size_b)
        if R is None: continue
        for k in ("obs", "exp", "enr"): out[f"{blab}|{K}|{k}"] = R[k]
        rows.append(dict(blast=blab, K=K, ncl=R["ncl"], med_age=R["med_age"],
                         med_size=R["med_size"], mean_size=R["mean_size"]))
        print(f"{K:>6}  {blab:5}{R['ncl']:>10,}{('E%.2f'%R['med_age']):>13}{R['med_size']:>10.0f}")
    log(f"  {blab} done [{time.time()-t0:.0f}s]")

np.savez_compressed(f"{HERE}/out/c4_sweep.npz", types=np.array(universe, object),
                    KS=np.array(KS), **out)
json.dump(rows, open(f"{HERE}/out/c4_sweep_meta.json", "w"), indent=1)
print(f"\nwrote out/c4_sweep.npz, out/c4_sweep_meta.json [{time.time()-t0:.0f}s]")
print("\nK -> median clade ancestor age (E-day), the developmental-time axis:")
for K in KS:
    a = [r for r in rows if r["K"] == K and r["blast"] == "B1"]
    b = [r for r in rows if r["K"] == K and r["blast"] == "B2"]
    if a and b:
        print(f"  K={K:>5}:  A E{a[0]['med_age']:.2f} ({a[0]['ncl']:>7,} clades)   "
              f"B E{b[0]['med_age']:.2f} ({b[0]['ncl']:>7,} clades)")
