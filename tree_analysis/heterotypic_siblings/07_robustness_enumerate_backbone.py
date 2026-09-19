#!/usr/bin/env python3
"""STEP 1 -- enumerate ALL tree-sibling pairs on the FULL placed tree, with the
label-independent evidence attached to each.

Tree siblings are defined exactly as in the preceding manuscript subsection ("Tree siblings
share fate far in excess of chance"): for every internal node take its DIRECT LEAF-CHILDREN
and enumerate all unordered pairs among them. Cherries (k=2) give 1 pair; terminal polytomies
(k>=3) give C(k,2); leaves whose only siblings are internal nodes ("loners", k=1) give none --
and give none under the permutation null either, so they are treated fairly.

ALL sibling pairs are kept, not just the currently-heterotypic ones: the null permutes cell-type
labels over the fixed tree, so any pair can become heterotypic under a permutation. The
per-pair evidence recorded here is label-INDEPENDENT and therefore valid for observed and
permuted labels alike:
  age    -- MRCA age = date of the shared parent node (E-days)
  disc   -- # DNA-typewriter sites called in both tips that DISAGREE ("edits separating them")
  shared -- # sites called in both that agree and are edited (>=1)
  ncomp  -- # sites called in both tips
  kgroup -- size of the terminal sibling group the pair came from

Output: out/h1_pairs.npz.  Read-only on inputs. Deterministic.
"""
import numpy as np, csv, gzip, os, sys, time
from itertools import combinations

HERE  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SUF   = os.environ.get("RUN_SUF","")
TREE  = os.environ.get("RUN_TREE",
        "/Users/jay.shendure/Dropbox/claude/top_to_bottom_phylogeny/out/mergedtree_dttpq_v8.npz")
META  = "/Users/jay.shendure/Dropbox/claude/current/final_push/data/cell_metadata.v8.txt"
GENO  = "/Users/jay.shendure/Dropbox/claude/top_to_bottom_phylogeny/out/newcode_v8_all.npz"
ROUTE = "/Users/jay.shendure/Dropbox/claude/current/final_push/fig6_v8/e3v8.routing_labels.tsv.gz"
PM    = "/Users/jay.shendure/Dropbox/claude/current/final_push/fig6_v8/celltype_postmitotic_v8.csv"
OUT   = f"{HERE}/out/h1_pairs{SUF}.npz"
t0 = time.time()
log = lambda *a: print(*a, file=sys.stderr, flush=True)

# ---- tree ---------------------------------------------------------------------------------
z = np.load(TREE, allow_pickle=True)
par, tm, isleaf, names = z["parent"], z["time"], z["is_leaf"], z["names"]
L = np.flatnonzero(isleaf)
lid = np.array([names[u] for u in L], object)
nL = len(L)
log(f"tree {os.path.basename(TREE)}: {len(par):,} nodes, {nL:,} leaves, tips at E{tm[isleaf].max():.2f}, "
    f"lineages@E6.5={int(((tm[par[par>=0]]<6.5)&(tm[par>=0]>=6.5)).sum())}")

# ---- per-leaf labels ----------------------------------------------------------------------
ct, mtj = {}, {}
with open(META) as f:
    r = csv.reader(f, delimiter="\t"); h = next(r)
    ci, cei, mi = h.index("cell_id"), h.index("celltype"), h.index("major_trajectory")
    for row in r:
        if len(row) > cei: ct[row[ci]] = row[cei]; mtj[row[ci]] = row[mi]
route = {}
with gzip.open(ROUTE, "rt") as f:
    f.readline()
    for line in f:
        p = line.rstrip("\n").split("\t"); route[p[0]] = p[1]
g = np.load(GENO, allow_pickle=True)
code = g["code"]
grow = {c: i for i, c in enumerate(g["cell_ids"])}

POST = [r["celltype"] for r in csv.DictReader(open(PM)) if r["classification"] == "post_mitotic"]
POSTS = set(POST)
types = sorted(set(ct.values()))
tix = {t: i for i, t in enumerate(types)}
leaf_ct = np.array([tix.get(ct.get(x, ""), -1) for x in lid], np.int16)
leaf_gr = np.array([grow.get(x, -1) for x in lid], np.int64)
leaf_bl = np.array([route.get(x, "NA") for x in lid], object)
type_is_pm = np.array([t in POSTS for t in types])
log(f"{len(types)} cell types ({len(POST)} post-mitotic); annotated leaves "
    f"{int((leaf_ct>=0).sum()):,}; genotyped {int((leaf_gr>=0).sum()):,}")

# ---- terminal sibling groups: direct leaf-children of each internal node -------------------
lp = par[L]
order = np.argsort(lp, kind="stable")
lp_s = lp[order]
bnd = np.flatnonzero(np.r_[True, lp_s[1:] != lp_s[:-1], True])
gsz = np.diff(bnd)
log(f"sibling groups {len(gsz):,}: k=1 {int((gsz==1).sum()):,} | k=2 {int((gsz==2).sum()):,} | "
    f"k>=3 {int((gsz>=3).sum()):,} (max k={gsz.max()})")

I, J, NODE, KG = [], [], [], []
for b in np.flatnonzero(gsz >= 2):
    mem = order[bnd[b]:bnd[b + 1]]
    node = lp_s[bnd[b]]; k = len(mem)
    for a, c in combinations(mem, 2):
        I.append(a); J.append(c); NODE.append(node); KG.append(k)
I = np.array(I, np.int32); J = np.array(J, np.int32)
NODE = np.array(NODE, np.int64); KG = np.array(KG, np.int16)
log(f"tree-sibling pairs: {len(I):,} ({int((KG==2).sum()):,} from cherries, "
    f"{int((KG>=3).sum()):,} from terminal polytomies) [{time.time()-t0:.0f}s]")

# ---- label-independent evidence -----------------------------------------------------------
AGE = tm[NODE].astype(np.float32)
ga, gb = leaf_gr[I], leaf_gr[J]
have = (ga >= 0) & (gb >= 0)
DISC = np.full(len(I), -1, np.int16); SHAR = np.full(len(I), -1, np.int16); NCMP = np.full(len(I), -1, np.int16)
CH = 200_000
for s in range(0, len(I), CH):
    e = min(s + CH, len(I)); m = have[s:e]
    if not m.any(): continue
    A = code[ga[s:e][m]].reshape(int(m.sum()), -1); B = code[gb[s:e][m]].reshape(int(m.sum()), -1)
    comp = (A >= 0) & (B >= 0)
    NCMP[s:e][m] = comp.sum(1); DISC[s:e][m] = (comp & (A != B)).sum(1)
    SHAR[s:e][m] = (comp & (A == B) & (A >= 1)).sum(1)
log(f"edit stats for {int(have.sum()):,}/{len(I):,} pairs [{time.time()-t0:.0f}s]")

blI = np.array([leaf_bl[i] for i in I], object)
blJ = np.array([leaf_bl[j] for j in J], object)
log(f"pairs whose members disagree on blastomere: {int((blI!=blJ).sum()):,}")

np.savez_compressed(OUT,
    i=I, j=J, node=NODE, kgroup=KG, age=AGE, disc=DISC, shared=SHAR, ncomp=NCMP,
    pair_blast=blI.astype(str),
    leaf_ct=leaf_ct, leaf_blast=leaf_bl.astype(str),
    types=np.array(types, object), type_is_pm=type_is_pm,
    post_types=np.array(POST, object))
log(f"wrote {OUT} [{time.time()-t0:.0f}s]")

# ---- how many pairs are heterotypic under the OBSERVED labels -----------------------------
pmI = type_is_pm[leaf_ct[I]]; pmJ = type_is_pm[leaf_ct[J]]
het = pmI ^ pmJ
log(f"\nobserved: {int(het.sum()):,} heterotypic sibling pairs "
    f"({100*het.mean():.1f}% of all sibling pairs); "
    f"{int((pmI&pmJ).sum()):,} post-mitotic|post-mitotic; {int((~pmI&~pmJ).sum()):,} prog|prog")
pmleaf = np.where(pmI, I, J)[het]
log(f"distinct post-mitotic cells with >=1 progenitor sibling: {len(np.unique(pmleaf)):,}")
