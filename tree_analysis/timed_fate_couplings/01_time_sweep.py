#!/usr/bin/env python3
"""STEP 7 -- the TIME SWEEP: clades defined by the age of their dated ancestor.

Used only for the sweep / timing analysis. The K=15 and K=200 analyses (Fig 6C-D) are unchanged
and stay in c1_cooccur.py.

For each developmental time t on a 0.25-day grid, take every lineage crossing t -- node u with
age(parent(u)) < t <= age(u) -- and let its clade be all descendant tips. This partitions the tips
exactly, every clade has ancestor age t by construction, and the clade count is the
lineages-through-time curve. Requires >= 3 sampled tips of the blastomere being analysed.
Co-occurrence and the exact closed-form tip-label-permutation null are unchanged.

DIVERGENCE AGE. For each cell-type pair, the earliest (deepest) time slice at which the pair is
still significantly co-enriched. This is a direct estimate of how far back the two fates share
lineage, replacing the indirect "date of the ancestral nodes at which the signal was concentrated".

The key test of whether the approach works is REPRODUCIBILITY: the divergence age is estimated
independently in blastomere A and in blastomere B and the two are compared.

Writes out/c7_timesweep.npz, out/c7_divergence_ages.csv, out/c7_carpet.png, out/c7_divmatrix.png,
out/c7_ab_divage.png.  Deterministic (no RNG).
"""
import gzip
import numpy as np, gzip, csv, collections, json, os, sys, time
from scipy.special import gammaln
from scipy.stats import spearmanr, pearsonr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm, Normalize
from scipy.cluster.hierarchy import linkage, leaves_list, fcluster
from scipy.spatial.distance import squareform

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
CK    = _os.path.join(_REPO, "support_data")
MIN_CELLS, MIN_CLADE, Z_THR, MIN_OBS = 100, 3, 3.0, 5
TIMES = [round(7.0 + 0.25 * i, 2) for i in range(26)]          # E7.00 .. E13.25
INK, INK2, MUTED = "#1f2328", "#4a5158", "#8b9198"
DIV = LinearSegmentedColormap.from_list("so", [
    (0.00, "#4dd2ff"), (0.25, "#1f6feb"), (0.50, "#000000"), (0.75, "#ff7a1a"), (1.00, "#ffd400")])
SEQ = LinearSegmentedColormap.from_list("age", [
    (0.00, "#4a1330"), (0.30, "#a63f2e"), (0.65, "#e0a13c"), (1.00, "#f2e6c9")])
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK})
log = lambda *a: print(*a, file=sys.stderr, flush=True)
t0 = time.time()

z = np.load(TREE, allow_pickle=True)
parent, node_time, is_leaf, names = z["parent"].astype(np.int64), z["time"], z["is_leaf"], z["names"]
N = len(parent); leaf_node = np.flatnonzero(is_leaf)
leaf_ids = np.array([names[u] for u in leaf_node], object)
pp = parent.copy(); pp[0] = 0
depth = np.zeros(N, np.int32); todo = np.flatnonzero(parent >= 0)
while True:
    nd = depth[parent[todo]] + 1
    if np.array_equal(nd, depth[todo]): break
    depth[todo] = nd
o_ = np.argsort(depth, kind="stable")
db = np.flatnonzero(np.r_[True, depth[o_][1:] != depth[o_][:-1], True])
levels = [o_[db[i]:db[i + 1]] for i in range(len(db) - 1)]

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
cA = collections.Counter(leaf_ct[leaf_bl == 1]); cB = collections.Counter(leaf_ct[leaf_bl == 2])
universe = sorted([t for t in cA if t != "NA" and cA[t] >= MIN_CELLS and cB.get(t, 0) >= MIN_CELLS])
tix = {t: i for i, t in enumerate(universe)}; T = len(universe)
leaf_type = np.array([tix.get(t, -1) for t in leaf_ct], np.int64)
log(f"{T} cell types; sweeping {len(TIMES)} time slices E{TIMES[0]}..E{TIMES[-1]}")

def subtree_size(m):
    sz = np.zeros(N, np.int64); sz[leaf_node[m]] = 1
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

def enrich(sel, bmask):
    covmask = bmask & (leaf_type >= 0)
    cov_node = leaf_node[covmask]; cov_type = leaf_type[covmask]
    anc = clade_of(sel); cr = anc[cov_node]; keep = sel[cr]
    ucl, cd = np.unique(cr[keep], return_inverse=True); tyk = cov_type[keep]
    ncl_cov = len(ucl); Ncov = int(cd.size)
    if ncl_cov < 2: return None
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
    return dict(obs=obs, exp=exp, enr=np.log2((obs + 1) / (exp + 1)),
                ncl=int(sel.sum()), ncl_cov=ncl_cov, tips=int(csize.sum()))

store = {}; meta = []
for bl, blab in [(1, "B1"), (2, "B2")]:
    bmask = leaf_bl == bl; size_b = subtree_size(bmask)
    for t in TIMES:
        cross = (parent >= 0) & (node_time[pp] < t) & (node_time >= t)
        sel = cross & (size_b >= MIN_CLADE)
        R = enrich(sel, bmask)
        if R is None: continue
        for k in ("obs", "exp", "enr"): store[f"{blab}|{t}|{k}"] = R[k]
        meta.append(dict(blast=blab, t=t, lineages=int(cross.sum()), clades=R["ncl"],
                         clades_cov=R["ncl_cov"], tips=R["tips"]))
    log(f"  {blab} done [{time.time()-t0:.0f}s]")
np.savez_compressed(f"{HERE}/out/c7_timesweep.npz", types=np.array(universe, object),
                    TIMES=np.array(TIMES), **store)
json.dump(meta, open(f"{HERE}/out/c7_timesweep_meta.json", "w"), indent=1)

TT = [t for t in TIMES if f"B1|{t}|enr" in store and f"B2|{t}|enr" in store]
iu = np.triu_indices(T, 1)
def sg(bl, t):
    o, x = store[f"{bl}|{t}|obs"], store[f"{bl}|{t}|exp"]
    return ((o - x) / np.sqrt(np.maximum(x, 1e-9)) >= Z_THR) & (o >= MIN_OBS)
SIG = {t: sg("B1", t) & sg("B2", t) for t in TT}
ENR = {t: (store[f"B1|{t}|enr"] + store[f"B2|{t}|enr"]) / 2 for t in TT}

print(f"\n{'t':>7}{'lineages':>11}{'clades>=3':>11}{'tips':>10}{'sig pairs':>11}")
for t in TT:
    m = [r for r in meta if r["blast"] == "B1" and r["t"] == t][0]
    print(f"{t:>7.2f}{m['lineages']:>11,}{m['clades']:>11,}{m['tips']:>10,}{int(SIG[t][iu].sum()):>11}")

# ---- divergence age per pair, and the A-vs-B reproducibility test ---------------------------
def deepest(sigfn):
    D = np.full((T, T), np.nan)
    for t in TT:
        S = sigfn(t)
        D[S & np.isnan(D)] = t          # TT ascending -> first hit is the earliest
    return D
Dboth = deepest(lambda t: SIG[t])
DA = deepest(lambda t: sg("B1", t)); DB = deepest(lambda t: sg("B2", t))
ok = ~np.isnan(DA) & ~np.isnan(DB)
pa, pb = DA[iu], DB[iu]; m = ~np.isnan(pa) & ~np.isnan(pb)
rho = spearmanr(pa[m], pb[m]).correlation; rr = pearsonr(pa[m], pb[m])[0]
med_abs = float(np.median(np.abs(pa[m] - pb[m])))
print(f"\n=== does it work? divergence age estimated independently in A and in B ===")
print(f"  pairs with an estimate in both halves: {int(m.sum())}")
print(f"  Spearman rho = {rho:.3f}   Pearson r = {rr:.3f}   median |A-B| = {med_abs:.2f} days")
within = [float(np.mean(np.abs(pa[m] - pb[m]) <= d)) for d in (0.25, 0.5, 1.0)]
print(f"  agreement within 0.25 d: {within[0]:.0%}   0.5 d: {within[1]:.0%}   1.0 d: {within[2]:.0%}")
# contiguity of significance in time (is the signal a clean band?)
pairs = [(a, b) for a, b in zip(*iu) if not np.isnan(Dboth[a, b])]
contig = 0
for a, b in pairs:
    s = np.array([SIG[t][a, b] for t in TT])
    idx = np.where(s)[0]
    if len(idx) and idx.max() - idx.min() + 1 == len(idx): contig += 1
print(f"  significance forms one contiguous time band for {contig}/{len(pairs)} pairs "
      f"({contig/max(len(pairs),1):.0%})")

gl = {}
for fn, okf in [(f"{CK}/germ_layer_map_validated.csv", lambda r: r["status"] in ("data-backed", "tree-resolved")),
                (f"{CK}/germ_layer_map_v6.csv", lambda r: r["confidence"] == "clear")]:
    for r in csv.DictReader(open(fn)):
        if r["celltype"] not in gl and r["germ_layer"] and okf(r): gl[r["celltype"]] = r["germ_layer"]
with open(f"{HERE}/out/c7_divergence_ages.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["celltype_1", "celltype_2", "germ_layer_1", "germ_layer_2", "same_germ_layer",
                "divergence_age_both", "divergence_age_A", "divergence_age_B",
                "n_slices_significant", "max_log2_enr"])
    for a, b in sorted(pairs, key=lambda p: Dboth[p[0], p[1]]):
        cur = [ENR[t][a, b] for t in TT]
        w.writerow([universe[a], universe[b], gl.get(universe[a], ""), gl.get(universe[b], ""),
                    int(gl.get(universe[a], "x") == gl.get(universe[b], "y")),
                    Dboth[a, b], DA[a, b] if not np.isnan(DA[a, b]) else "",
                    DB[a, b] if not np.isnan(DB[a, b]) else "",
                    int(sum(SIG[t][a, b] for t in TT)), round(float(max(cur)), 3)])
print(f"\n=== deepest couplings (earliest slice still significant in both halves) ===")
for a, b in sorted(pairs, key=lambda p: Dboth[p[0], p[1]])[:16]:
    print(f"  E{Dboth[a,b]:<6.2f} (A E{DA[a,b]:.2f} / B E{DB[a,b]:.2f})  {universe[a]} | {universe[b]}")

# ---- ordering + blocks reused from the K=15 analysis (Fig 6C), for comparability -------------
kc = np.load(f"{HERE}/out/c4_sweep.npz", allow_pickle=True)
ref = (kc["B1|15|enr"] + kc["B2|15|enr"]) / 2
Dm = 1 - np.corrcoef(ref); np.fill_diagonal(Dm, 0); Dm = (Dm + Dm.T) / 2
Z = linkage(squareform(Dm, checks=False), "ward"); order = leaves_list(Z)
cl = fcluster(Z, t=7, criterion="maxclust")
ordered = [universe[i] for i in order]
NAME = {"Microglia": "Haematopoiesis", "Thalamic neuronal precursors": "Anterior neuroectoderm",
        "Neural crest (PNS glia)": "Neural crest / PNS", "Spinal cord/r7/r8": "Posterior neuroectoderm",
        "Sclerotome": "Mesoderm / mesenchyme / endothelium",
        "Midgut/Hindgut epithelial cells": "Endoderm",
        "Olfactory epithelial cells": "Surface ectoderm / placodes"}
clof = {universe[i]: int(cl[i]) for i in range(T)}
bname = {clof[k]: v for k, v in NAME.items() if k in clof}
seq = [clof[t] for t in ordered]; bands = []; s = 0
for i in range(1, T + 1):
    if i == T or seq[i] != seq[i - 1]: bands.append((s, i - 1, bname.get(seq[s], "?"))); s = i

# ---- Fig: carpet over real developmental time -----------------------------------------------
curve = np.array([[ENR[t][a, b] for t in TT] for a, b in pairs])
sigm = np.array([[SIG[t][a, b] for t in TT] for a, b in pairs])
dv = np.array([Dboth[a, b] for a, b in pairs])
o1 = np.lexsort((-curve.max(1), dv))
lo, hi = np.percentile(curve, [2, 98])
norm = TwoSlopeNorm(vcenter=0.0, vmin=min(lo, -.1), vmax=max(hi, .1))
fig, ax = plt.subplots(figsize=(11.8, 12.2))
im = ax.imshow(curve[o1], cmap=DIV, norm=norm, aspect="auto", interpolation="nearest")
xt = [i for i, t in enumerate(TT) if abs(t * 2 - round(t * 2)) < 1e-6 and (round(t * 2) % 2 == 0)]
ax.set_xticks(xt); ax.set_xticklabels([f"E{TT[i]:.0f}" for i in xt])
ax.set_xlabel("developmental time of the clade ancestor  (earlier to the left)", labelpad=7)
ax.set_ylabel(f"cell-type pair (n = {len(pairs)}), ordered by divergence age")
lab = [f"{universe[pairs[i][0]]} | {universe[pairs[i][1]]}" for i in o1]
st = max(1, len(lab) // 58)
ax.set_yticks(range(0, len(lab), st)); ax.set_yticklabels([lab[i] for i in range(0, len(lab), st)], fontsize=4.7)
ys, xs = np.where(sigm[o1])
ax.plot(xs, ys, ls="", marker=".", ms=1.0, color="white", alpha=.5)
ax.set_title("Coupling across developmental time — clades defined by dated-ancestor age\n"
             "white dots = significant in both blastomeres", fontsize=10.5, pad=10)
for sp in ax.spines.values(): sp.set_visible(False)
fig.colorbar(im, ax=ax, shrink=.45, label="log₂ enrichment (mean of A, B)")
fig.tight_layout(); fig.savefig(f"{HERE}/out/c7_carpet.png", dpi=145, bbox_inches="tight")

# ---- Fig: divergence-age matrix -------------------------------------------------------------
P = Dboth[np.ix_(order, order)]
fig, ax = plt.subplots(figsize=(12.6, 11))
ax.set_facecolor("#f4f5f6")
im = ax.imshow(P, cmap=SEQ, aspect="equal", interpolation="nearest",
               norm=Normalize(vmin=np.nanmin(Dboth), vmax=np.nanmax(Dboth)))
for a, b, nm in bands:
    ax.add_patch(plt.Rectangle((a - .5, a - .5), b - a + 1, b - a + 1, fill=False, edgecolor=INK, lw=1.3))
    ax.text(b + 1.2, (a + b) / 2, nm, fontsize=7.4, va="center", ha="left", color=INK)
ax.set_xticks(range(T)); ax.set_xticklabels(ordered, rotation=90, fontsize=4.6)
ax.set_yticks(range(T)); ax.set_yticklabels(ordered, fontsize=4.6)
ax.set_xlim(-.5, T + 26); ax.tick_params(length=0)
ax.set_title("Divergence age of each cell-type coupling\n"
             "colour = earliest developmental time at which the pair is still significantly coupled; "
             "grey = never", fontsize=10.5, pad=10)
for sp in ax.spines.values(): sp.set_visible(False)
cb = fig.colorbar(im, ax=ax, shrink=.42, label="divergence age (E-day)"); cb.outline.set_visible(False)
fig.tight_layout(); fig.savefig(f"{HERE}/out/c7_divmatrix.png", dpi=145, bbox_inches="tight")

# ---- Fig: A vs B divergence age --------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.4, 6.2))
jit = (np.arange(int(m.sum())) % 5 - 2) * 0.022
ax.scatter(pa[m] + jit, pb[m] - jit, s=26, c="#2f6db3", alpha=.6, lw=0, zorder=3)
lim = [min(pa[m].min(), pb[m].min()) - .3, max(pa[m].max(), pb[m].max()) + .3]
ax.plot(lim, lim, ls="--", lw=.9, color=MUTED, zorder=1)
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("divergence age from blastomere A (E-day)")
ax.set_ylabel("divergence age from blastomere B (E-day)")
ax.set_title("Independent replication of the divergence-age estimate\n"
             f"Spearman ρ = {rho:.2f}, median |A−B| = {med_abs:.2f} d, "
             f"{within[1]:.0%} within 0.5 d  (n = {int(m.sum())})", fontsize=10)
ax.grid(alpha=.25, lw=.5); ax.set_axisbelow(True)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig(f"{HERE}/out/c7_ab_divage.png", dpi=150)
print("\nwrote c7_timesweep.npz, c7_divergence_ages.csv, c7_carpet.png, c7_divmatrix.png, c7_ab_divage.png")
