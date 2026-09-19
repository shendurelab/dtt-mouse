#!/usr/bin/env python3
"""STEP 10 -- package the three figures for the coupling-depth paragraph: draft plots plus the
underlying values, named to the manuscript panels.

  figXa_coupling_carpet          Fig. Xa   enrichment of every coupled pair at every time slice
  figXb_coupling_hierarchy       Fig. Xb   the dated hierarchy (average linkage over coupling depth)
  figSXa_coupling_depth_AB       Supp Xa   coupling depth estimated independently in A and in B

Writes out/ship_coupling/. Panel letters are placeholders (X) matching the draft text.
"""
import numpy as np, csv, os, shutil, json
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr, pearsonr

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CK = "/Users/jay.shendure/Dropbox/claude/penultimate_clade_k_analysis"
SHIP = f"{HERE}/out/ship_coupling"
CENSOR, Z_THR, MIN_OBS = 13.5, 3.0, 5
os.makedirs(SHIP, exist_ok=True)

d = np.load(f"{HERE}/out/c7_timesweep.npz", allow_pickle=True)
types = list(d["types"]); T = len(types)
TT = [float(t) for t in d["TIMES"] if f"B1|{float(t)}|enr" in d and f"B2|{float(t)}|enr" in d]
meta = json.load(open(f"{HERE}/out/c7_timesweep_meta.json"))
def sg(bl, t):
    o, x = d[f"{bl}|{t}|obs"], d[f"{bl}|{t}|exp"]
    return ((o - x) / np.sqrt(np.maximum(x, 1e-9)) >= Z_THR) & (o >= MIN_OBS)
SIG = {t: sg("B1", t) & sg("B2", t) for t in TT}
iu = np.triu_indices(T, 1)
def depth_matrix(fn):
    D = np.full((T, T), np.nan)
    for t in TT:
        S = fn(t); D[S & np.isnan(D)] = t
    return D
Dboth = depth_matrix(lambda t: SIG[t])
DA = depth_matrix(lambda t: sg("B1", t)); DB = depth_matrix(lambda t: sg("B2", t))
pairs = [(a, b) for a, b in zip(*iu) if not np.isnan(Dboth[a, b])]

gl = {}
for fn, okf in [(f"{CK}/germ_layer_map_validated.csv", lambda r: r["status"] in ("data-backed", "tree-resolved")),
                (f"{CK}/germ_layer_map_v6.csv", lambda r: r["confidence"] == "clear")]:
    for r in csv.DictReader(open(fn)):
        if r["celltype"] not in gl and r["germ_layer"] and okf(r): gl[r["celltype"]] = r["germ_layer"]
NICE = {"BLOOD": "Blood", "MESODERM": "Mesoderm", "NEURAL_CREST": "Neural crest",
        "NEUROECTODERM": "Neuroectoderm", "SURFACE_ECTO": "Surface ectoderm", "ENDODERM": "Endoderm"}
G = lambda c: NICE.get(gl.get(c, ""), "")

# Order each pair so the cell type appearing EARLIER in the Qiu et al. developmental graph is
# listed first (see c11_qiu_depth.py). Guidepost only; affects no statistic. Ties/unmapped ->
# alphabetical, so the ordering is deterministic.
qd = {}
for r in csv.DictReader(open(f"{HERE}/out/qiu_celltype_depth.csv")):
    if r["qiu_depth_from_oocyte"] != "": qd[r["celltype"]] = int(r["qiu_depth_from_oocyte"])
def ordpair(x, y):
    kx = (qd.get(x, 10**6), x); ky = (qd.get(y, 10**6), y)
    return (x, y) if kx <= ky else (y, x)

# ---- Fig Xa source data: long format, one row per (pair, time slice) -------------------------
order = sorted(pairs, key=lambda p: (Dboth[p[0], p[1]],
                                     -max(((d[f"B1|{t}|enr"][p] + d[f"B2|{t}|enr"][p]) / 2) for t in TT)))
with open(f"{SHIP}/figXa_coupling_carpet.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["row_in_figure", "celltype_1", "celltype_2", "germ_layer_1", "germ_layer_2",
                "clade_ancestor_age_E", "log2_enrichment_mean", "log2_enrichment_A",
                "log2_enrichment_B", "significant_both_blastomeres", "coupling_depth_E"])
    for n, (a, b) in enumerate(order, start=1):
        for t in TT:
            eA, eB = float(d[f"B1|{t}|enr"][a, b]), float(d[f"B2|{t}|enr"][a, b])
            c1, c2 = ordpair(types[a], types[b])
            w.writerow([n, c1, c2, G(c1), G(c2), t,
                        round((eA + eB) / 2, 4), round(eA, 4), round(eB, 4),
                        int(SIG[t][a, b]), Dboth[a, b]])

# ---- Fig Xb source data: the dated merges ---------------------------------------------------
keep = sorted({a for a, b in pairs} | {b for a, b in pairs})
kn = [types[i] for i in keep]
M = Dboth[np.ix_(keep, keep)].copy(); M[np.isnan(M)] = CENSOR
np.fill_diagonal(M, 0.0); M = (M + M.T) / 2
Z = linkage(squareform(M, checks=False), method="average")
n = len(kn); lab = {i: [kn[i]] for i in range(n)}
with open(f"{SHIP}/figXb_coupling_hierarchy.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["merge_order", "branch_point_age_E", "n_cell_types", "members"])
    for i, (x, y, h, cnt) in enumerate(Z):
        mem = lab[int(x)] + lab[int(y)]; lab[n + i] = mem
        w.writerow([i + 1, round(float(h), 3), int(cnt), "; ".join(sorted(mem))])
with open(f"{SHIP}/figXb_coupling_depth_matrix.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow([""] + kn)
    for i, r in enumerate(Dboth[np.ix_(keep, keep)]):
        w.writerow([kn[i]] + ["" if np.isnan(v) else round(float(v), 2) for v in r])

# ---- Supp Fig Xa source data ----------------------------------------------------------------
pa, pb = DA[iu], DB[iu]; m = ~np.isnan(pa) & ~np.isnan(pb)
rho = spearmanr(pa[m], pb[m]).correlation; rp = pearsonr(pa[m], pb[m])[0]
med = float(np.median(np.abs(pa[m] - pb[m])))
w05 = float(np.mean(np.abs(pa[m] - pb[m]) <= 0.5))
with open(f"{SHIP}/figSXa_coupling_depth_AB.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["celltype_1", "celltype_2", "germ_layer_1", "germ_layer_2",
                "coupling_depth_A_E", "coupling_depth_B_E", "abs_difference_days",
                "coupling_depth_both_E"])
    for a, b in zip(*iu):
        if np.isnan(DA[a, b]) or np.isnan(DB[a, b]): continue
        c1, c2 = ordpair(types[a], types[b])
        w.writerow([c1, c2, G(c1), G(c2), DA[a, b], DB[a, b],
                    round(abs(DA[a, b] - DB[a, b]), 2),
                    "" if np.isnan(Dboth[a, b]) else Dboth[a, b]])

# ---- the per-slice clade counts that define the time axis ------------------------------------
with open(f"{SHIP}/coupling_time_slices.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["clade_ancestor_age_E", "blastomere", "lineages_present",
                "clades_with_3plus_cells", "cells_in_those_clades", "significant_pairs_both"])
    for r in meta:
        ns = int(SIG[r["t"]][iu].sum()) if r["t"] in SIG else ""
        w.writerow([r["t"], "A" if r["blast"] == "B1" else "B", r["lineages"], r["clades"],
                    r["tips"], ns])

for src, dst in [("c8_carpet.png", "figXa_coupling_carpet_draft.png"),
                 ("c8_dendrogram.png", "figXb_coupling_hierarchy_draft.png"),
                 ("c7_ab_divage.png", "figSXa_coupling_depth_AB_draft.png")]:
    shutil.copyfile(f"{HERE}/out/{src}", f"{SHIP}/{dst}")

print(f"pairs {len(pairs)}; cell types in tree {len(kn)}; slices {len(TT)}")
print(f"A/B: Spearman rho {rho:.3f}, Pearson r {rp:.3f}, median |A-B| {med:.2f} d, "
      f"within 0.5 d {w05:.0%} (n={int(m.sum())})")
print("\n" + "\n".join(sorted(os.listdir(SHIP))))
