#!/usr/bin/env python3
"""STEP 4 -- emit the three shippable manuscript files and verify every numeric claim in the
subsection prose against the underlying data.

Shippable files (written to out/ship/), each CSV beside its visual:
  tableSX_postmitotic_celltypes.csv / .png            Table XX  -- the 27 post-mitotic types
  fig6B_captured_divisions.csv / _draft.png           Fig 6B    -- the 34 captured divisions
  figSX_captured_divisions_AB_replication.csv / _draft.png  Supp Fig -- all tested pairings, A vs B

Table XX is rendered as a table, not a chart: the quantities are identity and bookkeeping
(which types were evaluated, how many progenitors each received), and the counts span 0-6, so a
colour ramp would imply variation that isn't there. Zero counts are set in muted ink so the ten
types with no proposed progenitor read at a glance.
"""
import csv, os, shutil, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CK = "/Users/jay.shendure/Dropbox/claude/penultimate_clade_k_analysis"
CELL = "ageE12.5_disc1"
SHIP = f"{HERE}/out/ship"
os.makedirs(SHIP, exist_ok=True)

gl = {}
for fn, ok in [(f"{CK}/germ_layer_map_validated.csv", lambda r: r["status"] in ("data-backed", "tree-resolved")),
               (f"{CK}/germ_layer_map_v6.csv", lambda r: r["confidence"] == "clear")]:
    for r in csv.DictReader(open(fn)):
        if r["celltype"] not in gl and r["germ_layer"] and ok(r): gl[r["celltype"]] = r["germ_layer"]
GLNICE = {"BLOOD": "Blood", "MESODERM": "Mesoderm", "NEURAL_CREST": "Neural crest",
          "NEUROECTODERM": "Neuroectoderm", "SURFACE_ECTO": "Surface ectoderm", "ENDODERM": "Endoderm"}

post = [r for r in csv.DictReader(open(f"{CK}/celltype_postmitotic_v6.csv"))
        if r["classification"] == "post_mitotic"]
allp = list(csv.DictReader(open(f"{HERE}/out/h2_pairings_{CELL}.csv")))
cap = [r for r in allp if r["captured"] == "1"]
cap.sort(key=lambda r: -float(r["fold"]))

# ---- Table XX -----------------------------------------------------------------------------
# n_tips_in_tree is the operative denominator: the number of tips of that type in the
# 1,340,794-tip analysed tree. celltype_postmitotic_v6.csv's `n` counts cells in the
# transcriptome metadata instead, which is larger for 25 of the 27 types (up to 6x for
# primitive erythroid cells), so both are reported and the tree count is used for the visual.
import collections
_d = np.load(f"{HERE}/out/h1_pairs.npz", allow_pickle=True)
_types = list(_d["types"]); _tipn = collections.Counter(_d["leaf_ct"].astype(int).tolist())
ntip = lambda ct: _tipn.get(_types.index(ct), 0) if ct in _types else 0

TBROWS = sorted(post, key=lambda r: (GLNICE.get(gl.get(r["celltype"], "?"), "zz"), r["celltype"]))
with open(f"{SHIP}/tableSX_postmitotic_celltypes.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["celltype", "major_trajectory", "germ_layer", "n_tips_in_tree",
                "n_cells_metadata", "progenitors_proposed", "pairings_tested"])
    for r in TBROWS:
        ct = r["celltype"]
        w.writerow([ct, r["trajectory"], GLNICE.get(gl.get(ct, ""), ""), ntip(ct), r["n"],
                    sum(1 for c in cap if c["postmitotic"] == ct),
                    sum(1 for p in allp if p["postmitotic"] == ct)])

# ---- Fig 6B source data -------------------------------------------------------------------
F6 = ["progenitor", "postmitotic", "prog_germlayer", "pm_germlayer", "same_germlayer",
      "n_cells", "n_expected", "fold", "q", "n_cells_A", "fold_A", "q_A", "n_cells_B", "fold_B", "q_B"]
with open(f"{SHIP}/fig6B_captured_divisions.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(F6)
    for r in cap:
        w.writerow([r["progenitor"], r["postmitotic"], GLNICE.get(r["prog_germlayer"], ""),
                    GLNICE.get(r["pm_germlayer"], ""), r["same_germlayer"], r["obs"], r["exp"],
                    r["fold"], r["q"], r["obsA"], r["foldA"], r["qA"], r["obsB"], r["foldB"], r["qB"]])

# ---- Supp Fig source data -----------------------------------------------------------------
with open(f"{SHIP}/figSX_captured_divisions_AB_replication.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["progenitor", "postmitotic", "n_cells_A", "fold_A", "q_A",
                "n_cells_B", "fold_B", "q_B", "tested_in_both", "captured"])
    for r in sorted(allp, key=lambda r: -float(r["fold"])):
        tb = int(bool(r["foldA"]) and bool(r["foldB"]))
        w.writerow([r["progenitor"], r["postmitotic"], r["obsA"], r["foldA"], r["qA"],
                    r["obsB"], r["foldB"], r["qB"], tb, r["captured"]])

# ---- visuals, placed beside the CSVs ------------------------------------------------------
shutil.copyfile(f"{HERE}/out/h3_heatmap.png",    f"{SHIP}/fig6B_captured_divisions_draft.png")
shutil.copyfile(f"{HERE}/out/h3_ab_scatter.png", f"{SHIP}/figSX_captured_divisions_AB_replication_draft.png")

INK, INK2, MUTED, RULE = "#1f2328", "#4a5158", "#a0a6ac", "#d7dbdf"
tb = list(csv.DictReader(open(f"{SHIP}/tableSX_postmitotic_celltypes.csv")))
groups = []
for r in tb:
    if not groups or groups[-1][0] != r["germ_layer"]: groups.append((r["germ_layer"], []))
    groups[-1][1].append(r)
nrow = len(tb) + len(groups)
fig, ax = plt.subplots(figsize=(11.6, 0.265 * nrow + 1.9))
ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(-3.0, nrow + 3.4)
XCT, XTR, XTIP, XMETA, XPRO, XTEST = 0.022, 0.360, 0.640, 0.760, 0.884, 0.988
HEAD = [(XCT, "Post-mitotic cell type", "left"), (XTR, "Major trajectory", "left"),
        (XTIP, "Tips in\ntree", "right"), (XMETA, "Cells in\nmetadata", "right"),
        (XPRO, "Progenitors\nproposed", "right"), (XTEST, "Pairings\ntested", "right")]
y = nrow + 1.55
for x, lab, ha in HEAD:
    ax.text(x, y, lab, ha=ha, va="center", fontsize=8.3, color=INK2, weight="semibold",
            linespacing=1.45)
y -= 1.05
ax.plot([0.012, 0.995], [y, y], color=INK2, lw=1.0)
for glname, rows in groups:
    y -= 1.05
    ax.text(XCT, y, glname.upper(), ha="left", va="center", fontsize=8.1, color=INK, weight="bold")
    for r in rows:
        y -= 1.0
        npro, ntest = int(r["progenitors_proposed"]), int(r["pairings_tested"])
        ax.text(XCT + 0.016, y, r["celltype"], ha="left", va="center", fontsize=8.3, color=INK)
        ax.text(XTR, y, r["major_trajectory"].replace("_", " "), ha="left", va="center",
                fontsize=7.6, color=MUTED)
        ax.text(XTIP, y, f"{int(r['n_tips_in_tree']):,}", ha="right", va="center", fontsize=8.3, color=INK2)
        ax.text(XMETA, y, f"{int(r['n_cells_metadata']):,}", ha="right", va="center", fontsize=7.8, color=MUTED)
        ax.text(XPRO, y, str(npro), ha="right", va="center", fontsize=8.9,
                color=INK if npro else MUTED, weight="semibold" if npro else "normal")
        ax.text(XTEST, y, str(ntest), ha="right", va="center", fontsize=8.3, color=INK2 if ntest else MUTED)
    y -= 0.55
    ax.plot([0.012, 0.995], [y, y], color=RULE, lw=0.7)
nprop = sum(1 for r in tb if int(r["progenitors_proposed"]) > 0)
ncap = sum(int(r["progenitors_proposed"]) for r in tb)
ax.text(0.012, nrow + 3.0, "Table XX.  Post-mitotic cell types evaluated for "
        "differentiation-associated tree-sibling pairings", ha="left", va="center",
        fontsize=10.5, color=INK, weight="semibold")
ax.text(0.012, y - 1.05, f"{len(tb)} post-mitotic types evaluated; {nprop} received at least one "
        f"proposed progenitor ({ncap} captured divisions in total).", ha="left", va="center",
        fontsize=7.9, color=INK2)
ax.text(0.012, y - 1.95, "Tips in tree = cells of that type among the 1,340,794 tips analysed "
        "(the operative denominator); cells in metadata = all cells of that type in the "
        "transcriptome. Pairings tested =\nprogenitor→post-mitotic pairings seen ≥3 times among "
        "qualifying tree siblings (shared ancestor ≥E12.5, ≤1 discordant site).",
        ha="left", va="top", fontsize=7.3, color=MUTED, linespacing=1.5)
fig.savefig(f"{SHIP}/tableSX_postmitotic_celltypes.png", dpi=200, bbox_inches="tight",
            facecolor="white", pad_inches=0.28)
plt.close(fig)

# ---- verification of every numeric claim in the prose -------------------------------------
d = np.load(f"{HERE}/out/h1_pairs.npz", allow_pickle=True)
sw = {(r["age_min"], r["disc_max"]): r for r in csv.DictReader(open(f"{HERE}/out/h2_sweep.csv"))}
s = sw[("E12.5", "1")]
folds = [float(r["fold"]) for r in cap]
tb = [r for r in allp if r["foldA"] and r["foldB"]]
from scipy.stats import spearmanr
rho = spearmanr([float(r["foldA"]) for r in tb], [float(r["foldB"]) for r in tb]).correlation
ncg = [r for r in cap if r["prog_germlayer"] == "NEURAL_CREST"]
pns = [r for r in ncg if r["pm_germlayer"] == "NEURAL_CREST"]
mk = [r for r in cap if r["postmitotic"] == "Megakaryocytes"]
F = lambda a, b: next(float(r["fold"]) for r in cap if r["progenitor"] == a and r["postmitotic"] == b)

# v8 values (regenerated 2026-08-25 on the RT-cross-talk-corrected trees).
# These encode the prose numbers so the shipped data and the manuscript cannot drift apart.
Fq = lambda a, b: next((float(r["fold"]) for r in cap
                        if r["progenitor"] == a and r["postmitotic"] == b), None)
checks = [
    ("870,785 tree-sibling pairs",            len(d["i"]),                          870785),
    ("27 post-mitotic cell types",            len(post),                            27),
    ("82,300 post-mitotic cells",             int(s["cells_paired"]),               82300),
    ("640 pairings tested",                   int(s["tested_pool"]),                640),
    ("37 pairings captured",                  len(cap),                             37),
    ("max fold 29.2",                         round(max(folds), 1),                 29.2),
    ("min fold >= 3",                         min(folds) >= 3.0,                    True),
    ("13,549 post-mitotic cells captured",    sum(int(r["obs"]) for r in cap),      13549),
    ("Spearman rho 0.82",                     round(rho, 2),                        0.82),
    ("18 of 27 types w/ >=1 progenitor",      len({r["postmitotic"] for r in cap}), 18),
    ("24 distinct progenitor states",         len({r["progenitor"] for r in cap}),  24),
    ("35 of 37 within one germ layer",        sum(1 for r in cap if r["same_germlayer"] == "1"), 35),
]
print("=== captured divisions (v8), by fold ===")
for r in sorted(cap, key=lambda r: -float(r["fold"])):
    print(f"  {float(r['fold']):6.1f}x  {int(r['obs']):>5} cells  "
          f"{r['progenitor']} -> {r['postmitotic']}")
print()
print("=== verification of prose numbers ===")
bad = 0
for name, got, want in checks:
    ok = got == want
    bad += (not ok)
    print(f"  [{'OK ' if ok else 'FAIL'}] {name:<40} got={got!r:<10} expected={want!r}")
print(f"\n{len(checks)-bad}/{len(checks)} checks pass")
print(f"\nother claims in the prose (no single expected value):")
print(f"  osteoclasts from {sum(1 for r in cap if r['postmitotic']=='Osteoclasts')} progenitors; "
      f"primitive erythroid from {sum(1 for r in cap if r['postmitotic']=='Primitive erythroid cells')}; "
      f"granulocytes from {sum(1 for r in cap if r['postmitotic']=='Granulocytes')}")
print(f"  neural-crest-progenitor pairings total {len(ncg)} (of which {len(pns)} to peripheral neurons, "
      f"{len(ncg)-len(pns)} to CNS = cranial motor neurons)")
print(f"  pairings tested in both halves: {len(tb)}")
print(f"\nwrote {SHIP}/ (3 files)")
