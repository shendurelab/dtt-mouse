#!/usr/bin/env python3
"""Table S5 (v8) -- the 27 post-mitotic cell types, now naming the progenitors asserted for
each, and repeating the analysis columns on the backbone tree.

Replaces the v6 contents of the shipped sheet. What changed and why:

  * v6 -> v8. Every number is re-derived from the fixed final placement tree
    (mergedtree_dttpq_v8.npz, 1,281,141 tips) and the v8 metadata. Some heterotypic
    pairings did shift; the per-pairing sheet makes that inspectable.
  * COLUMN H names the progenitors. Previously the table gave only the count
    (`progenitors_proposed`); it now lists each asserted progenitor with its pooled
    fold-enrichment and BH q, strongest first.
  * COLUMNS I-L repeat D/F/G plus the progenitor list on the BACKBONE tree
    (mergedtree_bbfinal_v8.npz, 655,701 tips), the placement-free control. `n_cells_metadata`
    is deliberately NOT repeated: it counts cells in the transcriptome metadata and is
    identical for both trees.
  * COLUMNS M-N carry the testability bookkeeping behind the manuscript's
    "24 of 31 heterotypic pairings testable on both remained significant": per type, how many
    of its full-tree captured pairings are even testable on the backbone (>= MINOBS observed
    there), and how many of those stay significant.

Definitions are h2r.py's, at the manuscript's operative grid cell (MRCA age >= E12.5, at most
1 discordant monomer; global within-blastomere null, NPERM=500). A pairing is CAPTURED when
pooled q < 0.01 and pooled fold >= 3.0. A pairing is TESTED when it clears the obs >= 3 floor
that BH is applied over; pairings below it are neither tested nor counted.

Reads (all pre-computed, nothing re-run here):
  figures/v8/out/h1_pairs.npz  /  h1_pairs_bb.npz             tip counts per cell type
  figures/v8/out/h2_pairings_ageE12.5_disc1.csv                full tree,  global null
  figures/v8/out/h2_pairings_bb_global_ageE12.5_disc1.csv      backbone,   global null
  fig6_v8/celltype_postmitotic_v8.csv                          the 27 types + germ layer
  final_push/data/cell_metadata.v8.txt                         trajectory + metadata counts

Writes:
  tableS5_postmitotic_v8.csv     the table itself (27 rows)
  tableS5b_pairings_v8.csv       every tested pairing, both trees side by side
"""
import collections
import csv
import os

import numpy as np

V8 = "/Users/jay.shendure/Dropbox/claude/mouse_sprint/tape_pipeline/figures/v8/out/"
FP = "/Users/jay.shendure/Dropbox/claude/current/final_push/"
CELL = "ageE12.5_disc1"
MINOBS = 3          # h2r.py's BH floor; below it a pairing is not tested at all
GLNICE = {"BLOOD": "Blood", "MESODERM": "Mesoderm", "NEURAL_CREST": "Neural crest",
          "NEUROECTODERM": "Neuroectoderm", "SURFACE_ECTO": "Surface ectoderm",
          "ENDODERM": "Endoderm"}


def tipcounts(npz):
    d = np.load(V8 + npz, allow_pickle=True)
    types = list(d["types"])
    n = collections.Counter(d["leaf_ct"].astype(int).tolist())
    return {t: n.get(i, 0) for i, t in enumerate(types)}


def pairings(fn):
    return list(csv.DictReader(open(V8 + fn)))


def fmt_prog(rows):
    """'Progenitor (12.3x, q=4e-16); ...', strongest first."""
    rows = sorted(rows, key=lambda r: -float(r["fold"]))
    out = []
    for r in rows:
        q = float(r["q"])
        qs = "q<1e-300" if q == 0 else f"q={q:.1e}"
        out.append(f'{r["progenitor"]} ({float(r["fold"]):.1f}x, {qs})')
    return "; ".join(out)


# ---------------------------------------------------------------- inputs
post = [r for r in csv.DictReader(open(FP + "fig6_v8/celltype_postmitotic_v8.csv"))
        if r["classification"] == "post_mitotic"]

traj, ncells = {}, collections.Counter()
with open(FP + "data/cell_metadata.v8.txt") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        traj.setdefault(row["celltype"], row["major_trajectory"])
        ncells[row["celltype"]] += 1

tip_full, tip_bb = tipcounts("h1_pairs.npz"), tipcounts("h1_pairs_bb.npz")
P_full = pairings(f"h2_pairings_{CELL}.csv")
P_bb = pairings(f"h2_pairings_bb_global_{CELL}.csv")

cap_full = [r for r in P_full if r["captured"] == "1"]
cap_bb = [r for r in P_bb if r["captured"] == "1"]
key = lambda r: (r["progenitor"], r["postmitotic"])
bb_by_key = {key(r): r for r in P_bb}
bb_cap_keys = {key(r) for r in cap_bb}

by_pm = lambda rows: collections.defaultdict(list, {
    k: list(v) for k, v in collections.groupby(sorted(rows, key=lambda r: r["postmitotic"]),
                                               key=lambda r: r["postmitotic"])})
from itertools import groupby as _gb
def group(rows):
    d = collections.defaultdict(list)
    for r in rows:
        d[r["postmitotic"]].append(r)
    return d

g_cap_full, g_all_full = group(cap_full), group(P_full)
g_cap_bb, g_all_bb = group(cap_bb), group(P_bb)

# ---------------------------------------------------------------- the table
HDR = ["celltype", "major_trajectory", "germ_layer",
       "n_tips_in_tree", "n_cells_metadata", "progenitors_proposed", "pairings_tested",
       "progenitors_asserted",
       "n_tips_in_backbone", "progenitors_proposed_backbone", "pairings_tested_backbone",
       "progenitors_asserted_backbone",
       "pairings_testable_on_both", "pairings_retained_on_backbone"]

rows_out = []
tot_testable = tot_retained = 0
for r in sorted(post, key=lambda r: (GLNICE.get(r["germ_layer"], "zz"), r["celltype"])):
    ct = r["celltype"]
    cf, af = g_cap_full.get(ct, []), g_all_full.get(ct, [])
    cb, ab = g_cap_bb.get(ct, []), g_all_bb.get(ct, [])
    testable = [x for x in cf if key(x) in bb_by_key]
    retained = [x for x in testable if key(x) in bb_cap_keys]
    tot_testable += len(testable); tot_retained += len(retained)
    rows_out.append([ct, traj.get(ct, ""), GLNICE.get(r["germ_layer"], r["germ_layer"]),
                     tip_full.get(ct, 0), ncells[ct], len(cf), len(af), fmt_prog(cf),
                     tip_bb.get(ct, 0), len(cb), len(ab), fmt_prog(cb),
                     len(testable), len(retained)])

with open(FP + "tableS5_postmitotic_v8.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(HDR); w.writerows(rows_out)

# ---------------------------------------------------------------- per-pairing sheet
HDR2 = ["progenitor", "postmitotic", "prog_germlayer", "pm_germlayer", "same_germlayer",
        "obs", "n_expected", "fold", "q", "captured",
        "obs_A", "fold_A", "q_A", "obs_B", "fold_B", "q_B", "concordant_AB",
        "obs_backbone", "n_expected_backbone", "fold_backbone", "q_backbone",
        "captured_backbone", "backbone_status"]
rows2 = []
for r in sorted(P_full, key=lambda r: (-int(r["captured"]), -float(r["fold"]))):
    b = bb_by_key.get(key(r))
    if b is None:
        status = "not testable on backbone (<%d cells)" % MINOBS
        bcols = ["", "", "", "", ""]
    else:
        status = "significant" if b["captured"] == "1" else "tested, not significant"
        bcols = [b["obs"], b["exp"], b["fold"], b["q"], b["captured"]]
    rows2.append([r["progenitor"], r["postmitotic"], GLNICE.get(r["prog_germlayer"], ""),
                  GLNICE.get(r["pm_germlayer"], ""), r["same_germlayer"],
                  r["obs"], r["exp"], r["fold"], r["q"], r["captured"],
                  r["obsA"], r["foldA"], r["qA"], r["obsB"], r["foldB"], r["qB"],
                  r["concordant_AB"]] + bcols + [status])
with open(FP + "tableS5b_pairings_v8.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(HDR2); w.writerows(rows2)

# ---------------------------------------------------------------- checks against the prose
print("full placed tree (v8), MRCA >= E12.5, <=1 discordant monomer, global null")
print(f"  pairings tested              {len(P_full)}")
print(f"  captured (q<0.01, fold>=3)   {len(cap_full)}")
print(f"  max fold                     {max(float(r['fold']) for r in cap_full):.1f}x")
print(f"  distinct progenitor states   {len({r['progenitor'] for r in cap_full})}")
print(f"  post-mitotic types assigned  {len({r['postmitotic'] for r in cap_full})} / {len(post)}")
print(f"  post-mitotic cells in tested pairings  {sum(int(r['obs']) for r in P_full):,}")
print("\nbackbone tree (v8), same grid cell and null")
print(f"  pairings tested              {len(P_bb)}")
print(f"  captured                     {len(cap_bb)}")
print(f"  post-mitotic types assigned  {len({r['postmitotic'] for r in cap_bb})}")
print(f"\nof the {len(cap_full)} full-tree captured pairings:")
print(f"  testable on the backbone     {tot_testable}")
print(f"  still significant there      {tot_retained}")
print(f"  lost to power (<{MINOBS} cells)      {len(cap_full)-tot_testable}")
print(f"  tested but not significant   {tot_testable-tot_retained}")
print(f"\nwrote {FP}tableS5_postmitotic_v8.csv ({len(rows_out)} rows)")
print(f"wrote {FP}tableS5b_pairings_v8.csv ({len(rows2)} rows)")
