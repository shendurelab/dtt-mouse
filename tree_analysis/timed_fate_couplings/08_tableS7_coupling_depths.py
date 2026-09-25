import os as _os
_REPO = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))

def _results(*parts):
    """Intermediate/output dir for this analysis step, override with DTT_RESULTS."""
    base = _os.environ.get("DTT_RESULTS", _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "out"))
    p = _os.path.join(base, *parts) if parts else base
    _os.makedirs(_os.path.dirname(p) if _os.path.splitext(p)[1] else p, exist_ok=True)
    return p
#!/usr/bin/env python3
"""Table S6 (v8) -- timed cell-type couplings as a LIST, on the placed tree and the backbone,
under both the global tip-label null and the within-clone null.

Replaces the v6 "coupling-depth matrix across 76 cell types" (a 76x76 sparse grid holding 256
filled cells). A matrix that is ~4% full costs a page to say what a list says in 138 rows, and
it has no room for the per-blastomere dates, the evidence behind each call, a second tree, or a
second null.

CONTENTS. One row per coupled cell-type PAIR: the union of pairs coupled on the full placed
tree (129) and on the placement-free backbone tree (84), which overlap in 75 -- so 138 rows,
9 of them backbone-only. Ordered by coupling date, earliest first, so the table reads as a
developmental sequence rather than an alphabetical grid.

DEFINITIONS (unchanged from the c7 pipeline that produced the published numbers). A pair is
COUPLED when its co-occurrence within maximal terminal-proximate clades is significantly
enriched; the COUPLING DEPTH is the earliest time slice at which that enrichment holds, dated
in embryonic days. Depths are reported pooled and separately for the two independently
reconstructed half-embryos (blastomeres A and B). `n_slices_significant` and `max_log2_enr`
are the evidence behind the call.

THE TWO NULLS. The primary dates use the global null: tip labels permuted across all annotated
tips of the same blastomere. The within-clone null instead permutes labels only among cells
descended from the same E7.0 founder lineage (3,848 clones, median size 32), which is the
stricter test -- a pair that survives it is not explained by clonal composition alone. Each
tree therefore carries a clone-null depth and an outcome: `coupled` if the pair is still
coupled under that null, `not recovered` if it is not. The cell is blank where that tree did
not couple the pair under the global null in the first place, so there was nothing to retest.

TIES TO THE MANUSCRIPT, all recomputed and printed by this script rather than asserted:
  placed vs backbone, global null, over the 75 shared pairs
      median shift -0.75 d, Spearman rho = 0.838   ("deeper ... rank order rho = 0.84")
  backbone, global vs within-clone null
      55 of 84 survive, median +0.50 d, rho = 0.757  ("55/84 ... +0.50 days ... rho = 0.76")
  placed tree, global vs within-clone null (not quoted in the text)
      79 of 129 survive, median +0.50 d, rho = 0.822

Reads:  figures/v8/out/c7_divergence_ages.csv      placed tree, global null
        figures/v8/out/c7_divergence_ages_bb.csv   backbone tree, global null
        fig6_v8/coupling_full_cloneE7.csv          placed tree, within-clone null
        fig6_v8/coupling_bb_cloneE7.csv            backbone tree, within-clone null
Writes: tableS6_couplings_v8.csv        the list
        tableS6_couplings_v8_paste.csv  the same with the sheet's title block
"""
import csv
import statistics as st

from scipy.stats import spearmanr

V8 = _results() + _os.sep
F6 = _results() + _os.sep
FP = _results() + _os.sep
GLNICE = {"BLOOD": "Blood", "MESODERM": "Mesoderm", "NEURAL_CREST": "Neural crest",
          "NEUROECTODERM": "Neuroectoderm", "SURFACE_ECTO": "Surface ectoderm",
          "ENDODERM": "Endoderm"}

norm = lambda a, b: (a, b) if a <= b else (b, a)


def load_global(fn):
    out = {}
    for r in csv.DictReader(open(V8 + fn)):
        a, b = r["celltype_1"], r["celltype_2"]
        if a > b:
            # the pair is unordered, so normalise it and carry the germ layers with the names.
            # A and B are BLASTOMERES, not cell types -- those columns must not be swapped.
            a, b = b, a
            r["germ_layer_1"], r["germ_layer_2"] = r["germ_layer_2"], r["germ_layer_1"]
        r["celltype_1"], r["celltype_2"] = a, b
        out[(a, b)] = r
    return out


def load_clone(fn):
    return {norm(r["celltype_1"], r["celltype_2"]): float(r["coupling_depth_E"])
            for r in csv.DictReader(open(F6 + fn))}


F = load_global("c7_divergence_ages.csv")
B = load_global("c7_divergence_ages_bb.csv")
# the matched robustness series: one script (coupling_robust.py), nperm=50, global and
# within-clone nulls run identically on each tree, so the two can be compared to each other
FG = load_clone("coupling_valid_global50.csv")
BG = load_clone("coupling_backbone_global50.csv")
FC = load_clone("coupling_full_cloneE7.csv")
BC = load_clone("coupling_bb_cloneE7.csv")

num = lambda r, k: (float(r[k]) if r and r[k] not in ("", "nan", None) else None)
fmt = lambda v: "" if v is None else f"{v:g}"


def block(g, glob, clone, key):
    """Per tree: the five primary-dating columns, then the three matched-null columns."""
    gd, cd = glob.get(key), clone.get(key)
    if not g:
        # the pair is not coupled on this tree under the primary dating, so nothing to retest
        return ["", "", "", "", "", fmt(gd), fmt(cd), ""]
    return [fmt(num(g, "divergence_age_both")),
            fmt(num(g, "divergence_age_A")), fmt(num(g, "divergence_age_B")),
            g["n_slices_significant"], g["max_log2_enr"],
            fmt(gd), fmt(cd), "coupled" if cd is not None else "not recovered"]


rows = []
for k in sorted(set(F) | set(B)):
    f, b = F.get(k), B.get(k)
    src = f or b
    df, db = num(f, "divergence_age_both"), num(b, "divergence_age_both")
    where = "both trees" if (f and b) else ("placed tree only" if f else "backbone tree only")
    shift = "" if (df is None or db is None) else f"{db - df:+g}"
    rows.append({"sort": df if df is not None else db,
                 "row": [src["celltype_1"], src["celltype_2"],
                         GLNICE.get(src["germ_layer_1"], src["germ_layer_1"]),
                         GLNICE.get(src["germ_layer_2"], src["germ_layer_2"]),
                         src["same_germ_layer"], where]
                        + block(f, FG, FC, k) + block(b, BG, BC, k) + [shift]})

# earliest coupling first; ties broken by the pair name so the order is reproducible
rows.sort(key=lambda x: (x["sort"], x["row"][0], x["row"][1]))
body = [x["row"] for x in rows]

TREEBLOCK = ["coupling depth (E)", "depth in blastomere A (E)", "depth in blastomere B (E)",
             "# significant time slices", "max log2 enrichment",
             "global null: coupling depth (E)", "within-clone null: coupling depth (E)",
             "within-clone null: outcome"]
SUBBANNER = ["primary dating (global null)", "", "", "", "",
             "matched-null robustness series", "", ""]
HDR = (["cell type 1", "cell type 2", "germ layer 1", "germ layer 2", "same germ layer",
        "detected in"] + TREEBLOCK + TREEBLOCK + ["depth shift (backbone - placed, d)"])
BANNER = ([""] * 6 + ["placed tree"] + [""] * 7 + ["backbone tree"] + [""] * 7 + [""])
BANNER2 = [""] * 6 + SUBBANNER + SUBBANNER + [""]
assert len(HDR) == len(BANNER) == len(BANNER2) == len(body[0]), \
    (len(HDR), len(BANNER), len(BANNER2), len(body[0]))

with open(FP + "tableS6_couplings_v8.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(BANNER); w.writerow(BANNER2); w.writerow(HDR); w.writerows(body)

TITLE = ("Table S6. Timed cell-type couplings on the placed tree and the backbone tree, "
         "ordered by coupling date.")
SUB = ("One row per coupled cell-type pair: the union of pairs coupled on the full placed tree "
       "(v8, 129 pairs) and on the placement-free backbone tree (v8, 84 pairs), which overlap "
       "in 75. Coupling depth is the earliest time slice at which co-occurrence enrichment "
       "within maximal terminal-proximate clades holds, in embryonic days; A and B are the two "
       "independently reconstructed half-embryos. PRIMARY DATING gives the published depth for "
       "each tree. The MATCHED-NULL ROBUSTNESS SERIES re-dates the same pairs with the global "
       "tip-label null and with the stricter within-clone null -- which permutes labels only "
       "among cells descended from the same E7.0 founder lineage -- run identically on each "
       "tree so the two nulls are directly comparable; its outcome column reads `coupled` "
       "where the pair survives the within-clone null and is blank where that tree did not "
       "couple the pair under the primary dating. The manuscript quotes the BACKBONE column of "
       "this series.")
with open(FP + "tableS6_couplings_v8_paste.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow([TITLE] + [""] * (len(HDR) - 1))
    w.writerow([SUB] + [""] * (len(HDR) - 1))
    w.writerow([""] * len(HDR))
    w.writerow(BANNER); w.writerow(BANNER2); w.writerow(HDR); w.writerows(body)

# ---------------------------------------------------------------- checks
def compare(A, C, la, quoted):
    i = sorted(set(A) & set(C))
    a = [float(A[k]["divergence_age_both"]) if isinstance(A[k], dict) else A[k] for k in i]
    c = [C[k] for k in i]
    d = [y - x for x, y in zip(a, c)]
    print(f"  {la:<34} {len(i):3d} of {len(A):3d} survive | median {st.median(d):+.2f} d "
          f"| rho {spearmanr(a, c).statistic:.3f}   {quoted}")

both = sorted(set(F) & set(B))
af = [float(F[k]["divergence_age_both"]) for k in both]
ab = [float(B[k]["divergence_age_both"]) for k in both]
d = [y - x for x, y in zip(af, ab)]
print(f"placed tree couplings      {len(F)}")
print(f"backbone couplings         {len(B)}")
print(f"  coupled on both trees    {len(both)}   placed only {len(set(F)-set(B))}   "
      f"backbone only {len(set(B)-set(F))}")
print(f"table rows (the union)     {len(body)}   columns {len(HDR)}")
print(f"\nplaced vs backbone, global null, over the {len(both)} shared pairs:")
print(f"  median depth shift       {st.median(d):+.2f} d   (manuscript: -0.75)")
print(f"  Spearman rho             {spearmanr(af, ab).statistic:.3f}   (manuscript: 0.84)")
print(f"  deeper on backbone       {sum(1 for x in d if x < 0)} / {len(d)}")
print("\nprimary dating vs within-clone null (mixes pipelines):")
compare(B, BC, "backbone tree", "(old text: 55/84, +0.50, 0.76)")
compare(F, FC, "placed tree", "(not quoted)")
print("\nmatched-null series, one pipeline on both sides:")
compare(BG, BC, "backbone: global -> within-clone", "(new text: 55/86, +0.75, 0.77)")
compare(FG, FC, "placed: global -> within-clone", "(not quoted)")
compare(FG, BG, "placed -> backbone, global null", "(new text: rho 0.86)")
# outcome columns: 6 leading + 8-wide placed block -> index 13; backbone block -> index 21
OUT_P, OUT_B = 6 + 7, 6 + 8 + 7
assert HDR[OUT_P] == HDR[OUT_B] == "within-clone null: outcome", (HDR[OUT_P], HDR[OUT_B])
nsurv = lambda col: sum(1 for r in body if r[col] == "coupled")
print(f"\ncolumn check: 'coupled' rows  placed {nsurv(OUT_P)}   backbone {nsurv(OUT_B)}"
      f"   (expect 79 and 55)")
dates = [x["sort"] for x in rows]
print(f"coupling dates span        E{min(dates):g} - E{max(dates):g}")
print(f"distinct cell types        {len({c for x in rows for c in x['row'][:2]})}")
print(f"\nwrote {FP}tableS6_couplings_v8.csv")
print(f"wrote {FP}tableS6_couplings_v8_paste.csv")
