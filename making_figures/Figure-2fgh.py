#!/usr/bin/env python3
# Figure-2fgh.py
# Fig. 2F: time tree + TAPE-alignment heatmap, on the 655,701-tip backbone
# tree (same tree as single_fig_3_panels_de_analysis.R's panels J/K).
# Fig. 2G/H: two nested cell-type zoom-ins into that tree (~1000 -> ~100 cells).
#
# Needs (see the "zoom clade selection" section below): the two zoom clades
# (Fig. 2G, 2H) are picked by a diversity-maximizing search over the whole
# backbone tree that is NOT implemented here (kept in R, using ape/phangorn --
# see Figure-2fgh_select_zoom.R). Fig. 2F's tree+TAPE panel does not need it
# and always runs; 2G/2H need it and will stop with a clear message if it
# hasn't been run yet.
#
# Run from the repo root:
#   python3 making_figures/Figure-2fgh.py

import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans"]
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch, Rectangle
from ete3 import Tree

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
from tree_layout import layout_tree, tree_segments          # noqa: E402
from tape_alignment import N_SITES, read_integration_columns  # noqa: E402

#-----------------------------------------------------------------------------#
# 0. inputs (repo-root relative)                                              #
#-----------------------------------------------------------------------------#
TREE_NWK = "tree_building/results/3-dated-tree/merged_minB2h_lineage_constrained.nwk"
TAPE_B1  = "support_data/e3v8.B1_tape_consensus.tsv.gz"
TAPE_B2  = "support_data/e3v8.B2_tape_consensus.tsv.gz"
METADATA = "support_data/cell_metadata.v8.txt.gz"   # cell_id, major_trajectory, celltype, ...

for f in (TREE_NWK, TAPE_B1, TAPE_B2, METADATA):
    assert os.path.exists(f), f"missing required input: {f}"

OUTDIR = "making_figures/Figure-2fgh_output"
os.makedirs(OUTDIR, exist_ok=True)
ZDIR = os.path.join(OUTDIR, "zoom")
ZOOM = dict(
    zoom1_nwk=os.path.join(ZDIR, "zoom1.nwk"),
    zoom2_nwk=os.path.join(ZDIR, "zoom2.nwk"),
    zoom1_day=os.path.join(ZDIR, "zoom1_day_offset.txt"),
    zoom2_day=os.path.join(ZDIR, "zoom2_day_offset.txt"),
    zoom1_tips=os.path.join(ZDIR, "zoom1_tips.txt"),
    zoom2_tips=os.path.join(ZDIR, "zoom2_tips.txt"),
)
zoom_ready = all(os.path.exists(f) for f in ZOOM.values())
if not zoom_ready:
    print(f"[zoom] clade-selection artifacts not found under {ZDIR}/ -- "
          "Fig. 2F will draw without its nesting-cue highlight box, and "
          "Fig. 2G/2H will be skipped. Run the (not-yet-ported) zoom-clade "
          "selection step first -- see mouse_sprint's "
          "src/3-plot-tree/02_select_nested_zoom.R (ape/phangorn descendant-"
          "count search); not reimplemented here in Python to avoid a "
          "tie-break order picking a different clade.")

print(f"tree     : {TREE_NWK}")
print(f"outdir   : {OUTDIR}")

if zoom_ready:
    with open(ZOOM["zoom1_day"]) as fh:
        ZOOM1_DAY = float(fh.read().strip())
    with open(ZOOM["zoom2_day"]) as fh:
        ZOOM2_DAY = float(fh.read().strip())

#-----------------------------------------------------------------------------#
# Fig. 2F: full time tree + TAPE-alignment heatmap                            #
#-----------------------------------------------------------------------------#
F_PALETTE = [
    ("CTG",  "#e7298a"),
    ("GCC",  "#a6cee3"),
    ("AAG",  "#1f78b4"),
    ("CAC",  "#b2df8a"),
    ("GATG", "#33a02c"),
    ("CCC",  "#fb9a99"),
    ("GGC",  "#e31a1c"),
    ("ACG",  "#fdbf6f"),
    ("GAA",  "#ff7f00"),
    ("ACC",  "#cab2d6"),
    ("ACA",  "#6a3d9a"),
    ("CCG",  "#ffff99"),
    ("ACT",  "#b15928"),
    ("U",    "#e1e0d9"),  # unedited site -- an observed call, not missing data
]
F_SYMBOL_TO_IDX = {sym: i for i, (sym, _) in enumerate(F_PALETTE)}
F_IDX_OTHER = len(F_PALETTE)
F_IDX_NA = len(F_PALETTE) + 1
F_COLOR_OTHER = "#999999"
F_COLOR_NA = "#ffffff"
F_MAX_MISSING_TIP_FRAC = 0.01


def build_alignment_matrix(tape_tsvs, tip_order, integrations):
    frames = [pd.read_csv(p, sep="\t", dtype=str) for p in tape_tsvs]
    for p, frame in zip(tape_tsvs, frames):
        assert frame.columns[0] == "cell_id", f"{p}: first column must be cell_id"
    df = pd.concat(frames, ignore_index=True)
    dup_ids = df["cell_id"][df["cell_id"].duplicated()]
    assert dup_ids.empty, (
        f"{len(dup_ids)} cell_id(s) appear in more than one tape table "
        f"(e.g. {dup_ids.iloc[0]!r}) -- expected disjoint per-side tables")
    id_set = set(df["cell_id"])
    assert tip_order[0] in id_set, f"tip {tip_order[0]!r} not found in tape table -> id mismatch"
    df = df.set_index("cell_id").reindex(tip_order)

    n_rows = len(tip_order)
    n_cols = len(integrations) * N_SITES
    mat = np.full((n_rows, n_cols), F_IDX_NA, dtype=np.int16)
    n_calls = 0
    n_other = 0
    n_missing_cell = 0

    for r, cell_id in enumerate(tip_order):
        if cell_id not in id_set:
            n_missing_cell += 1
            continue
        for k, integ in enumerate(integrations):
            value = df.at[cell_id, integ]
            col0 = k * N_SITES
            if value is None or (isinstance(value, float) and np.isnan(value)) \
                    or value == "" or value == "NA" or value == "nan":
                continue
            sites = value.split("|")
            assert len(sites) == N_SITES, f"tape {integ} for {cell_id} has {len(sites)} sites: {value!r}"
            for s, sym in enumerate(sites):
                n_calls += 1
                idx = F_SYMBOL_TO_IDX.get(sym)
                if idx is None:
                    mat[r, col0 + s] = F_IDX_OTHER
                    n_other += 1
                else:
                    mat[r, col0 + s] = idx

    assert n_missing_cell <= F_MAX_MISSING_TIP_FRAC * n_rows, (
        f"{n_missing_cell}/{n_rows} tips absent from the tape table "
        f"(> {F_MAX_MISSING_TIP_FRAC:.0%}) -> likely cell_id mismatch, not coverage")

    total_sites = n_rows * n_cols
    n_na = int((mat == F_IDX_NA).sum())
    stats = dict(
        n_rows=n_rows, n_missing_cell=n_missing_cell,
        na_frac=n_na / total_sites,
        other_frac=(n_other / n_calls) if n_calls else 0.0,
        n_calls=n_calls, n_other=n_other,
    )
    return mat, stats


def tape_colormap():
    colors = [hexcol for _, hexcol in F_PALETTE] + [F_COLOR_OTHER, F_COLOR_NA]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(len(colors) + 1) - 0.5, len(colors))
    return cmap, norm


def add_tape_legend(fig, n_other, **legend_kwargs):
    handles = [
        Patch(facecolor=hexcol, edgecolor="none", label=("Unedited" if sym == "U" else sym))
        for sym, hexcol in F_PALETTE
    ]
    if n_other > 0:
        handles.append(Patch(facecolor=F_COLOR_OTHER, edgecolor="none", label="other"))
    fig.legend(handles=handles, frameon=False, **legend_kwargs)


F_TAG = "F"
F_STEM = f"panel_{F_TAG}_timetree_tapes"
F_OUT_PNG = os.path.join(OUTDIR, f"{F_STEM}.png")
F_OUT_PDF = os.path.join(OUTDIR, f"{F_STEM}.pdf")

print(f"\n[panel F] loading backbone tree from {TREE_NWK} ...")
f_tree = Tree(TREE_NWK, format=1)
f_tip_order = [leaf.name for leaf in f_tree.get_leaves()]
f_n_tips = len(f_tip_order)
print(f"[panel F] tips: {f_n_tips:,}")

f_x, f_y = layout_tree(f_tree, f_tip_order)
f_segs = tree_segments(f_tree, f_x, f_y)

# The merged root is a synthetic day-0 attachment point (both blastomere
# stems pinned to the same calibrated depth), so the default rendering draws
# a fork right at x=0 -- replace it with a single stem out to the calibrated
# split depth, matching how the split is actually dated.
f_root = f_tree
if len(f_root.children) == 2 and abs(f_x[f_root.children[0]] - f_x[f_root.children[1]]) < 1e-9:
    c0, c1 = f_root.children
    split_x = f_x[c0]
    y0, y1 = f_y[c0], f_y[c1]
    ylo, yhi = min(y0, y1), max(y0, y1)
    default_segs = [
        [(f_x[f_root], ylo), (f_x[f_root], yhi)],
        [(f_x[f_root], y0), (f_x[c0], y0)],
        [(f_x[f_root], y1), (f_x[c1], y1)],
    ]
    f_segs = [s for s in f_segs if s not in default_segs]
    f_segs.append([(f_x[f_root], f_y[f_root]), (split_x, f_y[f_root])])
    f_segs.append([(split_x, ylo), (split_x, yhi)])

print(f"[panel F] tree line segments: {len(f_segs):,}")

# integrations read off the first TSV's header, filtered to ACGT-only columns
# so the appended QC passthrough columns (n_loci, pass_qc, ...) are excluded
# without hardcoding an exclude-list that could drift.
f_integrations = [c for c in read_integration_columns(TAPE_B1) if re.fullmatch(r"[ACGT]+", c)]
f_mat, f_stats = build_alignment_matrix([TAPE_B1, TAPE_B2], f_tip_order, f_integrations)
f_n_cols = len(f_integrations) * N_SITES
print(f"[panel F] integrations: {len(f_integrations)}")
print(f"[panel F] tips missing from tape table: {f_stats['n_missing_cell']}")
print(f"[panel F] observed site-calls: {f_stats['n_calls']:,}")
print(f"[panel F] 'other' fraction: {f_stats['other_frac']:.4%} ({f_stats['n_other']} calls)")
print(f"[panel F] NA fraction: {f_stats['na_frac']:.2%}")

f_cmap, f_norm = tape_colormap()

F_SCALE = 0.736   # matches the other Figure-2 panels' historical on-page text size
f_body_fontsize = 9.0 / F_SCALE
f_title_fontsize = 10.5 / F_SCALE
f_tag_fontsize = 12.0 / F_SCALE

f_fig = plt.figure(figsize=(17, 16))
f_gs_right = 0.925
f_gs = f_fig.add_gridspec(1, 2, width_ratios=[1.0, 2.4], wspace=0.02,
                          left=0.045, right=f_gs_right, top=0.99, bottom=0.06)
f_ax_tree = f_fig.add_subplot(f_gs[0, 0])
f_ax_heat = f_fig.add_subplot(f_gs[0, 1])

f_lc = LineCollection(f_segs, colors="#000000", linewidths=0.15, rasterized=True)
f_ax_tree.add_collection(f_lc)
f_xmax = max(f_x.values())
f_ax_tree.set_xlim(0, f_xmax * 1.02)
f_ax_tree.set_ylim(f_n_tips - 0.5, -0.5)
f_ax_tree.set_yticks([])
f_regular_ticks = [t for t in np.arange(1.5, f_xmax, 2.0) if f_xmax - t > 0.6]
f_ticks = f_regular_ticks + [f_xmax]
f_ax_tree.set_xticks(f_ticks)
f_ax_tree.set_xticklabels([f"{t:.1f}" for t in f_regular_ticks] + [f"{f_xmax:.1f}"],
                          fontsize=f_body_fontsize)
f_ax_tree.set_xlabel("Time (days)", fontsize=f_title_fontsize)
for spine in ("top", "right", "left"):
    f_ax_tree.spines[spine].set_visible(False)

# Nesting cue: red box around zoom1's tip rows, showing where Fig. 2G zooms
# in (a matching box is drawn across the heatmap below). Only when ready.
f_highlight_rows = None
if zoom_ready:
    f_hl = set(l.strip() for l in open(ZOOM["zoom1_tips"]) if l.strip())
    f_rows = [i for i, cid in enumerate(f_tip_order) if cid in f_hl]
    if f_rows:
        f_highlight_rows = (min(f_rows), max(f_rows))
        r0, r1 = f_highlight_rows
        f_ax_tree.add_patch(Rectangle((ZOOM1_DAY, r0 - 0.5), f_xmax * 1.02 - ZOOM1_DAY, (r1 - r0 + 1),
                                      fill=False, edgecolor="#d1495b", linewidth=1.6, zorder=5))
        print(f"[panel F] highlight: {len(f_rows)} tips, rows {r0}-{r1}, box left E{ZOOM1_DAY:.2f}")

f_ax_heat.imshow(
    f_mat, aspect="auto", interpolation="nearest", cmap=f_cmap, norm=f_norm,
    origin="upper", extent=[0, f_n_cols, f_n_tips - 0.5, -0.5], rasterized=True,
)
for k in range(1, len(f_integrations)):
    f_ax_heat.axvline(k * N_SITES, color="white", linewidth=1.1)
f_tick_centers = [k * N_SITES + N_SITES / 2 for k in range(len(f_integrations))]
f_ax_heat.set_xticks(f_tick_centers)
f_ax_heat.set_xticklabels([integ[:6] for integ in f_integrations], rotation=90, fontsize=f_body_fontsize)
f_ax_heat.set_yticks([])
f_ax_heat.set_ylim(f_n_tips - 0.5, -0.5)
for spine in f_ax_heat.spines.values():
    spine.set_visible(False)

if f_highlight_rows is not None:
    r0, r1 = f_highlight_rows
    f_ax_heat.add_patch(Rectangle((0, r0 - 0.5), f_n_cols, (r1 - r0 + 1),
                                  fill=False, edgecolor="#d1495b", linewidth=1.6, zorder=5))

F_LEGEND_OUT = os.path.join(OUTDIR, "tape_legend.tsv")
with open(F_LEGEND_OUT, "w") as fh:
    for sym, hexcol in F_PALETTE:
        fh.write(f"{'Unedited' if sym == 'U' else sym}\t{hexcol}\n")
    if f_stats["n_other"] > 0:
        fh.write(f"other\t{F_COLOR_OTHER}\n")
print(f"[panel F] wrote tape-key -> {F_LEGEND_OUT}")
add_tape_legend(f_fig, f_stats["n_other"], loc="center left", ncol=1,
                fontsize=f_body_fontsize, handlelength=1.1,
                labelspacing=0.55, bbox_to_anchor=(f_gs_right + 0.012, 0.5))

f_ax_tree.annotate(F_TAG, xy=(0, 1), xycoords="axes fraction", xytext=(-2, 6),
                   textcoords="offset points", fontsize=f_tag_fontsize,
                   fontweight="bold", ha="left", va="bottom", color="black")

f_fig.savefig(F_OUT_PNG, dpi=300, bbox_inches="tight")
f_fig.savefig(F_OUT_PDF, dpi=300, bbox_inches="tight")
print(f"[panel F] wrote {F_OUT_PNG}\n[panel F] wrote {F_OUT_PDF}")

#-----------------------------------------------------------------------------#
# Fig. 2G, 2H: nested cell-type zoom-ins                                       #
#-----------------------------------------------------------------------------#
GH_TYPE_COLUMN = "major_trajectory"
GH_COLOR_UNKNOWN = "#cccccc"
GH_TREE_LINE_COLOR = "#111111"
# Same palette as making_figures/Figure-2.R's major_trajectory_color_plate
# (Fig. 2E) -- duplicated here since that file defines it inline in R rather
# than in a sourceable module; keep the two in sync if it changes there.
GH_TRAJECTORY_COLORS = {
    "Neuroectoderm_and_glia":            "#f96100",
    "Intermediate_neuronal_progenitors": "#2e0ab7",
    "Eye_and_other":                     "#00d450",
    "Ependymal_cells":                   "#b75bff",
    "CNS_neurons":                       "#e5c000",
    "Mesoderm":                          "#bb46c5",
    "Definitive_erythroid":              "#dc453e",
    "Epithelium":                        "#af9fb6",
    "Endothelium":                       "#00a34e",
    "Muscle_cells":                      "#ffa1f5",
    "Hepatocytes":                       "#185700",
    "White_blood_cells":                 "#7ca0ff",
    "Neural_crest_PNS_glia":             "#fff167",
    "Adipocytes":                        "#7f3e39",
    "Primitive_erythroid":               "#ffa9a1",
    "Neural_crest_PNS_neurons":          "#b5ce92",
    "T_cells":                           "#ff9d47",
    "Lung_and_airway":                   "#02b0d1",
    "Intestine":                         "#ff007a",
    "B_cells":                           "#01b7a6",
    "Olfactory_sensory_neurons":         "#e6230b",
    "Cardiomyocytes":                    "#643e8c",
    "Oligodendrocytes":                  "#916e00",
    "Mast_cells":                        "#005361",
    "Megakaryocytes":                    "#3f283d",
    "Testis_and_adrenal":                "#585d3b",
}


def pretty(name):
    if "_PNS_" in name:
        base, tail = name.split("_PNS_")
        return f"{base.replace('_', ' ')} (PNS {tail.replace('_', ' ')})"
    return name.replace("_and_", " & ").replace("_", " ")


def global_category_colors(meta_txt):
    categories = list(GH_TRAJECTORY_COLORS.keys()) + ["unknown"]
    colors = list(GH_TRAJECTORY_COLORS.values()) + [GH_COLOR_UNKNOWN]
    meta = pd.read_csv(meta_txt, sep="\t", dtype=str)
    return categories, colors, dict(zip(meta["cell_id"], meta[GH_TYPE_COLUMN]))


def day_ticks(day_offset, xmax, step=1.0):
    tip_day = xmax + day_offset
    first = np.ceil(day_offset - 0.5) + 0.5
    regular = [t for t in np.arange(first, tip_day, step) if tip_day - t > 0.4]
    ticks_abs = regular + [tip_day]
    return [t - day_offset for t in ticks_abs], [f"{t:.1f}" for t in ticks_abs]


def plot_subclade_zoom(tree_nwk, day_offset, tag, label, legend_out=None,
                       highlight_tips=None, highlight_day=None):
    tree = Tree(tree_nwk, format=1)
    leaves = tree.get_leaves()
    tip_order = [leaf.name for leaf in leaves]
    n_tips = len(tip_order)
    x, y = layout_tree(tree, tip_order)
    segs = tree_segments(tree, x, y)
    xmax = max(x.values())

    categories, colors, type_of = global_category_colors(METADATA)
    code_of = {c: i for i, c in enumerate(categories)}
    unknown_code = code_of["unknown"]
    codes = np.array([code_of.get(type_of.get(cid, "unknown"), unknown_code)
                      for cid in tip_order], dtype=np.int16)
    present = sorted(set(codes.tolist()))
    n_present = sum(1 for c in present if categories[c] != "unknown")
    print(f"[panel {tag}] {label}: {n_tips} tips, {n_present} cell types, clade root E{day_offset:.2f}")

    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(len(colors) + 1) - 0.5, len(colors))

    fig = plt.figure(figsize=(3.6, 9.16))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.07], wspace=0.02,
                          left=0.02, right=0.90, top=0.965, bottom=0.055)
    ax_tree = fig.add_subplot(gs[0, 0])
    ax_strip = fig.add_subplot(gs[0, 1])

    lw = 0.25 if n_tips > 400 else 0.5
    ax_tree.add_collection(LineCollection(segs, colors=GH_TREE_LINE_COLOR, linewidths=lw, rasterized=True))
    ax_tree.set_xlim(0, xmax * 1.02)
    ax_tree.set_ylim(n_tips - 0.5, -0.5)
    ax_tree.set_yticks([])
    tl, labs = day_ticks(day_offset, xmax)
    ax_tree.set_xticks(tl)
    ax_tree.set_xticklabels(labs, fontsize=7)
    ax_tree.set_xlabel("Time (days)", fontsize=8)
    for sp in ("top", "right", "left"):
        ax_tree.spines[sp].set_visible(False)
    ax_tree.annotate(tag, xy=(0, 1), xycoords="axes fraction", xytext=(-2, 6),
                     textcoords="offset points", fontsize=12, fontweight="bold",
                     ha="right", va="bottom")

    if highlight_tips and highlight_day is not None:
        hl = set(l.strip() for l in open(highlight_tips) if l.strip())
        rows = [i for i, cid in enumerate(tip_order) if cid in hl]
        if rows:
            r0, r1 = min(rows), max(rows)
            x0 = highlight_day - day_offset
            ax_tree.add_patch(Rectangle((x0, r0 - 0.5), xmax * 1.02 - x0, (r1 - r0 + 1),
                                        fill=False, edgecolor="#d1495b", linewidth=1.2, zorder=5))
            print(f"[panel {tag}]   highlight: {len(rows)} tips, rows {r0}-{r1}, box left E{highlight_day:.2f}")

    ax_strip.imshow(codes.reshape(-1, 1), aspect="auto", interpolation="nearest",
                    cmap=cmap, norm=norm, origin="upper",
                    extent=[0, 1, n_tips - 0.5, -0.5], rasterized=True)
    ax_strip.set_xticks([])
    ax_strip.set_yticks([])
    ax_strip.set_ylim(n_tips - 0.5, -0.5)
    for sp in ax_strip.spines.values():
        sp.set_visible(False)

    if legend_out:
        with open(legend_out, "w") as fh:
            for c in present:
                if categories[c] == "unknown":
                    continue
                fh.write(f"{pretty(categories[c])}\t{colors[c]}\n")
        print(f"[panel {tag}]   wrote shared-key categories -> {legend_out}")

    stem = f"panel_{tag}_{label}"
    out_png = os.path.join(OUTDIR, f"{stem}.png")
    out_pdf = os.path.join(OUTDIR, f"{stem}.pdf")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")

    out_csv = os.path.join(OUTDIR, f"{stem}.csv")
    tip_x = [x[leaf] for leaf in leaves]
    pd.DataFrame({
        "tip_label": tip_order,
        "row": range(n_tips),
        "time_since_clade_root_days": tip_x,
        "absolute_embryonic_day": [v + day_offset for v in tip_x],
        "major_trajectory": [categories[c] for c in codes],
    }).to_csv(out_csv, index=False)

    print(f"[panel {tag}]   wrote {out_png}")
    return fig


if zoom_ready:
    GH_LEGEND_OUT = os.path.join(ZDIR, "zoom_legend.tsv")
    p_g = plot_subclade_zoom(
        ZOOM["zoom1_nwk"], ZOOM1_DAY, tag="G", label="zoom1",
        legend_out=GH_LEGEND_OUT,
        highlight_tips=ZOOM["zoom2_tips"], highlight_day=ZOOM2_DAY,
    )
    p_h = plot_subclade_zoom(
        ZOOM["zoom2_nwk"], ZOOM2_DAY, tag="H", label="zoom2",
    )
else:
    print("\n[panels G/H] skipped -- zoom-clade selection not yet run (see note above)")

print(f"\ndone -> {OUTDIR}")
