"""Per-integration lineage trees.

Each TAPE integration (barcode) gets a prefix-trie over its ordered edit chains
(tape.editchain): the root->tip paths are the lineage relationships, tip depth =
number of sites written. This module lays that trie out with matplotlib -- a
6-site nomenclature strip and a weight bar per tip -- and renders a grid over all
integrations. It replaces the two divergent renderers (lineage_trees_v2 /
gen_cell_trees) with one function.
"""
from __future__ import annotations
import math
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from config import N_SITES, PALETTE, UNEDITED_COLOR, OTHER_COLOR, TREE_TOPK
from tape.editchain import build_trie, nodes_per_level


def insertion_colors(per_bc: dict) -> dict:
    """Assign colors to the genome-wide most abundant insertions (stable order)."""
    tot = Counter()
    for counts in per_bc.values():
        for pat, c in counts.items():
            for s in pat:
                if s != "U":
                    tot[s] += c
    top = [k for k, _ in tot.most_common(len(PALETTE))]
    return {k: PALETTE[i] for i, k in enumerate(top)}


def color_of(symbol, colors) -> str:
    if symbol == "U":
        return UNEDITED_COLOR
    return colors.get(symbol, OTHER_COLOR)


def _pad(chain) -> tuple:
    """Edit chain -> full 6-site pattern for the nomenclature strip."""
    return tuple(chain) + ("U",) * (N_SITES - len(chain))


def draw_integration_tree(ax, bc: str, counts: dict, colors: dict,
                          topk: int = TREE_TOPK, weight_name: str = "reads"):
    """Draw one integration's lineage trie on `ax`.

    counts: {pattern: weight} for this barcode (already de-noised). The top `topk`
    tips by weight are drawn; node widths per level are annotated from the full set.
    """
    full_root = build_trie(counts)
    npl = nodes_per_level(full_root, N_SITES)

    top_items = sorted(counts.items(), key=lambda x: -x[1])[:topk]
    shown = sum(w for _, w in top_items)
    total = sum(counts.values()) or 1
    root = build_trie(dict(top_items))

    # DFS tip order -> y rows; internal node y = mean of its subtree tips
    rows = []
    def dfs(node):
        if node.tip > 0:
            rows.append(node)
        for _, ch in sorted(node.children.items(), key=lambda kv: -kv[1].subtree):
            dfs(ch)
    dfs(root)
    for i, node in enumerate(rows):
        node.y = float(i)

    def set_y(node):
        ys = [node.y] if node.tip > 0 else []
        for ch in node.children.values():
            ys += set_y(ch)
        node.y = float(np.mean(ys)) if ys else 0.0
        return ys
    set_y(root)

    n = len(rows)
    wmax = max((nd.subtree for nd in rows), default=1)
    tipw = [nd.tip for nd in rows]
    lo = math.log10(max(1, min(tipw))) if tipw else 0
    hi = math.log10(max(tipw)) if tipw else 1

    # edges, line width ~ log weight of the child subtree
    def walk(node):
        for ch in node.children.values():
            ax.plot([node.depth, ch.depth], [node.y, ch.y], "-", color="#0C6B74",
                    lw=0.5 + 2.2 * (math.log10(max(2, ch.subtree)) / math.log10(max(2, wmax))),
                    alpha=.65, solid_capstyle="round", zorder=1)
            walk(ch)
    walk(root)

    sx = N_SITES + 0.5
    cellw = 0.34
    barx = sx + N_SITES * cellw + 0.25
    for nd in rows:
        pat = _pad(nd.chain())
        ax.scatter([nd.depth], [nd.y], s=8, color="#16242E", zorder=3)
        for i, s in enumerate(pat):
            ax.add_patch(Rectangle((sx + i * cellw, nd.y - 0.4), cellw * 0.92, 0.8,
                                   facecolor=color_of(s, colors), edgecolor="white",
                                   lw=0.4, zorder=2))
        frac = (math.log10(max(1, nd.tip)) - lo) / max(0.1, (hi - lo))
        ax.add_patch(Rectangle((barx, nd.y - 0.34), 0.2 + 2.6 * frac, 0.68,
                               facecolor="#16242E", edgecolor="none", zorder=2))
        label = f"{int(nd.tip):,}"
        ax.text(barx + 0.25 + 2.6 * frac + 0.05, nd.y, label, va="center",
                fontsize=5.0, color="#445")

    ax.set_xlim(-0.3, barx + 3.4)
    ax.set_ylim(n - 0.3, -0.7)
    ax.set_title(f"{bc}  -  top {min(topk, len(counts))} of {len(counts)} lineages "
                 f"({shown / total * 100:.0f}% of {weight_name})",
                 fontsize=8.0, loc="left")
    ax.set_xticks(range(N_SITES + 1))
    ax.set_xticklabels(["root\n(1)"] + [f"site {i}\n({npl[i]} nodes)"
                                        for i in range(1, N_SITES + 1)], fontsize=6.0)
    ax.set_yticks([])
    ax.tick_params(length=2)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)


def render_all(per_bc_clean: dict, out_path, weight_name: str = "reads",
               title: str = "Per-integration TAPE lineage trees"):
    """Render a grid of per-integration trees (one panel per barcode) + legend."""
    bcs = sorted(per_bc_clean, key=lambda b: -sum(per_bc_clean[b].values()))
    colors = insertion_colors(per_bc_clean)

    nrow = (len(bcs) + 1) // 2
    fig, axes = plt.subplots(nrow, 2, figsize=(15, 3.0 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for k, bc in enumerate(bcs):
        draw_integration_tree(axes[k], bc, per_bc_clean[bc], colors, weight_name=weight_name)
    for j in range(len(bcs), len(axes)):
        axes[j].axis("off")

    handles = [plt.Line2D([0], [0], marker="s", ls="", mfc=colors[k], mec="white",
                          ms=9, label=k) for k in colors]
    handles.append(plt.Line2D([0], [0], marker="s", ls="", mfc=UNEDITED_COLOR,
                              mec="white", ms=9, label="unedited"))
    axes[-1].legend(handles=handles, loc="center", ncol=2, fontsize=8,
                    title="insertion (monomer barcode)", frameon=False)
    axes[-1].axis("off")

    fig.suptitle(title + " (tip depth = # sites edited)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    fig.savefig(str(out_path), dpi=140)
    plt.close(fig)
    return out_path
