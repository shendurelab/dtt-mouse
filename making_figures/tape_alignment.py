#!/usr/bin/env python3
# tape_alignment.py
# Shared TAPE-alignment-matrix + symbol-palette helpers for the tree+heatmap
# panels (plot_tree_tapes.py, plot_subclade_timetree_tapes.py) -- kept in one
# place for the same reason tree_layout.py exists (see that module's own
# docstring): so a fix to the palette or the NA-handling contract can be
# applied once instead of maintained in two copies that can silently drift
# apart. Each caller still owns its OWN integration list (v1's are hardcoded
# to a specific reference order; v2's are read off the TSV header, since the
# two versions' integration ID strings differ in length) -- this module only
# owns the parts that must stay identical between panels.

import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

N_SITES = 6  # monomer sites per integration

# ---- symbol palette (hex values eyedropped from sc_fig_nj_clone.png; tab20) ----
# Ordered exactly as the reference legend (dark then light within each hue pair).
# Index into this list is the value stored in the heatmap matrix.
PALETTE = [
    ("GCC",   "#1f77b4"),  # blue dark
    ("AAG",   "#aec7e8"),  # blue light
    ("CCC",   "#ff7f0e"),  # orange dark
    ("CAC",   "#ffbb78"),  # orange light
    ("ACT",   "#2ca02c"),  # green dark
    ("GATG",  "#98df8a"),  # green light
    ("ACA",   "#d62728"),  # red dark
    ("CTG",   "#ff9896"),  # red light
    ("ACC",   "#9467bd"),  # purple dark
    ("GGC",   "#c5b0d5"),  # purple light
    ("ACG",   "#8c564b"),  # brown dark
    ("CCG",   "#c49c94"),  # brown light
    ("ATT",   "#e377c2"),  # pink dark
    ("GAA",   "#f7b6d2"),  # pink light
    ("TGATG", "#7f7f7f"),  # gray dark
    ("TGAT",  "#c7c7c7"),  # gray light
    ("U",     "#e3e3e3"),  # unedited, very light gray
]
SYMBOL_TO_IDX = {sym: i for i, (sym, _) in enumerate(PALETTE)}

# Two extra reserved codes appended after the palette entries:
IDX_OTHER = len(PALETTE)      # any symbol not in the palette (rare tail)
IDX_NA = len(PALETTE) + 1     # integration not observed / missing site
COLOR_OTHER = "#606060"       # neutral gray, distinct from TGATG/TGAT/U
COLOR_NA = "#ffffff"          # white background (as in the reference)

# Maximum fraction of tips we tolerate being absent from the tape table before
# treating it as an ID-scheme mismatch (a bug) rather than genuine low coverage.
MAX_MISSING_TIP_FRAC = 0.01


def read_integration_columns(tape_tsv):
    # Integration IDs straight off the TSV header (cheap: nrows=0 reads only
    # the header row) -- the version-proof way to get the list, since v1's
    # IDs are 13 characters and v2's are 12; hardcoding either breaks the
    # other. Callers that need a specific fixed reference order (v1's
    # existing panel) keep their own hardcoded list instead of calling this.
    header = pd.read_csv(tape_tsv, sep="\t", nrows=0, dtype=str)
    assert header.columns[0] == "cell_id"
    return list(header.columns[1:])


def build_alignment_matrix(tape_tsv, tip_order, integrations):
    # 0. Read the TAPE table bare (all strings so tokens are never coerced).
    df = pd.read_csv(tape_tsv, sep="\t", dtype=str)
    # 1. Assert cell_id matches tip label format (spot-check the first tip).
    assert df.columns[0] == "cell_id"
    id_set = set(df["cell_id"])
    assert tip_order[0] in id_set, (
        f"tip {tip_order[0]!r} not found in tape table -> id mismatch")
    # 2. Subset+join to our tips only, in row order == tip_order.
    df = df.set_index("cell_id").reindex(tip_order)  # missing tips -> all-NaN row

    # 3. Parse each integration string into an N_SITES-column integer-code
    #    block; codes index PALETTE / OTHER / NA.
    n_rows = len(tip_order)
    n_cols = len(integrations) * N_SITES
    mat = np.full((n_rows, n_cols), IDX_NA, dtype=np.int16)  # default = missing
    n_calls = 0        # observed (non-NA) site calls seen
    n_other = 0        # calls sent to "other" (symbol not in palette)
    n_missing_cell = 0  # tips absent from the tape table entirely

    for r, cell_id in enumerate(tip_order):
        if cell_id not in id_set:
            n_missing_cell += 1
            continue
        for k, integ in enumerate(integrations):
            value = df.at[cell_id, integ]
            col0 = k * N_SITES
            # a whole-cell NA = integration NOT observed -> N_SITES missing sites.
            # (Do NOT split the literal "NA" string on '|'.)
            # We read the TSV with dtype=str, so in practice a not-observed cell
            # arrives as the literal string "NA"; the None/NaN/""/"nan" branches
            # are cheap belt-and-suspenders guards for other loaders/edge rows.
            if value is None or (isinstance(value, float) and np.isnan(value)) \
                    or value == "" or value == "NA" or value == "nan":
                continue  # leave those N_SITES columns at IDX_NA
            sites = value.split("|")
            # data contract: an observed tape has exactly N_SITES sites.
            assert len(sites) == N_SITES, \
                f"tape {integ} for {cell_id} has {len(sites)} sites: {value!r}"
            for s, sym in enumerate(sites):
                n_calls += 1
                idx = SYMBOL_TO_IDX.get(sym)
                if idx is None:
                    mat[r, col0 + s] = IDX_OTHER
                    n_other += 1
                else:
                    mat[r, col0 + s] = idx

    # Fail fast if a large share of tips is missing from the tape table: that is
    # an ID-scheme mismatch (a bug), not real low coverage, and it would silently
    # render as fully-"unedited" rows. A tiny miss count is tolerated.
    assert n_missing_cell <= MAX_MISSING_TIP_FRAC * n_rows, (
        f"{n_missing_cell}/{n_rows} tips absent from the tape table "
        f"(> {MAX_MISSING_TIP_FRAC:.0%}) -> likely cell_id mismatch, not coverage")

    total_sites = n_rows * n_cols
    n_na = int((mat == IDX_NA).sum())
    stats = dict(
        n_rows=n_rows,
        n_missing_cell=n_missing_cell,
        na_frac=n_na / total_sites,
        other_frac=(n_other / n_calls) if n_calls else 0.0,
        n_calls=n_calls,
        n_other=n_other,
    )
    return mat, stats


def tape_colormap():
    # Discrete colormap: palette colors + OTHER + NA. BoundaryNorm maps
    # integer code v -> color slot v exactly.
    colors = [hexcol for _, hexcol in PALETTE] + [COLOR_OTHER, COLOR_NA]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(len(colors) + 1) - 0.5, len(colors))
    return cmap, norm


def add_tape_legend(fig, **legend_kwargs):
    # The palette symbols + "other" (NA/white is background, not legended,
    # matching the reference figure). Layout (bbox_to_anchor/fontsize/ncol/...)
    # differs per panel, so it's forwarded rather than fixed here.
    handles = [Patch(facecolor=hexcol, edgecolor="none", label=sym)
               for sym, hexcol in PALETTE]
    handles.append(Patch(facecolor=COLOR_OTHER, edgecolor="none", label="other"))
    fig.legend(handles=handles, frameon=False, **legend_kwargs)
