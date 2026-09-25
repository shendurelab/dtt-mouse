#!/usr/bin/env python3
# tree_layout.py
# Shared tree-geometry helpers for the tree panels in making_figures/.
#
# These functions encode the ONE invariant both panels must obey: the tree
# lines and the heatmap/strip rows have to register vertically, pixel for pixel.
# They live here (rather than copy-pasted into each plotter) so that invariant
# lives in a single place and the two panels can never silently drift apart.
# Everything else stays per-script; this module holds only the shared geometry.


import numpy as np
from ete3 import Tree


def load_tree_and_order(subtree_nwk, tiporder_txt):
    # 1. Load ONLY the subtree artifact (never the 100k tree).
    tree = Tree(subtree_nwk, format=1)
    # 2. Load the authoritative top->bottom tip order (heatmap/strip row order).
    with open(tiporder_txt) as fh:
        tip_order = [ln.strip() for ln in fh if ln.strip()]
    return tree, tip_order


def layout_tree(tree, tip_order):
    # Assign each node an (x, y):
    #   x = root-to-node distance (summed branch lengths) -> "NJ distance from root"
    #   y = leaf row index for leaves; internal nodes = mean of their children's y.
    # Leaves are pinned to their row in tip_order so the tree lines up with the
    # heatmap/strip rows exactly.
    row_of = {name: i for i, name in enumerate(tip_order)}

    # 1. x by accumulating branch lengths from the root down.
    x = {}
    for node in tree.traverse("preorder"):
        x[node] = (0.0 if node.up is None else x[node.up] + node.dist)

    # 2. y bottom-up: leaves get their fixed row, internals average their kids.
    y = {}
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            y[node] = row_of[node.name]
        else:
            y[node] = float(np.mean([y[c] for c in node.children]))
    return x, y


def tree_segments(tree, x, y):
    # Build the LineCollection segments for a classic rectangular cladogram:
    #   - one HORIZONTAL segment per node (from its parent's x to its own x)
    #   - one VERTICAL connector per internal node spanning its children's y.
    # Relies on x/y being populated for every node (as layout_tree does).
    segs = []
    for node in tree.traverse():
        if node.up is not None:
            segs.append([(x[node.up], y[node]), (x[node], y[node])])  # horizontal
        if not node.is_leaf():
            ys = [y[c] for c in node.children]
            segs.append([(x[node], min(ys)), (x[node], max(ys))])     # vertical
    return segs
