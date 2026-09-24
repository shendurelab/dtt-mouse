#!/usr/bin/env python3
# tree_layout.py
# Shared tree-geometry helpers for issue #2's two plotting scripts
# (plot_tree_tapes.py = Plot a, plot_tree_types.py = Plot b).
#
# These three functions encode the ONE invariant both panels must obey: the tree
# lines and the heatmap/strip rows have to register vertically, pixel for pixel.
# They live here (rather than copy-pasted into each plotter) so that invariant
# lives in a single place and the two panels can never silently drift apart.
# Everything else stays per-script; this module holds only the shared geometry.
#
# layout_tree_polar()/tree_segments_polar() (added for the v2 full-tree circular
# panel) reuse layout_tree()'s rectangular (depth, leaf-rank) coordinates and
# only remap them into polar space -- the leaf-rank -> angle map is affine, so
# it commutes with the existing mean-of-children recipe for internal-node y
# unchanged.

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


def layout_tree_polar(tree, tip_order, theta_gap_deg=20.0):
    # Circular counterpart of layout_tree(): same depth (-> radius) and
    # leaf-rank (-> angle) recipe, just remapped into polar space. theta_gap_deg
    # is left OPEN (no tips drawn there) for a radial scale bar, matching the
    # reference figure's style -- tips fill [theta_gap_deg/2, 360-theta_gap_deg/2].
    # Returns (theta, r) dicts keyed by ete3 node, in RADIANS and raw depth units
    # respectively -- same node-keying convention as layout_tree()'s (x, y).
    x, y = layout_tree(tree, tip_order)
    n_tips = len(tip_order)
    theta_range_deg = 360.0 - theta_gap_deg
    theta_start_deg = theta_gap_deg / 2.0
    denom = max(n_tips - 1, 1)  # avoid /0 for a degenerate single-tip tree
    # y -> theta is affine (a + b*y), so it commutes with layout_tree()'s
    # mean-of-children recipe for internal nodes: no need to recompute means.
    theta = {
        node: np.radians(theta_start_deg + (yv / denom) * theta_range_deg)
        for node, yv in y.items()
    }
    r = x
    return theta, r


def tree_segments_polar(tree, theta, r, n_arc_pts=8):
    # Polar counterpart of tree_segments(): for each node, one RADIAL segment
    # (fixed angle, r from parent's r to its own r) plus, for each internal
    # node, an ARC connecting its children's angles at its own radius. A
    # LineCollection only draws straight segments, so the arc is discretized
    # into n_arc_pts points (a polyline) rather than drawn as a true curve --
    # cheap and visually indistinguishable at the small angular spans typical
    # of a many-tip tree, and most nodes near the tips span far less than a
    # pixel anyway.
    segs = []
    for node in tree.traverse():
        if node.up is not None:
            t = theta[node]
            r0, r1 = r[node.up], r[node]
            segs.append([(r0 * np.cos(t), r0 * np.sin(t)),
                         (r1 * np.cos(t), r1 * np.sin(t))])
        if not node.is_leaf():
            thetas = [theta[c] for c in node.children]
            t_min, t_max = min(thetas), max(thetas)
            if t_max > t_min:
                rc = r[node]
                arc_theta = np.linspace(t_min, t_max, n_arc_pts)
                xs = rc * np.cos(arc_theta)
                ys = rc * np.sin(arc_theta)
                segs.extend(
                    [(xs[i], ys[i]), (xs[i + 1], ys[i + 1])]
                    for i in range(len(arc_theta) - 1)
                )
    return segs
