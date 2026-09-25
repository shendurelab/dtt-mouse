#!/usr/bin/env python3
"""Build mergedtree_dttpq_v8.npz from the Newick tree that ships in support_data.

Several tree_analysis scripts read the merged, placed, dated tree as a flat
array bundle rather than re-parsing Newick every time. That bundle is ~93 MB,
which is a lot to carry in a git repository for something fully derivable from
support_data/merged_full_placed.nwk (56 MB), so it is generated here instead.

    python3 tools/make_merged_tree_npz.py

writes support_data/mergedtree_dttpq_v8.npz with:

    parent   int32   (n_nodes,)  index of each node's parent; -1 at the root
    blen     float64 (n_nodes,)  branch length to the parent
    time     float64 (n_nodes,)  root-to-node distance (the dated tree's node time)
    is_leaf  bool    (n_nodes,)  True for tips
    names    object  (n_nodes,)  node label; "" for unlabelled internal nodes

Nodes are in the order the Newick is parsed, which is the order the analysis
scripts assume: a node always appears after its parent.
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DEFAULT_IN = os.path.join(REPO, "support_data", "merged_full_placed.nwk")
DEFAULT_OUT = os.path.join(REPO, "support_data", "mergedtree_dttpq_v8.npz")


def parse_newick(text):
    """Iterative Newick parser. Returns (parent, blen, is_leaf, names).

    Iterative rather than recursive because the tree is ~2.4M nodes deep enough
    to blow the default recursion limit.
    """
    text = text.strip()
    if text.endswith(";"):
        text = text[:-1]

    parent, blen, is_leaf, names = [], [], [], []
    stack = []
    i, n = 0, len(text)

    def add(par):
        parent.append(par)
        blen.append(0.0)
        is_leaf.append(True)
        names.append("")
        return len(parent) - 1

    def read_label(j):
        k = j
        while k < n and text[k] not in "(),:;":
            k += 1
        return text[j:k], k

    def read_len(j):
        if j < n and text[j] == ":":
            k = j + 1
            while k < n and (text[k].isdigit() or text[k] in ".eE+-"):
                k += 1
            try:
                return float(text[j + 1:k]), k
            except ValueError:
                return 0.0, k
        return 0.0, j

    cur = None
    while i < n:
        c = text[i]
        if c == "(":
            node = add(stack[-1] if stack else -1)
            if stack:
                is_leaf[stack[-1]] = False
            stack.append(node)
            cur = None
            i += 1
        elif c == ",":
            i += 1
            cur = None
        elif c == ")":
            cur = stack.pop()
            i += 1
            lab, i = read_label(i)
            if lab:
                names[cur] = lab
            bl, i = read_len(i)
            blen[cur] = bl
        else:
            lab, i = read_label(i)
            bl, i = read_len(i)
            if lab or bl:
                node = add(stack[-1] if stack else -1)
                if stack:
                    is_leaf[stack[-1]] = False
                names[node] = lab
                blen[node] = bl
                cur = node
            else:
                i += 1
    return (np.asarray(parent, dtype=np.int32),
            np.asarray(blen, dtype=np.float64),
            np.asarray(is_leaf, dtype=bool),
            np.asarray(names, dtype=object))


def node_times(parent, blen):
    """Root-to-node distance. Relies on parents preceding children."""
    t = np.zeros(parent.shape[0], dtype=np.float64)
    for i in range(parent.shape[0]):
        p = parent[i]
        if p >= 0:
            t[i] = t[p] + blen[i]
    return t


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="inp", default=DEFAULT_IN)
    ap.add_argument("--out", dest="out", default=DEFAULT_OUT)
    args = ap.parse_args()

    if not os.path.exists(args.inp):
        sys.exit(f"missing input tree: {args.inp}")

    sys.stderr.write(f"reading {args.inp}\n")
    parent, blen, is_leaf, names = parse_newick(open(args.inp).read())
    sys.stderr.write(f"  {parent.size:,} nodes  ({int(is_leaf.sum()):,} leaves)\n")
    time = node_times(parent, blen)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    np.savez_compressed(args.out, parent=parent, blen=blen, time=time,
                        is_leaf=is_leaf, names=names)
    sys.stderr.write(f"wrote {args.out} "
                     f"({os.path.getsize(args.out) / 1e6:.1f} MB)\n")


if __name__ == "__main__":
    main()
