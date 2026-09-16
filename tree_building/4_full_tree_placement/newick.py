#!/usr/bin/env python3
"""Minimal, dependency-free Newick I/O tuned for these trees.

These lineage trees are huge (hundreds of thousands of tips) and extremely
unbalanced ("caterpillars"), so recursion would overflow the stack. Every parser
here is iterative and single-pass. No ete3 / dendropy.

Two readers, for two different needs:
  parse_structure(text) -> full node arrays (parent, name, kids, blen), used by the
      tree BUILDER which must expand backbone leaves in place.
  iter_sibling_pairs(path) / tips_of(path) -> the leaf/sibling view.
"""


def parse_structure(text):
    """Parse a Newick string into flat node arrays.

    Returns (parent, name, kids, blen) where each is indexed by node id (0 = root):
      parent[i]  parent node id, -1 for the root
      name[i]    label string ("" for unnamed internal nodes)
      kids[i]    list of child node ids
      blen[i]    branch length to parent as a string ("" if none) — kept verbatim so
                 we can re-emit backbone lengths without float round-trips.
    """
    parent = [-1]
    name = [""]
    kids = [[]]
    blen = [""]
    cur = 0
    i = 0
    n = len(text)
    stack = []

    def new_node(par):
        j = len(parent)
        parent.append(par)
        name.append("")
        kids.append([])
        blen.append("")
        kids[par].append(j)
        return j

    while i < n:
        c = text[i]
        if c == '(':
            cur = new_node(cur)
            stack.append(parent[cur])
            i += 1
        elif c == ',':
            cur = new_node(stack[-1])
            i += 1
        elif c == ')':
            cur = stack.pop()
            i += 1
        elif c == ';':
            break
        elif c == ':':
            # branch length token: read to the next structural char, store verbatim
            j = i + 1
            while j < n and text[j] not in ',():;':
                j += 1
            blen[cur] = text[i + 1:j]
            i = j
        else:
            # a label token (leaf name or internal support/name)
            j = i
            while j < n and text[j] not in ',():;':
                j += 1
            name[cur] = text[i:j]
            i = j
    return parent, name, kids, blen


def iter_sibling_pairs(path):
    """Yield structure for the leaf/sibling view in one pass.

    Returns (leaves, sib_pairs, n_internal, n_cherries):
      leaves      list of leaf labels (order = first appearance)
      sib_pairs   list of (a, b) leaf pairs sharing a parent (a polytomy of k leaves
                  contributes C(k, 2) pairs; a bifurcating node contributes its cherry)
      n_internal  count of internal nodes
      n_cherries  count of strict cherries (exactly two children, both leaves)
    Subtrees are collapsed to a marker as they close, so memory is O(depth), not O(n).
    """
    with open(path) as fh:
        s = fh.read()

    INTERNAL = object()
    stack = [[]]
    leaves = []
    sib_pairs = []
    n_internal = 0
    n_cherries = 0
    tok = []
    just_closed = False

    def take():
        if not tok:
            return ""
        nm = "".join(tok).split(":", 1)[0].strip()
        tok.clear()
        return nm

    for ch in s:
        if ch == '(':
            take()
            stack.append([])
            tok.clear()
            just_closed = False
        elif ch in ',);':
            if just_closed:
                take()                       # discard closed node's own label
                just_closed = False
            else:
                nm = take()
                if nm:
                    stack[-1].append(nm)
                    leaves.append(nm)
            if ch == ')':
                children = stack.pop()
                n_internal += 1
                leaf_kids = [c for c in children if c is not INTERNAL]
                if len(children) == 2 and len(leaf_kids) == 2:
                    n_cherries += 1
                for a in range(len(leaf_kids)):
                    for b in range(a + 1, len(leaf_kids)):
                        sib_pairs.append((leaf_kids[a], leaf_kids[b]))
                stack[-1].append(INTERNAL)
                just_closed = True
            elif ch == ';':
                break
        else:
            tok.append(ch)

    return leaves, sib_pairs, n_internal, n_cherries


def tips_of(path):
    """Set of leaf labels in a Newick file."""
    leaves, _, _, _ = iter_sibling_pairs(path)
    return set(leaves)
