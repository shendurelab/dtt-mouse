"""The edit-chain trie -- the lineage abstraction shared by bulk and single-cell.

Because TAPE editing is strictly ordered and monotonic, a lineage's genotype at
one integration is fully described by the ORDERED CHAIN of edit symbols from the
5' end up to the first unedited site. A shallower chain is always a PREFIX of a
deeper one, so all chains at a barcode nest into a trie whose root->tip paths ARE
the lineage relationships.

Two consumers, one structure:
  - single cell: build a trie per (cell, integration) from UMI/molecule counts,
    then `consensus()` collapses carry-over prefixes onto the deepest dominant
    tip and flags a competing allele only at a real branch -> ONE genotype/cell.
  - bulk: build a trie per integration from read counts; the trie itself is the
    per-integration lineage tree. `fold_dropout()` optionally reassigns shallow
    interior tips that look like 3'-truncation onto their deeper parent lineage.
"""
from __future__ import annotations

from config import DROPOUT_FOLD_RATIO, BRANCH_MIN_FRAC, BRANCH_MIN_READS


def edit_chain(pattern) -> tuple:
    """Ordered edit symbols of a genotype, up to the first unedited/unresolved site.

    Accepts either a full site pattern (e.g. ('GCC','GCC','U','U','U','U')) or an
    already-trimmed chain (('GCC','GCC',...)). Trailing 'U'/None are dropped.
    """
    chain = []
    for s in pattern:
        if s == "U" or s is None:
            break
        chain.append(s)
    return tuple(chain)


class TrieNode:
    __slots__ = ("symbol", "depth", "parent", "children", "tip", "subtree", "y")

    def __init__(self, symbol=None, depth=0, parent=None):
        self.symbol = symbol        # the edit symbol on the edge into this node
        self.depth = depth          # number of edits from the root (0 = root)
        self.parent = parent
        self.children = {}          # symbol -> TrieNode
        self.tip = 0.0              # weight of chains ENDING exactly here
        self.subtree = 0.0          # weight of this node's whole subtree (tip + all descendants)
        self.y = 0.0                # layout coordinate (set by the tree renderer)

    def child(self, symbol) -> "TrieNode":
        c = self.children.get(symbol)
        if c is None:
            c = TrieNode(symbol, self.depth + 1, self)
            self.children[symbol] = c
        return c

    @property
    def through(self) -> float:
        """Weight that passes THROUGH this node into deeper lineages."""
        return self.subtree - self.tip

    def chain(self) -> tuple:
        """The edit chain (root -> this node)."""
        out, n = [], self
        while n.parent is not None:
            out.append(n.symbol)
            n = n.parent
        return tuple(reversed(out))


def build_trie(counts: dict) -> TrieNode:
    """Build an edit-chain trie from {genotype: weight}. `weight` is reads (bulk)
    or molecules (single cell). Genotypes may be full patterns or chains."""
    root = TrieNode()
    for pat, w in counts.items():
        node = root
        for sym in edit_chain(pat):
            node = node.child(sym)
        node.tip += w
    _compute_subtree(root)
    return root


def _compute_subtree(node: TrieNode) -> float:
    total = node.tip
    for c in node.children.values():
        total += _compute_subtree(c)
    node.subtree = total
    return total


def is_branch(node: TrieNode,
              min_frac: float = BRANCH_MIN_FRAC,
              min_reads: float = BRANCH_MIN_READS) -> bool:
    """True if >=2 children are each well-supported -- two different edit symbols
    competing at the same next site (a doublet/mixture, or a real lineage split)."""
    strong = [c for c in node.children.values()
              if c.subtree >= min_reads and c.subtree >= min_frac * max(node.through, 1e-9)]
    return len(strong) >= 2


def consensus(root: TrieNode,
              min_frac: float = BRANCH_MIN_FRAC,
              min_reads: float = BRANCH_MIN_READS):
    """Single-cell consensus: descend the dominant child, folding every prefix's
    carry-over onto the chain, until a leaf or a real branch.

    Returns (chain, branched, dominance):
      chain     -- consensus edit chain (deepest dominant tip)
      branched  -- True if a competing allele (branch) was hit (flag as doublet)
      dominance -- weight on the chosen path / total weight below the split point
    """
    node = root
    while node.children:
        if is_branch(node, min_frac, min_reads):
            return node.chain(), True, _dominance(node)
        top = max(node.children.values(), key=lambda c: c.subtree)
        node = top
    return node.chain(), False, _dominance(node.parent) if node.parent else 1.0


def _dominance(node: TrieNode) -> float:
    if node is None or not node.children:
        return 1.0
    kids = sorted(node.children.values(), key=lambda c: -c.subtree)
    denom = sum(c.subtree for c in kids) + node.tip
    return kids[0].subtree / denom if denom else 1.0


def dominant_lineage(root: TrieNode):
    """Descend by max subtree to the deepest DOMINANT tip and score it against the
    WHOLE locus. Returns (chain, n_support, dominance, branched):

      n_support -- molecules on the dominant root->tip path = the dominant genotype's
                   own UMIs PLUS every prefix (shorter chain) on that path, since a
                   prefix is 3' dropout of the same lineage (= total - off-path).
      dominance -- n_support / total_molecules. Genotypes that BRANCH off the path (a
                   different edit at some site = a second lineage) are the only thing
                   that lowers it; dropout prefixes do not.
      branched  -- True if any off-path branch carried >= 2 molecules (a real competitor).

    Used by the single-cell caller: call the dominant genotype only when it is both
    well-supported and clearly dominant, else leave the locus missing.
    """
    total = root.subtree
    if total <= 0:
        return (), 0.0, 0.0, False
    off = 0.0
    branched = False
    node = root
    while node.children:
        top = max(node.children.values(), key=lambda c: c.subtree)
        for c in node.children.values():
            if c is not top:
                off += c.subtree
                if c.subtree >= 2:
                    branched = True
        node = top
    n_support = total - off
    return node.chain(), n_support, n_support / total, branched


def fold_dropout(root: TrieNode, ratio: float = DROPOUT_FOLD_RATIO) -> None:
    """Bulk dropout folding (in place). An interior node's tip weight is treated
    as 3'-truncation of its dominant deeper child -- not a genuine early stop --
    and reassigned to that child's tip, but ONLY when the child's subtree
    outweighs the interior tip by >= `ratio`. Conservative by design: a shallow
    lineage with its own substantial read support survives as a real stop.

    ratio <= 0 disables folding (every observed depth kept as a distinct tip).
    """
    if ratio <= 0:
        return
    # bottom-up so folds cascade correctly
    for node in _postorder(root):
        if node.tip <= 0 or not node.children:
            continue
        top = max(node.children.values(), key=lambda c: c.subtree)
        if top.subtree >= ratio * node.tip:
            # attribute truncated reads to the DEEPEST dominant descendant lineage
            dest = top
            while dest.children:
                dest = max(dest.children.values(), key=lambda c: c.subtree)
            dest.tip += node.tip
            node.tip = 0.0
    _compute_subtree(root)


def _postorder(root: TrieNode):
    out = []
    def rec(n):
        for c in n.children.values():
            rec(c)
        out.append(n)
    rec(root)
    return out


def tips(root: TrieNode):
    """Yield (chain, weight) for every surviving lineage tip (tip weight > 0)."""
    for node in _postorder(root):
        if node.tip > 0:
            yield node.chain(), node.tip


def nodes_per_level(root: TrieNode, n_levels: int) -> list[int]:
    """Count distinct lineage states (trie nodes) at each depth 0..n_levels."""
    cnt = [0] * (n_levels + 1)
    for node in _postorder(root):
        if node.depth <= n_levels:
            cnt[node.depth] += 1
    return cnt
