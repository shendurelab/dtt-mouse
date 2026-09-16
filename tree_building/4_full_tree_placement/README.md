# Placing the remaining cells onto the dated backbone

The NJ backbone (`1_build_nj_backbone`) is built only from the `ge7_founderok`
subset (~640k cells). This step attaches every other cell in the full per-side
consensus (~1.3M cells total) by DTT-ancestor nearest-neighbour search: each cell
attaches to its single closest cell (a backbone anchor or another placed cell),
pendant length = its DTT edit distance to their inferred common ancestor, scaled
to days by that side's LSD2 clock rate.

Run per side:

```
./run_side.sh B1
./run_side.sh B2
```

Then merge both sides -- `dated_placed_{SIDE}.nwk` is a dated tree just like
`3_date_tree`'s output, so the same merge script applies unchanged:

```
LEFT_NWK=results/4-full-tree-placement/B1/dated_placed_B1.nwk \
RIGHT_NWK=results/4-full-tree-placement/B2/dated_placed_B2.nwk \
OUT_NWK=results/4-full-tree-placement/merged/merged_placed.nwk \
  Rscript 3_date_tree/merge_dated_subtrees.R
```

## How it works

1. `bestmatch.py` -- encode genotypes, build a rarity-capped inverted index over
   all cells, then for every non-backbone cell find its DTT-closest cell overall
   and its DTT-closest backbone cell (`dtt_match.py`).
2. `finalize_dated.py` -- resolve each cell's attachment point (`build_tree.py`,
   breaking any mutual-nearest-neighbour cycles so every chain reaches the
   backbone), compute its pendant length in DTT units (`dtt_lengths.py`), convert
   to days with the clock rate, and graft it onto the dated backbone
   (`dated_tree.py`). Also flags attachments that violate edit irreversibility
   (`validity.py`, diagnostic only -- flagged, not rejected).

Not every cell is guaranteed a place: a cell with no shared edit token with
anything, or one whose nearest-neighbour chain runs into an ungroundable cycle,
is dropped (see `build_tree.py`'s `resolve_parents` docstring). Check the run's
`held N,NNN` summary line to see how many.

## Tests

```
python3 -m pytest tests/
```

`tests/test_pendant_matches_R.py` cross-checks the pendant computation against
`1_build_nj_backbone/dtt_distance.R` on real data; it skips if `Rscript` isn't
on PATH.

## Not included here

- The cosine-similarity placement variant -- an earlier method, superseded by
  the DTT-ancestor search used here (validated to change topology, not just
  branch lengths: only ~12% of attachments agree between the two).
- Extending the placement to non-`pass_qc` cells (a second "freeze-and-add"
  round) -- not part of the tree this repo's `CURRENT_TREE.txt`-equivalent
  deliverable points at.
