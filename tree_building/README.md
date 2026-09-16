# Building the  NJ backbone and full tree

Four steps, run in order, from this directory as the repo root:

1. `1_build_nj_backbone/` — build the per-side (B1/B2) distance matrix and NJ trees
2. `2_root_tree/` — root each tree and drop the synthetic root tip
3. `3_date_tree/` — date the rooted tree with LSD2, under a literature cell-count constraint
4. `4_full_tree_placement/` — place the remaining cells onto the dated backbone

## 1. Build the distance matrix + NJ tree

`run_b1.sh` / `run_b2.sh` are the exact commands used for the real build.
Each points `run_full_distance.R` at that side's `ge7_founderok.tsv.gz` (see
`processed_data/README.md`), which parses the tapes, builds the distance
matrix, appends a synthetic founder-genotype root for the next step, and
hands the matrix straight to DecentTree's RapidNJ in-process due to large memory
and a lot of time spent writing and reading the distance matrix.
file, no re-parse.

**Needs a big-memory machine**: B1 (401,902 tips) ~3.0 TB RAM, B2 (238,112 tips)
~1.06 TB RAM. Not runnable on a normal workstation.

## 2. Root the tree

`root_tree.R` roots the raw NJ tree on its synthetic-root tip and drops it --
self-contained, just `ape::root`/`drop.tip`.

## 3. Date the tree

`3_date_tree/` corrects negative NJ branch lengths, dates the tree with LSD2,
and re-dates under a literature cell-count constraint so early splits can't
imply more coexisting lineages than the embryo actually had cells. See
`3_date_tree/README.md`.

## 4. Place the remaining cells

`4_full_tree_placement/` attaches every cell with at least 4 tapes that is not already in the NJ backbone onto the dated backbone by nearest-neighbour DTT distance, grafting each with a time-scaled pendant. See
`4_full_tree_placement/README.md`.

## Tests

```
Rscript 1_build_nj_backbone/tests/test_full_distance_pipeline.R
python3 -m pytest 4_full_tree_placement/tests/
```

## External prerequisites

- R + `Rcpp`, `RcppParallel`, `data.table`, `ape`, `BAT`
- Python 3 + `numpy`
- DecentTree source, for the in-process NJ build:
  `git clone --recursive https://github.com/iqtree/decenttree.git decenttree`
- LSD2 source, for tree dating:
  `git clone https://github.com/tothuhien/lsd2.git lsd2 && (cd lsd2/src && make)`
- `pigz` on PATH

