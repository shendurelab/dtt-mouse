# Dating the rooted NJ tree

Two steps, run in order, per side (B1/B2) then merged:

1. `run_dating.sh` — BAT-correct negative NJ branches, LSD2-date, merge B1+B2
2. `run_constrained_dating.sh` — re-date under a literature min-age ladder, so
   early splits can't imply more coexisting lineages than the embryo had cells

Both call `run_lsd2.sh` for the actual LSD2 invocation (and the day-unit
newick it emits via `time_tree.R`); constrained re-dating just adds a
minimum-age datefile (`build_min_age_datefile.R`'s output) on top.

## 1. Unconstrained dating

```
ROOTED_TREE_B1=results/2-rooted-nj/B1/nj_rooted_ingroup.nwk \
ROOTED_TREE_B2=results/2-rooted-nj/B2/nj_rooted_ingroup.nwk \
STAGEDIR=results/3-dated-tree MERGE_TOKEN=ge7 \
  bash 3_date_tree/run_dating.sh
```

Takes each side's rooted tree from `2_root_tree/`, corrects negative branch
lengths (an NJ artefact -- LSD2 needs non-negative substitution counts), runs
LSD2 with root/tip calibration dates, then merges B1+B2 on a shared day-0
zygote root.

## 2. Constrained re-dating

Plain LSD2 reads "zero edits" as "~zero elapsed time" and can pack more early
splits into the first few days than the embryo actually had cells. This step
builds an LSD2 minimum-age constraint per node from a literature cell-count
ceiling (`processed_data/sample_matched_ceiling_sourced.csv`) and re-dates,
iterated 3x (each pass re-ranks nodes on the previous pass's output):

```
STAGEDIR=results/3-dated-tree MERGE_TOKEN=ge7 \
RANK_TREE_B1=results/3-dated-tree/B1/lsd2/nj99478/minB2h/time_tree.nwk \
RANK_TREE_B2=results/3-dated-tree/B2/lsd2/nj99478/minB2h/time_tree.nwk \
  bash 3_date_tree/run_constrained_dating.sh
```

This is the tree used downstream: `results/3-dated-tree/merged/merged_time_tree_ge7_attempt1_minage_sourced_minB2h_l0.01.nwk`.

## 3. Sensitivity analysis: asymmetric B1/B2 split

Step 2 assumes each side's cell-count ceiling is half the whole-embryo
ceiling (`SIDE_FRAC_B1`/`SIDE_FRAC_B2`, both default `0.5`). To test
sensitivity to that assumption -- e.g. an asymmetric 0.63/0.37 split -- rerun
step 2 with those env vars overridden, staged under a separate `STAGEDIR` so
the primary run's outputs are untouched:

```
# a. stage the existing divergence trees (2_root_tree's output) into the
#    STAGEDIR layout run_constrained_dating.sh expects
mkdir -p results/3-dated-tree-sensitivity/B1/lsd2/nj99478
mkdir -p results/3-dated-tree-sensitivity/B2/lsd2/nj99478
cp results/2-rooted-nj/divergence_nonneg_B1.nwk results/3-dated-tree-sensitivity/B1/lsd2/nj99478/tree_nonneg.nwk
cp results/2-rooted-nj/divergence_nonneg_B2.nwk results/3-dated-tree-sensitivity/B2/lsd2/nj99478/tree_nonneg.nwk

# b. run the constrained re-dating with the asymmetric split
STAGEDIR=results/3-dated-tree-sensitivity MERGE_TOKEN=ge7 \
RANK_TREE_B1=results/3-dated-tree/perside_B1_minB2h_unconstrained.nwk \
RANK_TREE_B2=results/3-dated-tree/perside_B2_minB2h_unconstrained.nwk \
VARIANT=attempt1_minage_sourced_minB2h_l0.01_sidefrac63-37 \
SIDE_FRAC_B1=0.63 SIDE_FRAC_B2=0.37 \
  bash 3_date_tree/run_constrained_dating.sh

# c. copy the final outputs to flat, self-describing names (mirrors the
#    perside_*/merged_* naming already used under results/3-dated-tree/)
cp results/3-dated-tree-sensitivity/B1/lsd2/nj99478/attempt1_minage_sourced_minB2h_l0.01_sidefrac63-37/time_tree.nwk \
   results/3-dated-tree-sensitivity/perside_B1_minB2h_lineage_constrained_sidefrac63-37.nwk
cp results/3-dated-tree-sensitivity/B2/lsd2/nj99478/attempt1_minage_sourced_minB2h_l0.01_sidefrac63-37/time_tree.nwk \
   results/3-dated-tree-sensitivity/perside_B2_minB2h_lineage_constrained_sidefrac63-37.nwk
cp results/3-dated-tree-sensitivity/merged/merged_time_tree_ge7_attempt1_minage_sourced_minB2h_l0.01_sidefrac63-37.nwk \
   results/3-dated-tree-sensitivity/merged_minB2h_lineage_constrained_sidefrac63-37.nwk
```

This reruns LSD2 constrained dating 3 iterations x 2 sides on trees with
~400k (B1) and ~240k (B2) tips -- budget time/resources accordingly.

## External prerequisites

- R + `ape`, `BAT`
- LSD2 source, built from `src/`:
  `git clone https://github.com/tothuhien/lsd2.git lsd2 && (cd lsd2/src && make)`
