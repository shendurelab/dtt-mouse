# External inputs

Three inputs the analysis uses are not redistributed here — two because they belong to
other publications or archives, one because of size. Fetch them into this directory (or
point the listed environment variable somewhere else) before running the scripts that
need them.

## 1. Qiu et al. 2024 mouse developmental atlas, Supplementary Table 4

Used by the clade-k / timed-fate-coupling analyses as the reference cell-type hierarchy.

- Paper: Qiu, C. *et al.* "A single-cell time-lapse of mouse prenatal development from gastrula
  to birth." *Nature* **626**, 1084–1093 (2024). <https://doi.org/10.1038/s41586-024-07069-w>
- File: Supplementary Table 4, `41586_2024_7069_MOESM4_ESM.xlsx`, from the paper's
  Supplementary Information.

Save it here as `41586_2024_7069_MOESM4_ESM.xlsx`, or set:

```bash
export QIU2024_SUPP_XLSX=/path/to/41586_2024_7069_MOESM4_ESM.xlsx
```

Not redistributed: it is Springer Nature supplementary material, not ours to relicense.

## 2. Mouse atlas `.h5ad` objects

Used by the cell-type integration and ancestral-state imputation steps
(`scRNAseq_processing/step4_integration_atlas_E1275_E1425.R`,
`tree_analysis/step5_*`). Downloaded per developmental day from:

<https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax/mm39_version/>

as `adata.{day_id}.h5ad`. The scripts already carry this URL in a comment at the top.

## 3. Raw tape calls

`tape_calls.tsv.gz` (~250 MB) is an intermediate of the single-cell tape pipeline, produced
from the raw reads. Raw sequencing data are on GEO under **GSE341627**; regenerate the
intermediate by running `sc_tape/run_sc.sh` with `TAPE_RAW_DIR` pointing at the downloaded
reads. Scripts that need it read `DTT_TAPE_CALLS`:

```bash
export TAPE_RAW_DIR=/path/to/raw
export DTT_TAPE_CALLS=/path/to/tape_calls.tsv.gz
```

---

# Generated, not shipped

`support_data/mergedtree_dttpq_v8.npz` is a flat array encoding of the merged, placed,
dated tree that several `tree_analysis` scripts read instead of re-parsing Newick. It is
~93 MB and fully derivable from `support_data/merged_full_placed.nwk`, which already ships,
so it is generated rather than committed:

```bash
python3 tools/make_merged_tree_npz.py
```

Takes a few seconds and reproduces the original bundle exactly (verified array-by-array:
`parent`, `blen`, `time`, `is_leaf`, `names` over all 2,391,035 nodes / 1,281,141 leaves).
Override the location with `DTT_MERGED_TREE_NPZ`.

# Environment variables

| Variable | Purpose |
| --- | --- |
| `DTT_CELL_METADATA` | override `support_data/cell_metadata.v8.txt.gz` |
| `DTT_ROUTING_LABELS` | override `support_data/e3v8.routing_labels.tsv.gz` |
| `DTT_MERGED_TREE_NPZ` | override the generated merged-tree bundle |
| `DTT_NEWCODE_NPZ` | override `support_data/newcode_v8_all.npz` |
| `QIU2024_SUPP_XLSX` | location of the Qiu et al. 2024 supplementary table |
| `DTT_TAPE_CALLS` | location of `tape_calls.tsv.gz` |
| `TAPE_RAW_DIR` | root of the raw reads downloaded from GEO |
