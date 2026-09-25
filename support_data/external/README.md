# External inputs

Several inputs the analysis uses are not redistributed here — either because they
belong to other publications or archives, or because of size. Fetch them into this
directory (or point the listed environment variable somewhere else) before running the
scripts that need them.

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

## 3. Embryo #3 `.h5ad` object

`Figure-S2AB_S6.ipynb` reads the embryo #3 object via `DTT_ADATA_EMBRYO3` (default:
`adata_embryo3.h5ad` in the working directory). It's hosted separately, already under
that name:

<https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/tapemouse/adata_embryo3.h5ad>

Save it as `making_figures/adata_embryo3.h5ad`, or point the environment variable at
wherever you put it:

```bash
export DTT_ADATA_EMBRYO3=/path/to/adata_embryo3.h5ad
```

## 4. Raw tape calls

`tape_calls.tsv.gz` (~250 MB) is an intermediate of the single-cell tape pipeline, produced
from the raw reads. Raw sequencing data are on GEO under **GSE341627**; regenerate the
intermediate by running `sc_tape/run_sc.sh` with `TAPE_RAW_DIR` pointing at the downloaded
reads. Scripts that need it read `DTT_TAPE_CALLS`:

```bash
export TAPE_RAW_DIR=/path/to/raw
export DTT_TAPE_CALLS=/path/to/tape_calls.tsv.gz
```

---

## 5. Ancestral-state node inference

`pd_nodes_infer_200.txt` (~145 MB) is an output of the ancestral-state step, consumed by
`tree_analysis/traceback_paths/02_tableS8_qualifying_paths.py`. Too large to ship;
regenerate it with the `tree_analysis/ancestral_state/` scripts, or point `DTT_PD_NODES`
at an existing copy.

---

# The tree is read straight from Newick

Earlier versions of these scripts read the tree from a ~93 MB `mergedtree_dttpq_v8.npz`
bundle. That file was a flat array re-encoding of `support_data/merged_full_placed.nwk`,
which already ships, so nothing needs to be generated or committed: `tools/tree_io.py`
parses the Newick directly into the same five arrays (`parent`, `blen`, `time`, `is_leaf`,
`names`).

Parsing 2.4M nodes takes ~8 s, so the result is cached under
`support_data/.tree_cache/` (gitignored); repeat loads take ~0.4 s. Delete that directory
any time and it rebuilds. Point `DTT_TREE` at a different Newick to use another tree.

Verified array-by-array against the original bundle: all five arrays identical across
2,391,035 nodes / 1,281,141 leaves.

# Environment variables

| Variable | Purpose |
| --- | --- |
| `DTT_CELL_METADATA` | override `support_data/cell_metadata.v8.txt.gz` |
| `DTT_ROUTING_LABELS` | override `support_data/e3v8.routing_labels.tsv.gz` |
| `DTT_TREE` | override the tree Newick (default `support_data/merged_full_placed.nwk`) |
| `DTT_PD_NODES` | location of `pd_nodes_infer_200.txt` |
| `DTT_RESULTS` | where analysis steps write intermediates (default `<step>/out/`) |
| `DTT_NEWCODE_NPZ` | override `support_data/newcode_v8_all.npz` |
| `QIU2024_SUPP_XLSX` | location of the Qiu et al. 2024 supplementary table |
| `DTT_TAPE_CALLS` | location of `tape_calls.tsv.gz` |
| `DTT_ADATA_EMBRYO3` | location of the embryo-3 `.h5ad` |
| `TAPE_RAW_DIR` | root of the raw reads downloaded from GEO |
