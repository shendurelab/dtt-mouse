# DTT Mouse analysis

Code and processed data for *A DNA Typewriter records the cell lineage history of a mouse, from zygote to late organogenesis*.

## Abstract

Mammalian biology unfolds over time, within tissues and organs opaque to our eyes and instruments. We applied DNA Typewriter, a sequential molecular recorder, to record the cell lineage of a mouse over nearly two weeks of development. From one embryo, we reconstruct a time-calibrated, parsimony-supported, zygote-rooted phylogeny of 1.28M transcriptionally-profiled cells. A burst of editing unequivocally marks the daughters of the first cleavage, which serve as inline replicates. We quantify clonal dominance arising during gastrulation. Tree siblings share cell type far above chance; heterotypic siblings mark terminal differentiations. Temporal sweeps of clade co-occurrence recover a dated hierarchy of cell-type couplings; imputed labels for internal nodes recapitulate known state paths. A lineage-anchored ontogeny of mouse development, long out of reach, is coming into view.

---

## Where to start

Most readers will not want to re-run the pipeline. The three most useful entry points:

| If you want to… | Go to |
| --- | --- |
| Look at the tree | `support_data/merged_full_placed.nwk` (1,281,141 tips) or `merged_minB2h_lineage_constrained.nwk` (the 655,701-tip dated backbone) |
| Reproduce a figure | `making_figures/Figure-<N>.R`, reading from `figures_data/` |
| Understand how a Methods sentence was computed | `bulk_tape/METHODS.md` and `sc_tape/METHODS.md`, which map Methods prose to specific functions |

Raw sequencing data are on GEO under **GSE341627**. Cell metadata, tape-consensus
matrices, routing labels, founder genotypes and the trees themselves ship here in
`support_data/`, so the analysis and figure code runs without the raw data.

---

## Repository map

The pipeline runs left to right: raw reads become tape genotypes, tape genotypes become a
tree, and the tree is analysed and plotted.

```
raw fastq ──┬─> bulk_tape/          bulk amplicon TAPE
            ├─> sc_tape/            single-cell circtape ──┐
            └─> scRNAseq_processing/ transcriptome  ───────┤
                                                           v
                                              tree_building/  (NJ → root → date → place)
                                                           v
                                              tree_analysis/  (the biology)
                                                           v
                                              making_figures/ (panels)
```

| Directory | What it holds |
| --- | --- |
| **`bulk_tape/`** | Bulk amplicon TAPE: reads → integration barcodes and insertion vocabulary → parsing → de-noising → copy number → per-integration lineage trees. Self-contained, with `METHODS.md`, `config.py`, `example_data/` and `tests/`. |
| **`sc_tape/`** | Single-cell circtape: read merging → parsing and de-noising → per-cell consensus → blastomere routing → founder-state filter, ending in the two per-blastomere matrices that feed tree building. Same layout as `bulk_tape/`. |
| **`scRNAseq_processing/`** | sci-RNA-seq3 transcriptome processing: alignment and gene counts, RT-barcode correction, doublet removal, clustering and cell-type annotation by transfer from the mouse developmental atlas. |
| **`tree_building/`** | Four ordered stages — `1_build_nj_backbone/`, `2_root_tree/`, `3_date_tree/`, `4_full_tree_placement/` — plus `lib/` and the inputs in `processed_data/`. |
| **`tree_analysis/`** | Everything computed on the finished tree: editing rate, sibling analyses, clade co-occurrence, clonal dominance, ancestral-state imputation, timed fate couplings, traceback paths, dating QC, held-out placement accuracy. |
| **`making_figures/`** | One script per figure panel set, reading from `figures_data/`. |
| **`figures_data/`** | Tidy CSVs, one or more per panel, so figures can be redrawn without re-running any analysis. |
| **`support_data/`** | The shared inputs: trees, cell metadata, tape-consensus matrices, routing labels and reports, founder genotypes, insertion sites, cell-type tables. |

---

## Running it

### Environment

- **Python 3.9+** (developed on 3.9.6 and 3.12) — `numpy`, `scipy`, `matplotlib`; see
  `bulk_tape/requirements.txt` and `sc_tape/requirements.txt`
- **R** — `Rcpp`, `RcppParallel`, `data.table`, `ape`, `BAT`, plus `ggplot2` and friends
  for the figure scripts
- For tree building only, two external sources:
  - [DecentTree](https://github.com/iqtree/decenttree) for the in-process NJ build
  - [LSD2](https://github.com/tothuhien/lsd2) for molecular-clock dating
  - `pigz` on `PATH`

### Tape processing

```bash
# bulk, one embryo, stages 01–06
cd bulk_tape && ./run_bulk.sh DTTz_3_S3

# single cell, embryo #3, raw fastqs → per-blastomere matrices
cd sc_tape && TAPE_RAW_DIR=/path/to/fastq ./run_sc.sh
```

Both need the raw data from GEO. Every stage after the first is deterministic and
re-runnable, and both packages ship tests that run in seconds against `example_data/`
without any raw data:

```bash
python3 -m pytest bulk_tape/tests/ sc_tape/tests/
```

### Tree building

Run the four stages in order from `tree_building/`; see `tree_building/README.md` for
details and `3_date_tree/README.md` and `4_full_tree_placement/README.md` for the last two.

**This stage needs a big-memory machine.** Neighbour-joining peaks at roughly `20·n²`
bytes: about 3.0 TB for blastomere A (410,925 tips) and 1.06 TB for blastomere B
(244,776 tips). It is not runnable on a workstation. The resulting trees are checked in
under `tree_building/results/` and copied to `support_data/`, so everything downstream
works without repeating this step.

### Analysis and figures

`tree_analysis/` scripts read the trees from `support_data/` and are grouped by analysis;
the numbered `stepN_*.R` scripts at the top level run in order, and the subdirectories hold
the per-analysis code (`clade_cooccurrence/`, `clonal_dominance/`, `ancestral_state/`,
`timed_fate_couplings/`, `traceback_paths/`, `placement_accuracy/`, `dating_qc/`, …).

Figure scripts read tidy CSVs from `figures_data/` and can each be run on their own. They
use repo-relative paths like `./figures_data/…`, so **run them from the repository root**:

```bash
Rscript making_figures/Figure-1.R
```

---

## Figures

| Figure | Script |
| --- | --- |
| 1 | `making_figures/Figure-1.R` |
| 2 | `making_figures/Figure-2.R`; panels F–H via `Figure-2fgh.py` and `Figure-2fgh_select_zoom.R` |
| 4 | `making_figures/Figure-4.R` |
| 5 | `making_figures/Figure-5.R` |
| 6 | `making_figures/Figure-6.R` |
| S1 | `making_figures/Figure-S1.R` |
| S2, S6 | `making_figures/Figure-S2.R`, `Figure-S2AB_S6 .ipynb` |
| S3 | `making_figures/Figure-S3.R` |
| S4 | `making_figures/Figure-S4.R` |
| S5 | `making_figures/Figure-S5.R` |
| S7 | `making_figures/Figure-S7.R` |
| S8 | `making_figures/Figure-S8.R`; panel B from `tree_analysis/heldout_tape/fig_s8b_heldout_tape_distance_correlation.R` |
| S9 | `making_figures/Figure-S9.R` |
| S10 | `making_figures/Figure-S10.R` |
| S11 | `making_figures/Figure-S11.R` |
| S12 | `making_figures/Figure-S12.R`, scored by `tree_analysis/placement_accuracy/fig_s12_placement_accuracy_analysis.R` |
| S13 | `making_figures/Figure-S13.R` |
| S14 | `making_figures/Figure-S14.R` |
| S15 | `making_figures/Figure-S15.R` |
| S17 | `making_figures/Figure-S17.R` |
| S18 | `making_figures/Figure-S18.R` |

`tree_layout.py` and `tape_alignment.py` are shared geometry helpers imported by the
plotting scripts rather than run directly.

---

## Data availability

| What | Where |
| --- | --- |
| Raw sequencing | GEO **GSE341627** |
| Trees, metadata, tape matrices, founder genotypes | `support_data/` in this repository |
| Per-panel source data | `figures_data/` |
| Archived release | Zenodo — see the Data and code availability statement in the paper |

## License

GPL-3.0. See `LICENSE`.
