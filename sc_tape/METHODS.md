# Methods → code map (single-cell circtape)

Every claim the Methods section "Primary analysis of single-cell circtape" makes,
mapped onto the function that implements it and the constant that parameterises
it. All constants live in [`config.py`](config.py). Constants marked *(shared)*
have the same value in the bulk package.

## Input sequencing data

| Methods statement | Where |
|---|---|
| amplicons from the same indexed cDNA as the transcriptome libraries, both sci-RNA-seq3 experiments | `config.SC_SAMPLES` — `seq8` covers both experiments' cells, `seq8_long` the first only |
| all single-cell data derive from embryo #3 | `config.EMBRYOS` holds only `DTTz_3_S3` |
| embryo #2 molecules are discarded at barcode correction | `barcode.correct` with `BC_CORRECT_MAX` = 2; embryo #2's 17 barcodes are ≥4 mismatches from every embryo-#3 barcode |
| single-end 35/10/10/183 and paired 250/10/10/250 | `SC_SAMPLES[...].layout`; the amplicon is always the R2 read |
| header carries the cell identifier and an 8-bp UMI | `sc.parse_header`; `UMI_LEN` = 8 |

## Read merging (`tape/merge.py`)

| Methods statement | Function | Constant | Value |
|---|---|---|---|
| R2 is the amplicon read, R1 its 3′/poly-A mate | `merge_sample` docstring | — | — |
| fastp `--merge --correction`, paired-end adapter detection, quality- and length-filtering disabled | `run_fastp` | `FASTP_BIN` | `fastp` |
| minimum overlap 15 bp | `run_fastp` | `FASTP_OVERLAP_LEN_REQUIRE` | 15 |
| ≤8 mismatches in the overlap | `run_fastp` | `FASTP_OVERLAP_DIFF_LIMIT` | 8 |
| ≤30% mismatched bases in the overlap | `run_fastp` | `FASTP_OVERLAP_DIFF_PCT_LIMIT` | 30 |
| merge-with-fallback: merged read where a pair merged, else R2 | `combine_oriented` | — | returns the merge rate (22.6% on the published run) |
| merged reads reverse-complemented to amplicon orientation | `combine_oriented` → `_revcomp` | — | fastp emits merged reads in R1 orientation |
| both runs pooled by (cell, integration, UMI) in one pass | `scripts/sc_consensus_par.py::main` | — | one `mols[(cell, bc)][umi][pattern]` dict across samples |

## Reference sets, parsing and de-noising

Driver `scripts/sc_consensus_par.py` (or `sc_consensus.py`), using `tape/sc.py`
with the shared `tape/parser.py` and `tape/denoise.py`.

| Methods statement | Function | Constant | Value |
|---|---|---|---|
| whitelist and vocabulary taken from the bulk analysis | `io.load_whitelist`, `io.load_vocab` on `refs/DTTz_3_S3.*` | — | 11 barcodes, 13 symbols |
| reads begin at the barcode; barcode is the 12 bp 5′ of the constant G before CTCTGG | `sc.extract_bc_sc` | `BC_RIGHT`, `BC_RIGHT_SPACER`, `BC_LEN` *(shared)* | `CTCTGG`, `G`, 12 |
| corrected to the whitelist as in the bulk pipeline | `barcode.correct` | `BC_CORRECT_MAX` *(shared)* | 2, unique nearest only |
| array resolved and validated identically to bulk | `parser.walk_robust`, then the vocabulary test in `_parse_batch` | `N_SITES`, `GAP_MIN`, `GAP_MAX`, `ANCHOR`, `TERM`, `TERM_LONG` *(shared)* | 6, 13, 24, … |
| reads tallied per cell, per integration, per UMI | `_parse_batch` → `mols[(cell, bc)][umi][pattern]` | — | — |
| a UMI within one mismatch of a more abundant UMI merges into it | `sc.collapse_umis` | `UMI_MERGE_DIST` | 1 |
| a merged UMI with <2 reads is discarded | `sc.collapse_umis` | `MIN_READS_PER_UMI` | 2 |
| each molecule takes the dominant genotype among its reads | `sc.collapse_umis` | — | ties break by insertion order |
| error-collapse of single-mismatch neighbours | `denoise.collapse` | `COLLAPSE_RATIO` *(shared)* | 2.0, rarest-first |
| prefix/suffix chimera removal | `denoise.dechimera` | — | fully-edited genotypes only |
| genotypes with <2 molecules discarded | `denoise.denoise_barcode` | `MIN_READS` *(shared)* | 2 (molecules here, reads in bulk) |
| the cross-barcode swap filter is **not** applied per cell | `denoise_barcode` is called, not `denoise_all` | `CROSS_BC_DEPTH` | 3, unused here |

## Per-cell consensus (`tape/sc.py` + `tape/editchain.py`)

| Methods statement | Function | Constant | Value |
|---|---|---|---|
| one genotype per cell and integration, edit-chain model | `sc.consensus_locus` → `editchain.dominant_lineage` | — | — |
| a shallower chain on the same path is a prefix (carry-over), not a distinct allele | `editchain.build_trie`; prefixes lie on the root→tip path | — | editing is monotonic |
| the consensus descends the most-supported chain and reports the deepest dominant genotype | `dominant_lineage` — descends by max subtree to a leaf | — | — |
| two symbols competing at the same next monomer define a branch | `dominant_lineage` — any non-top child of a node on the path | — | — |
| a locus is flagged a putative doublet when an off-path branch carries ≥2 molecules | `dominant_lineage` → `branched` | — | 2, **no fraction test** |
| branching does not truncate the call | `dominant_lineage` never stops at a branch | — | contrast `editchain.consensus`, which does and is unused |
| dominance = path molecules / all molecules = 1 − branching/all | `dominant_lineage` | — | — |
| called lineage needs ≥3 molecules | `consensus_locus` (#2) | `SC_MIN_MOLECULES` | 3 |
| dominance must be ≥0.9, else the locus is left missing | `consensus_locus` (#3) | `SC_DOMINANCE_MIN` | 0.9 |
| per-cell metrics n_loci / n_doublet_loci / mean_dominance | `sc_consensus_par.main` output writer | — | over **called** loci only |

`BRANCH_MIN_FRAC` and `BRANCH_MIN_READS` parameterise `editchain.is_branch` and
`editchain.consensus`, which the single-cell caller does **not** use. Reading the
0.25 fraction as part of the doublet rule is the one easy mistake here.

## Blastomere routing (`scripts/blastomere_route.py`, commit `6ee057f`)

| Methods statement | Function | Constant | Value |
|---|---|---|---|
| defining alleles taken at monomer 1 | `encode_founder` | `FOUNDER_SITE` | 0 (site 1) |
| well-covered cells only, so shared missingness cannot drive the split | `derive_defining` | `CONF_MIN_LOCI` | 8 recovered integrations |
| leading principal component of rarity-weighted monomer-1 alleles | `confident_axis` | — | one-hot × √(−log allele frequency), L2-normalised per cell |
| split at the threshold maximising between-group variance | `confident_axis` | — | Otsu / balance-weighted between-class variance |
| allele present in ≥85% of recovered calls on one side | `derive_defining` | `FIX_THRESH` | 0.85 (observed ≥0.969) |
| and enriched by ≥0.30 over the other | `derive_defining` | `DISC_MIN` | 0.30 (observed ≤0.071 on the other side) |
| — | `derive_defining` | `MIN_REC` | 20 recovered cells per side for a locus to be considered |
| routed by majority vote, ≥1 own-side and ≤1 opposite-side allele | `route` | — | `w ≥ 1 & w > l & l ≤ 1` |
| the single opposite-side call is ambient bleed, that integration set to missing | `route` → `mask_g`, written by `write_matrix` | `MASK_TOKEN` | `NA`; 58,711 cells |
| set aside on no founder, a tie, or ≥2 opposite-side alleles | `route` | — | 5,308 / 8,451 / 2,476 |
| A = the larger routed group | `orient` | — | B1 = larger; A = B1, B = B2 |
| pass_qc = ≥4 loci, ≤1 doublet locus, mean_dominance ≥0.95 | `cell_pass_qc` | CLI `--min-n-loci`, `--max-doublet-loci`, `--min-mean-dom` | 4, 1, 0.95 |
| all routed cells remain in the matrices; the distance step selects the passing set | `write_matrix` writes every kept cell with a `pass_qc` column | — | non-destructive |

### The last filter before tree building

| Methods statement | Where | Value |
|---|---|---|
| the QC-pass set divides by tape count into the two tree inputs | `n_loci >= 7` / `4 <= n_loci <= 6` on the released matrices | 663,785 + 633,014 = 1,296,799 |
| at both stages, no tape genotype contradicting its blastomere's founder state | `scripts/founder_state_filter.py` (**reconstruction** — see below) | removes 8,084 (1.2%) and 7,574 (1.2%) |
| …giving the 655,701-cell backbone and the 625,440 placed cells | the published trees | 655,701 + 625,440 = 1,281,141 |

"founder state" is the near-fixed prefix of edits a blastomere's founder wrote
before the first cleavage: the monomer-1 defining allele extended one monomer at a
time while the next monomer is edited and near-fixed among the cells matching the
chain so far. `tables/e3v8.founder_chains.tsv` holds the 20 derived chains, up to
six monomers deep (blastomere A, `ATATCAAATTGA`: `GATG-ACT-ACT-AAG-ACT-AAG`) —
the 2-cell editing burst, read directly off the matrix.

**The filter that produced the published trees is not in this package.** The tree
pipeline consumes files named `*_tape_consensus.ge7_founderok.tsv.gz`, produced
upstream by code held elsewhere. `founder_state_filter.py` is a reconstruction,
validated against the published tip sets: it flags 4,195/4,195 (A) and 3,729/3,889
(B) of the backbone-eligible cells absent from the backbone and 4,776/4,833 (A)
and 2,622/2,741 (B) of the placement-eligible cells absent from the placed set,
with **zero false positives** — it never flags a cell that is in a published tree.
The 336 it misses (2% of 15,658) are presumably a slightly looser near-fixed
threshold than the 0.99 used here.

## What this package does not cover

Everything downstream of the two per-blastomere matrices — the distance metric,
neighbour-joining, dating, and placement — is a separate R pipeline and is
described in "Reconstruction of the backbone tree".

Nor does it derive the reference sets: the whitelist and insertion vocabulary are
the bulk ones, and their derivation (`tape/reference.py`, `SAMPLE_N`,
`BC_MIN_FREQ`, `BC_MIN_SEP`, `VOCAB_PEAK_MIN`) lives in the bulk package.
