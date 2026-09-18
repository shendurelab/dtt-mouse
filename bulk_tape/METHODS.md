# Methods → code map

Every claim the Methods section makes about the bulk tape data, mapped onto the
function that implements it and the constant that parameterises it. All
constants live in [`config.py`](config.py).

## Input sequencing data

| Methods statement | Where |
|---|---|
| one amplicon FASTQ per embryo | `config.EMBRYOS` |
| embryos #2 and #3 single-end; #6 paired-end with the amplicon on R2 | `config.EMBRYOS["DTTz_6_S6"]` — the reverse complement of R2 is a drop-in single-end amplicon; R1 is a short 3′ mate with no barcode or terminal landmark and is unused |
| each read spans the integration barcode, the six-monomer array and the 3′ terminal landmark | `config` read-architecture diagram; `tape/parser.py` |

## Integration barcode and insertion vocabulary (`tape/reference.py`)

| Methods statement | Function | Constant | Value |
|---|---|---|---|
| 12-bp tag, `NNNNNAANNNNN`, between the `GGCATG` and `CTCTGG` flanks | `barcode.extract_bc` | `BC_LEFT`, `BC_RIGHT`, `BC_RIGHT_SPACER`, `BC_LEN` | `GGCATG`, `CTCTGG`, `G`, 12 |
| whitelist derived from ~2M sampled reads | `reference.derive` pass 1 | `SAMPLE_N` | 2,000,000 |
| barcodes above ≥0.2% frequency | `reference.derive` pass 1 | `BC_MIN_FREQ` | 0.002 |
| separated by ≥3 mismatches, excluding single-mismatch shadows | `reference.derive` pass 1 | `BC_MIN_SEP` | 3 |
| candidate insertions counted per (barcode × monomer position) **across all reads** | `reference.derive` pass 2 | — | pass 2 iterates the full fastq, unlike pass 1 |
| admitted at ≥10,000 reads at some single (barcode × position) | `reference.derive` | `VOCAB_PEAK_MIN` | 10,000 |

The barcode whitelist is frequency-based, so a sample suffices. The vocabulary
rule is an **absolute** per-(barcode × position) count, which scales with
sequencing depth, so it is measured over all reads — on a lightly edited library
the peaks of a 2M-read sample would fall below the fixed threshold.

Observed separation in the shipped references, lowest admitted peak against
highest rejected peak at any single (barcode × position):

| embryo | admitted | lowest admitted | highest rejected | gap |
|---|---|---|---|---|
| #2 | 11 | 14,589 | 2,350 | 6.2× |
| #3 | 13 | 72,667 | 5,446 | 13.3× |
| #6 | 1 | 28,813 | 2,997 | 9.6× |

`tests/test_pipeline.py` re-derives all of these from `refs/` and asserts that
the 10,000 threshold falls inside each gap, so these figures cannot go stale.

## Read parsing (`tape/parser.py`, driver `scripts/02_parse.py`)

| Methods statement | Function | Constant | Value |
|---|---|---|---|
| barcode corrected to the nearest whitelist entry within 2 mismatches, only when unambiguous | `barcode.correct` | `BC_CORRECT_MAX` | 2 |
| array resolved forward from the 5′ anchor, using ordered editing to infer that everything after an unedited monomer is unedited | `parser.walk_robust`, forward pass | `ANCHOR`, `GAP_MIN`, `GAP_MAX` | `GGTGAG`, 13, 24 |
| and backward from the 3′ terminal landmark, recovering distal monomers of heavily edited molecules | `parser.walk_robust`, backward pass; `parser.find_terminal` | `TERM`, `TERM_LONG` | `GGTGAGCCAC`, `GGTGAGCCACAAGGCAAT` (fuzzy, ≤2 mismatches) |
| reads with an edit following an unedited monomer are discarded | `parser.walk_robust` ordering check | — | — |
| valid = whitelisted barcode **and** all six monomers resolved **and** every edit in the vocabulary | `scripts/02_parse.py` | `N_SITES` | 6 |
| valid reads tallied per barcode over their 6-monomer genotype | `scripts/02_parse.py` → `patterns_full.pkl` | — | — |

**Open defect in this step:** the forward and backward passes number the sites
from opposite ends of the anchor list, so a read that loses its first anchor is
parsed one site out of frame and ends up with site 1 dropped and the last symbol
duplicated. See "Known issues" in the README, the block above the backward pass
in `tape/parser.py`, and `tests/test_known_issue_frame_shift.py`.

The backward pass runs only when the forward pass did not reach an unedited
monomer — once a monomer is unedited the ordered-editing rule already fixes every
monomer downstream of it.

## De-noising (`tape/denoise.py`, driver `scripts/03_denoise.py`)

Applied in this order, per integration barcode:

| Step | Function | Constant | Value |
|---|---|---|---|
| barcodes with <2,000 total valid reads are not analysed | `scripts/03_denoise.py` | `MIN_BC_READS` | 2,000 |
| (i) error-collapse into a single-mismatch neighbour ≥2-fold more abundant | `denoise.collapse`, `denoise.neighbors` | `COLLAPSE_RATIO` | 2.0 |
| (ii) chimera removal of fully-edited genotypes explained as prefix(A)+suffix(B) | `denoise.dechimera` | — | — |
| (iii) minimum of 2 supporting reads | `denoise.denoise_all` | `MIN_READS` | 2 |
| (iv) cross-barcode filter on genotypes with ≥3 edits whose read support is highest in another barcode | `denoise.denoise_all` | `CROSS_BC_DEPTH` | 3 |

Two details the Methods summarises: a **neighbour** in step (i) is a single-base
substitution *inside one monomer's insertion symbol* — an unedited monomer is
never mutated into an edited one, since it carries no sequence to misread. In
step (ii) the walk is most-abundant-first and only already-seen, more abundant
genotypes can serve as parents; partially-edited genotypes are exempt because
their unedited tails are legitimately shared by many real lineages.

## Integration copy number (`tape/copynumber.py`)

`estimate_copies` takes per-barcode total reads (over the de-noised genotypes)
and returns integer copies, the single-copy unit and their sum. The unit starts
as the median over all barcodes and is then re-taken, up to three times, as the
median of the barcodes that round to one copy — robust when several barcodes are
multi-copy. Copies are `max(1, round(reads / unit))`; total integrations is their
sum; cells are genomic equivalents ÷ copies.

Result: embryo #3 is 11 barcodes all at one copy (11 integrations); embryo #2 is
17 barcodes carrying 35 integrations (one at 8 copies, two at 6, one at 2).

## Cell-number estimation (`tape/depth.py`)

| Methods statement | Function | Constant | Value |
|---|---|---|---|
| λ = mode of the per-genotype read-count distribution in log space | `depth.estimate_lambda` — Gaussian KDE over log10 counts, peak of a 400-point grid | `LAMBDA_KDE_BW` | 0.3 |
| … over the genotypes carrying enough reads to define the mode | `depth.estimate_lambda` | `LAMBDA_MIN_READS` | 200 |
| cells carrying a genotype = round(reads / λ) | `depth.cells_per_pattern` | — | — |
| genomic equivalents per integration = their sum | `depth.genomic_equivalents` | — | — |

Fewer than five genotypes above the floor falls back to the peak of a smoothed
40-bin log10 histogram over all genotypes (`depth._hist_mode`). Embryo #3 gives
179–366 genomic equivalents per integration (mean 276).

## Per-integration lineage trees (`tape/editchain.py`, `tape/tree.py`)

| Methods statement | Function | Constant | Value |
|---|---|---|---|
| prefix tree over the ordered edit symbols; root-to-tip paths are the editing history; tip depth = monomers written | `editchain.edit_chain`, `editchain.build_trie` | `N_SITES` | 6 |
| 3′-truncated reads attributed to the deeper lineage, only when it is ≥4× more abundant | `editchain.fold_dropout` | `DROPOUT_FOLD_RATIO` | 4.0 |
| one lineage tree per integration barcode | `tree.draw_integration_tree`, `tree.render_all` | `TREE_TOPK`, `PALETTE` | 28, 12 colours |

`fold_dropout` compares the dominant child's whole subtree weight against the
interior node's own tip weight and, when it clears the ratio, reassigns that tip
weight to the deepest dominant descendant; it runs bottom-up so folds cascade.
Setting the ratio to 0 disables folding. For embryo #3 this takes 4,989 clean
lineages to 4,542.

## Rarefaction (`tape/rarefaction.py`)

Hurlbert interpolation up to the observed number of cells, Chao1 extrapolation
beyond it (`RAREFACTION_EXTRAP` = 5×). Not part of the Methods paragraph above;
included because it consumes the same per-integration cell counts.

## Values reproduced by `tests/test_pipeline.py`

| | embryo #2 | embryo #3 | embryo #6 |
|---|---|---|---|
| whitelisted integration barcodes | 17 | 11 | 5 |
| insertion vocabulary | 11 | 13 | 1 |
| barcodes ≥ `MIN_BC_READS` | 17 | 11 | 5 |
| clean lineages | 1,562 | 4,989 | 12 |
| after edit-chain folding | — | 4,542 | — |
| genomic copy number | 8/6/6/2 + 13×1 | all 1 | all 1 |
| total integrations | 35 | 11 | 5 |
| monomer-6 edit rate | 3.20% | 53.23% | 0.08% |
