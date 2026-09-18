# Single-cell circtape analysis pipeline

Reference implementation of the **primary analysis of single-cell circtape**
described in the Methods: from raw sci-RNA-seq3 tape amplicon reads through a
cell × integration character matrix, to the two per-blastomere matrices that are
the input to tree building. Embryo #3, the only embryo with single-cell data.

Every threshold lives in [`config.py`](config.py) and nothing else hard-codes
one; [`METHODS.md`](METHODS.md) maps each sentence of the Methods onto the
function and the constant that implement it.

This is the single-cell companion to the **bulk tape** package. The read body
downstream of the integration barcode is identical in the two assays, so
`tape/parser.py`, `tape/denoise.py` and `tape/editchain.py` are the same code,
and the embryo-#3 barcode whitelist and insertion vocabulary in `refs/` are the
bulk-derived reference sets, shipped verbatim. What is single-cell-specific:
barcode extraction (the read starts at the barcode), the UMI layer, the per-cell
consensus, and blastomere routing.

---

## What the recorder is (the minimum you need)

A tape is a tandem array of six monomer "sites". A prime editor writes a short
insertion (a 3-bp *symbol* plus a constant `GGA` key) into the single active
site, which shifts the write-head one site downstream, so editing is **strictly
ordered and unidirectional**: site *k* is written only after sites 1…*k*−1. A
molecule's genotype is the ordered list of insertions up to the first unedited
site. Each genomic **integration** carries its own 12-bp *integration barcode*;
embryo #3 has 11.

Because editing is ordered, a shorter chain is always a **prefix** of a longer
one. Within one cell that means a molecule with fewer edits is 3′ dropout or
transcript carry-over of the cell's own deeper genotype, **not** a competing
allele — the fact the consensus model is built on.

```
one single-cell amplicon read (R2):

5' [integration barcode, 12 bp] G CTCTGG ATGAT
   a1 site1 a2 site2 a3 site3 a4 site4 a5 site5 a6 site6  TERM  3'-const...

   monomer anchor (a1..a6) = GGTGAGCACG      terminal monomer = GGTGAGCCAC
   unedited junction       = ...TGAT         edited junction  = <NNN><GGA>TGAT
```

A single-cell read **starts at the barcode** — there is no 5′ `GGCATG` flank as
in the bulk amplicon — so the barcode is the 12 bp immediately 5′ of the constant
`G` that precedes `CTCTGG`.

## Data flow

```
raw fastq.gz  (one single-end run; one paired run)
  │  00_merge_pairs.py        fastp overlap-merge, R2 fallback  [paired run only]
  ▼
data/seq8_long.sc_reads.fastq.gz     amplicon-oriented reads
  │  sc_consensus_par.py      BOTH runs in one pass:
  │                             parse + whitelist + vocabulary  (shared with bulk)
  │                             tally per (cell, integration, UMI)
  │                             UMI merge (Hamming 1) -> >=2 reads per molecule
  │                             de-noise molecule counts (collapse, chimera, >=2)
  │                             per-cell consensus (>=3 molecules, dominance >=0.9)
  ▼
data/e3v8.cell_consensus.tsv   cell x 11 integrations + n_loci / n_doublet_loci /
data/e3v8.sc_qc.tsv            mean_dominance                    (1,559,419 cells)
  │  blastomere_route.py      founder-allele routing to blastomere A/B, bleed mask,
  │                           set-aside, and the pass_qc flag
  ▼
data/e3v8.B1_tape_consensus.tsv   889,849 cells   ->  distance / NJ / dating
data/e3v8.B2_tape_consensus.tsv   653,335 cells
  + set_aside.tsv, routing_labels.tsv, defining_sites.tsv, routing_report.txt
```

`sc_consensus_par.py` is the parallel drop-in for `sc_consensus.py`: identical
output (verified with `--limit`), only the per-read parse is fanned across
processes. The published matrices were made with `sc_consensus_par.py --procs 8`.

## Run it

```bash
export TAPE_RAW_DIR=/path/to/fastqs      # large; not in this repo
export TAPE_DATA_DIR=/path/to/workdir    # intermediates are GB-scale
./run_sc.sh
```

Inputs expected under `TAPE_RAW_DIR`:

| file | what |
|---|---|
| `Tape_merged.seq8.R2.fastq.gz` | single-end run, 183-bp amplicon read; **both** sci-RNA-seq3 experiments (25 PCR plates) |
| `Tape_merged.seq8_long.R{1,2}.fastq.gz` | paired run, 2×250; **first** experiment only, re-sequenced at greater length |
| `cell_metadata.v8.txt` | cell → celltype / trajectory / UMAP; the recovery denominator (1,584,848 nuclei) |

`fastp` must be on `PATH` for stage 00 (`FASTP_BIN` to override). `pigz` is used
when present. Python ≥3.9; `numpy` for routing, nothing else (see
`requirements.txt`).

## Verify a checkout — no raw data, no network, ~20 s

```bash
python3 tests/test_pipeline.py                 # parse, thresholds, routing, published values
python3 tests/test_known_issue_frame_shift.py  # the pinned parser defect + its suppression
python3 check_sync_with_tape_pipeline.py       # no drift vs the working tree (in-repo only)
```

`tests/test_pipeline.py` re-runs **blastomere routing** on
`example_data/e3v8_subsample.cell_consensus.tsv.gz`, a 1/40 random cell
subsample (38,732 cells) of the published callset. It derives the **identical 20
defining founder alleles** as the full run and reproduces every published rate —
99.0% routed, 1.04% set aside, 84.0% passing QC, 85.3% / 82.2% per blastomere.

## Published values

| quantity | value |
|---|---|
| nuclei passing transcriptome QC | 1,584,848 |
| cells with ≥1 called tape genotype | 1,559,419 (98.4%) |
| median quality-passing UMIs (molecules) per cell | 51 |
| matrix fill | 52.9% (9,217,411 of 11 × 1,584,848) |
| per-integration recovery | 48.5–69.3% for nine, 21.8% and 3.7% for two |
| defining founder alleles | 20, over all 11 integrations |
| routed to A / B | 889,849 (57.1%) / 653,335 (41.9%) |
| bleed-masked cells (one integration → NA) | 58,711 |
| set aside (no founder / tie / ≥2 bleed) | 16,235 (1.04%) = 5,308 / 8,451 / 2,476 |
| passing QC (≥4 loci, ≤1 doublet, dominance ≥0.95) | 1,296,799 (84.0%) |
| → backbone-eligible (≥7 of 11) / placement-eligible (4–6) | 663,785 / 633,014 |
| lost to the founder-state filter | 8,084 (1.2%) backbone-eligible, 7,574 (1.2%) placement-eligible |
| backbone tips / cells placed / full tree | 655,701 / 625,440 / 1,281,141 |

## Known issues and caveats

- **Parser frame shift (open, not fixed).** A read that loses its *first* monomer
  anchor is parsed one site out of frame; monomer 1 is lost and the last symbol
  duplicated. The parser is shared verbatim with bulk, so the defect is
  identical, but the single-cell funnel suppresses it: a shifted genotype must
  win its molecule's read consensus, then recur in ≥2 independent molecules of
  the same (cell, integration) to survive the de-noising floor, then win a
  ≥3-molecule / ≥0.9-dominance contest. Measured upper bound in the published
  callset: **≤0.16% of calls** (`scripts/frame_shift_bound.py`;
  `tests/test_known_issue_frame_shift.py` pins the behaviour and the fix
  direction).
- **`n_loci` is written before bleed masking.** It therefore overstates the
  number of non-`NA` genotypes by one in each of the 58,711 bleed-masked cells,
  and the `pass_qc` flag inherits that. For true recovery, count non-`NA` tape
  columns.
- **Set-aside cells are not fate-neutral** — primitive erythroid 7.6×,
  cardiomyocytes 4.2×, endothelium 3.2× over-represented, CNS neurons depleted
  0.70×. They are handed off in `set_aside.tsv`, not silently dropped.
- **The founder-state filter shipped here is a reconstruction.** The published
  trees were built from `*.ge7_founderok.tsv.gz` inputs produced by code that is
  not in this package; `scripts/founder_state_filter.py` re-derives the founder
  chains and reproduces the published cell losses with zero false positives, but
  misses 336 of 15,658 (2%). Use it to audit, not to substitute.
- **The cross-barcode swap filter is not applied per cell.** It is a property of
  the bulk population; a single cell has no such comparison to make.
- **Embryo #2 was co-processed** in the first sci-RNA-seq3 experiment. Its 17
  integration barcodes are ≥4 mismatches from every embryo-#3 barcode, twice the
  correction radius, so embryo-#2 molecules are dropped at barcode correction;
  they are ~2.5% of barcode-extractable reads.

## Provenance — important

- The published matrices were generated with **`blastomere_route.py` at commit
  `6ee057f`**, which is the version shipped here. A later revision of that script
  replaces the `multibleed` rule with a stricter side-private
  `founder_conflict` veto; it is a defensible improvement but is **not** what the
  paper used (on the same matrix it yields 1,525,793 routed and 2.16% set aside).
  `check_sync_with_tape_pipeline.py` compares against `6ee057f` deliberately.
- Thresholds are the published ones: `MIN_READS_PER_UMI=2`,
  `SC_MIN_MOLECULES=3`, `SC_DOMINANCE_MIN=0.9`.
- `refs/DTTz_3_S3.*` are the bulk-derived reference sets, unchanged.
