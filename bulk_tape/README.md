# Bulk TAPE analysis pipeline

Reference implementation of the **primary analysis of bulk tape** described in the
Methods: from raw amplicon reads through per-integration lineage trees, for the
DNA Typewriter (TAPE) recorder in mouse embryos.

Everything the Methods section asserts about the bulk data is produced by the six
numbered stages below. Every threshold lives in [`config.py`](config.py) and
nothing else hard-codes one; [`METHODS.md`](METHODS.md) maps each sentence of the
Methods onto the function and the constant that implement it.

---

## What TAPE is (the minimum you need)

TAPE is a tandem array of six monomer "sites". A prime editor writes a short
insertion (a 3-bp *symbol* plus a constant `GGA` key) into the single active
site, which shifts the write-head one site downstream. Editing is therefore
**strictly ordered and unidirectional**: site *k* is written only after sites
1…*k*−1. A molecule's genotype is the ordered list of insertions up to the first
unedited site.

Each genomic **integration** of the recorder carries its own 12-bp *integration
barcode* (distinct from the monomer insertions); one integration is one clonal
lineage origin.

```
5' ...GGCATG [integration barcode, 12 bp] G CTCTGG ATGAT
   a1 site1 a2 site2 a3 site3 a4 site4 a5 site5 a6 site6  TERM  3'-const...

   monomer anchor (a1..a6) = GGTGAGCACG      terminal monomer = GGTGAGCCAC
   unedited junction       = ...TGAT         edited junction  = <NNN><GGA>TGAT
```

## Data flow

```
raw fastq.gz
  │  01_derive_reference.py    sample reads -> barcode whitelist + insertion vocabulary
  ▼
refs/<tag>.tapebc_whitelist.tsv, refs/<tag>.nnn_vocab.tsv         (checked in)
  │  02_parse.py               robust 5'-forward + 3'-backward parse; keep valid reads
  ▼
data/<tag>.patterns_full.pkl   {barcode: {6-site genotype: read count}}
  + per-site / per-barcode edit rates, read funnel        (checked in: example_data/)
  │  03_denoise.py             collapse -> chimera -> >=2 reads -> cross-barcode
  ▼
data/<tag>.clean_patterns.pkl  + de-noising funnel
  │  04_depth_model.py         lambda, genomic equivalents, integration copy number
  ▼
data/<tag>.depth_model.tsv
  │  05_integration_trees.py   edit-chain trie per barcode -> lineage trees
  │  06_rarefaction.py         Hurlbert interpolation + Chao1 extrapolation
  ▼
data/figures/<tag>.integration_trees.png, tables/<tag>.rarefaction_*.csv
```

Each stage reads only the previous stage's output, so you can re-run from any
point. **Stages 03–06 need no raw data** — they run from the committed
`example_data/*.patterns_full.pkl`.

## Install and verify

```bash
python3 -m venv .venv && ./.venv/bin/pip install -r requirements.txt
./.venv/bin/python tests/test_pipeline.py
```

`tests/test_pipeline.py` re-runs stages 03–06 from the committed inputs and
asserts the published values (embryo #3: 11 integrations, 4,989 clean lineages,
4,542 after edit-chain folding; embryo #2: 35 integrations across 17 barcodes).
It needs no raw data and no network.

## Running it

```bash
# stages 03-06 from the committed intermediates (no raw data required)
mkdir -p data && cp example_data/* data/
python3 scripts/03_denoise.py DTTz_3_S3
python3 scripts/04_depth_model.py DTTz_3_S3
python3 scripts/05_integration_trees.py DTTz_3_S3            # edit-chain folded (default)
python3 scripts/05_integration_trees.py DTTz_3_S3 --no-fold  # raw genotype view
python3 scripts/05_integration_trees.py DTTz_3_S3 --cells    # weight tips by cells, not reads
python3 scripts/06_rarefaction.py DTTz_3_S3

# full run from raw fastq (see "Raw data" below)
export TAPE_RAW_DIR=/path/to/fastqs
./run_bulk.sh DTTz_3_S3               # stages 01 -> 05
./run_bulk.sh DTTz_3_S3 --skip-ref    # reuse the checked-in refs/, skip stage 01

# figure / export tables
python3 scripts/bulk_figure_tables.py DTTz_3_S3 DTTz_2_S2
python3 scripts/site_insertion_spectrum.py DTTz_3_S3
python3 scripts/barcode_rank_abundance.py DTTz_1_S1 DTTz_2_S2 DTTz_3_S3 DTTz_4_S4
python3 scripts/export_bulk_lineage_genotypes.py DTTz_3_S3
```

Set `TAPE_RAW_DIR` (raw fastqs) and `TAPE_DATA_DIR` (generated intermediates) to
keep large files outside the repository.

## Layout

```
config.py            All constants, thresholds, paths, per-embryo settings.
                     The single source of truth -- no magic numbers elsewhere.
tape/                Core library.
  io.py              fastq iteration; whitelist/vocabulary loaders; pickle helpers
  barcode.py         integration-barcode extraction, Hamming distance, correction
  parser.py          junction grammar + robust 5'-forward / 3'-backward array parser
  reference.py       derive the whitelist + insertion vocabulary
  denoise.py         collapse + chimera + read floor + cross-barcode filter
  depth.py           lambda (log-mode) -> genomic equivalents / cells
  copynumber.py      per-barcode genomic copy number -> total number of integrations
  editchain.py       the edit-chain trie: the lineage model
  rarefaction.py     Hurlbert interpolation + Chao1 extrapolation
  tree.py            trie layout + per-integration lineage-tree rendering
scripts/             Stage drivers (01-06) and the figure/export scripts.
refs/                Checked-in whitelists + insertion vocabularies, one pair per embryo.
example_data/        Committed stage-02 outputs (per-barcode genotype counts plus the
                     read funnel and per-site edit rates), so stages 03-06 and the
                     figure scripts run out of the box.
tests/               Reproducibility check (no raw data needed).
run_bulk.sh          Driver for stages 01 -> 05.
```

## The embryos

| tag | embryo | sequencing | integrations | site-6 edit rate |
|---|---|---|---|---|
| `DTTz_2_S2` | #2 | single-end | 35 across 17 barcodes | 3.2% ("cold" recorder) |
| `DTTz_3_S3` | #3 | single-end | 11, all single-copy | 53.2% ("hot"; the lineage-analysis focus) |
| `DTTz_6_S6` | #6 | paired-end (see below) | 5 | 0.08% |
| `DTTz_1_S1`, `DTTz_4_S4` | #1, #4 | single-end | — | 0.12% / 0.00% (no recording) |

Embryo #6 was sequenced paired-end. Its amplicon is carried on R2, which is
reverse-complemented into a single-end amplicon fastq and then processed
identically to embryos #2 and #3; the R1 mate is a short 3′ read carrying no
barcode or terminal landmark and is unused.

## Raw data

The raw amplicon fastqs are not distributed with the code (they are large). Point
`TAPE_RAW_DIR` at them and add an `Embryo()` entry in `config.py` for any new
library. The committed `refs/` and `example_data/` are exactly what stages 01 and
02 produce from them, so the rest of the pipeline reproduces without them.

## Shared with the single-cell pipeline

The read body downstream of the integration barcode is identical in the bulk and
single-cell assays, so `parser.py`, `denoise.py` and `editchain.py` are shared
verbatim between them, and the `refs/` derived here (bulk has deep, clean
per-integration coverage) are what the single-cell pipeline consumes for the same
embryo. Only the barcode read-off differs. The single-cell front-end is released
separately.

The edit-chain model is the unifying abstraction. A lineage is an ordered chain
of edit symbols; because editing is monotonic, a shallower chain is a *prefix* of
a deeper one. A read or molecule ending at an interior node is either a genuine
early stop or 3′ truncation of a deeper lineage. Bulk keeps many lineages per
integration (the trie *is* the tree) and folds only prefixes that look like
truncation, conservatively; single cell folds all prefixes onto the deepest
dominant tip to call one genotype per cell.

## Known issues

- **Backward-pass frame shift (open).** `walk_robust` numbers the monomer sites
  from anchors[0] going forward and from the terminal landmark going backward.
  A read that loses its *first* anchor is parsed one site out of frame by the
  forward pass and correctly registered by the backward pass; the two then read
  the same junction into adjacent slots, dropping site 1 and duplicating the
  last symbol. The read passes every downstream check, and de-noising cannot
  catch it (a frame shift is not a one-symbol neighbour, not a prefix/suffix
  join, and keeps its own integration barcode). On embryo #3 the signature
  accounts for roughly 0.5% of fully-edited reads above a chance-match
  background, in genotypes at a few percent of their parent's read depth.
  Documented in `tape/parser.py`; reproduced and pinned by
  `tests/test_known_issue_frame_shift.py`.

## Caveats

- **Bulk has no UMIs.** Read-level de-noising (collapse + chimera) is the
  substitute; it cannot fully separate genuine early-stop lineages from PCR
  chimeras or 3′ truncation. The single-cell assay resolves this with UMIs.
- **Array length is fixed at six in the parser.** A minority of embryo-#3 arrays
  are physically expanded beyond six monomers; the parser reads the first six.
  Variable array length is not modelled.
- **λ and cell counts are estimates.** They are stable across the eleven
  single-copy integrations of embryo #3 but depend on the shape of the read-depth
  distribution, and are least reliable for the low-depth multi-copy barcodes of
  embryo #2.

## License

GPL-3.0. See [LICENSE](LICENSE).
