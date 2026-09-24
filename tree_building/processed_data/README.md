# Input files for tree building

Files here are DERIVED within this repo from support_data/ (not raw
tape_pipeline deliverables -- those live in the repo root's support_data/).

## Files and where they're used downstream

The following files are used as input for building the backbone tree (~656k
cells), which are the input tapes for Blastomere A (B1) and Blastomere B (B2).
Derived from support_data/ by 1_build_nj_backbone/00_filter_highqual_consensus.R
(pass_qc & n_loci>=7 & founder-consistent):

- e3v8.B1_tape_consensus.ge7_founderok.tsv.gz
- e3v8.B2_tape_consensus.ge7_founderok.tsv.gz

The sites that define the blastomere-specific edits (used to build the
synthetic founder-genotype root outgroup) ship directly in support_data/:
- ../support_data/e3v8.defining_sites.tsv

Each integration's blastomere-founder genotype, used both by
00_filter_highqual_consensus.R (founder-consistency filter) and downstream by
4_full_tree_placement (to exclude founder-shared edits from placement QC
metrics). Reshaped from support_data/e3v8.supp_table1_founder_genotypes.csv
(long format) into wide format by
1_build_nj_backbone/00a_reshape_v8_founder_genotypes.R:
- e3v8.supp_table1_founder_genotypes.wide.csv

Also used for tree dating -- a literature-derived embryo cell-count table,
not a tape_pipeline output:
- cell_counts_full_kojima_cao.csv
- sample_matched_ceiling_sourced.csv
