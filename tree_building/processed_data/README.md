# Input files for tree building

All files under this directory are **outputs of the tape_pipeline code**

## Files and where they're used downstream

The following files are used as input for building the backbone tree (~600k cells), which are the input tapes for Blastomere A (B1) and Blastomere B (B2):

- e3v5v6.B1_tape_consensus.ge7_founderok.tsv.gz
- e3v5v6.B2_tape_consensus.ge7_founderok.tsv.gz

plus the sites that define the blastomere specific edits:
- e3v5v6.defining_sites.tsv


The following files are used to build the full tree (~1.3 M cells):
- e3v5v6.B1_tape_consensus.tsv.gz
- e3v5v6.B2_tape_consensus.tsv.gz

plus each integration's blastomere-founder genotype, to exclude founder-shared
edits from the placement QC metrics:
- e3v5v6.supp_table1_founder_genotypes.csv

Also used for tree dating -- a literature-derived embryo cell-count table,
not a tape_pipeline output:
- sample_matched_ceiling_sourced.csv
