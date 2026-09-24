

### Correction of RT primer cross-contamination 
### We identified a contamination issue affecting a subset of the RT primer wells used during the sci-RNA-seq3 experiments. Briefly, an unusually high number of cell pairs shared both ligation and PCR barcodes but not RT barcodes. These pairs, which collectively impacted ~10% of cells, overwhelmingly: (i) derived from a few dozen RT well pairs, (ii) shared cell type, and (iii) had lower UMI counts. We recognized this as the signature of “cell splitting” (i.e., a single cell or nucleus manifesting in the data as two cells because its molecules are split across two distinct RT barcodes), effectively the opposite of “cell doublets”. To resolve this, we merged reads associated with the splitting artifact into single cells, across 54 RT wells forming 23 confusable groups (2-6 wells per group; cross-contamination occurred pairwise within each group), and reprocessed the data from raw reads.

### RT_barcode_pairs.txt is available at "support_data" folder
### https://github.com/shendurelab/dtt-mouse/tree/main/support_data

import os, gzip, re, sys

RT_barcode_correction = {}
with open(f"{work_path}/RT_barcode_pairs.txt") as f:
    for line in f:
        a, b = line.rstrip().split('\t')
        RT_barcode_correction[a] = b

batch_id = int(sys.argv[1])

experiment_id = sys.argv[2]

if experiment_id == "experiment2_20260713_seq2_XY":
    target_barcode = "embryo3"
else:
    target_barcode = "DTT_Z_858923_E2"

in_dir  = f"{work_path}/{experiment_id}/nobackup/output_{batch_id}/UMI_attach"
out_dir = f"{work_path}/{experiment_id}/nobackup/output_{batch_id}/UMI_attach_correct"

for fastq in os.listdir(in_dir):
    if not fastq.endswith(".fastq.gz"):
        continue

    with gzip.open(f"{in_dir}/{fastq}", "rt") as fin, \
         gzip.open(f"{out_dir}/{fastq}", "wt") as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            seq  = fin.readline()
            plus = fin.readline()
            qual = fin.readline()

            fields = header.rstrip().split(',')
            if fields[-1] != target_barcode:
                continue

            left, well = fields[0].split('_RT-')
            if well in RT_barcode_correction:
                well = RT_barcode_correction[well]
            fields[0] = f"{left}_RT-{well}"

            fout.write(','.join(fields) + '\n')
            fout.write(seq)
            fout.write(plus)
            fout.write(qual)

