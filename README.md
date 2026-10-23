Dataset:
Reviewed bacterial proteins with EC annotations downloaded from UniProtKB/Swiss-Prot.

Cleaning:
Run code/clean_ec3.py on the original UniProt dataset.

Final cleaned files:
data/proteins_ec3.tsv
data/proteins_ec3.fasta

Splitting:
Run code/make_splits_ec3.py.

MMseqs2 is required for 50% sequence-identity clustering.
Random seed: 355.
Minimum training examples per EC3 class: 10.