Dataset:
Reviewed bacterial proteins with EC annotations downloaded from UniProtKB/Swiss-Prot.

Cleaning:
Run code/clean_ec3.py on the original UniProt dataset.

Final cleaned files:
data/proteins_ec3.tsv
data/proteins_ec3.fasta

Splitting:
Run code/make_splits_ec3_final.py to create splits; do not rerun it to reproduce
the nearest-training analysis of the existing frozen splits.

MMseqs2 is required for 50% sequence-identity clustering.
Random seed: 355.
Minimum training examples per EC3 class: 10.

Nearest-training similarity:
Run `python code/nearest_train_identity.py` against the existing final splits.
See [the analysis report](analysis/nearest_train_identity/report.md) and
[methodological details](analysis/nearest_train_identity/README.md).
Regression checks: `python code/test_nearest_train_identity.py`.

Final analysis results, search logs, and validation records are versioned.
FASTA files and large raw hit tables use Git LFS (`git lfs install` before
cloning or pulling). MMseqs2 databases, temporary files, and preflight runs
remain local and can be regenerated.

Sequence-similarity strata:
Run `python code/make_similarity_bins.py` to create a new
`analysis/similarity_bins/` directory (existing outputs are never overwritten).
See [the binning decision and validation report](analysis/similarity_bins/report.md).
The common strata are <30%, 30-40%, 40-50%, 50-60%, 60-70%, and 70-100%,
with `no_hit` kept separate. Each numerical stratum has at least 100 proteins
per test group in every experiment. Regression checks:
`python code/test_make_similarity_bins.py`.
No classifier has been trained at this stage.
