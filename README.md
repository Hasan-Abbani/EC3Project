# EC3 enzyme-function prediction

Reviewed bacterial proteins with EC annotations were downloaded from
UniProtKB/Swiss-Prot (recorded release 2026_03). The dataset contains public
protein sequences and annotations. See [UniProt attribution and data licensing](DATA_LICENSE.md)
for the source, CC BY 4.0 terms, citation, and description of our modifications.


Project code is released under the MIT License. UniProt-derived data remain subject to the licensing terms described in [DATA_LICENSE.md](DATA_LICENSE.md).


## Setup and available artifacts

From the project root, install the analysis dependencies in a Python environment:

```text
python -m pip install numpy pandas matplotlib scikit-learn
```

Recorded analysis versions were NumPy 2.4.3, pandas 3.0.2, Matplotlib 3.10.8,
and scikit-learn 1.9.1. The embedding run used Python 3.12.10, PyTorch
2.11.0+cu126, and Transformers 4.57.6; the other analysis records include Python
3.13.12. These are the recorded run environments, rather than a single locked
environment for every stage.

FASTA files and large TSV files use Git LFS. Install Git LFS and run
`git lfs install` and `git lfs pull` after cloning. Frozen predictions, result
tables, figures, split assignments, and validation records are included.
Embeddings, fitted model arrays, dependency caches, and temporary databases are
kept locally and are excluded from Git. They must be regenerated or supplied
separately to repeat model training. The original UniProt download and cleaning
logs are not included.
The local `manuscript/` folder is also not included in the repository.

The current MMseqs2 scripts use Windows/WSL. Set `MMSEQS_WSL` in
`code/make_splits_ec3_final.py` and `MMSEQS` in `code/nearest_train_identity.py`
to your installed executable. `code/run_logreg.py` reads embeddings from
`analysis/esm2_embeddings/` (the default output of
`code/generate_esm2_embeddings.py`) and labels from `data/proteins_ec3.tsv`;
if yours are elsewhere, change `EMB_DIR` and `LABELS` near the top of the
script. The MMseqs2 path settings must be adjusted for another machine.

## Dataset construction

Run code/clean_ec3.py on the original UniProt dataset.

Final cleaned files:
data/proteins_ec3.tsv
data/proteins_ec3.fasta

## Frozen splits and similarity analysis

Run code/make_splits_ec3_final.py to create splits; do not rerun it to reproduce
the nearest-training analysis of the existing frozen splits.

MMseqs2 is required for 50% sequence-identity clustering.
Seen-test cluster-selection seed: 355.
Validation cluster-selection seed: 356 (`SEED + 1`).
Minimum training examples per EC3 class: 10.

Validate the existing final splits without changing files:

```text
python code/nearest_train_identity.py --validate-only
```

The full nearest-training search is performed by `code/nearest_train_identity.py`.
Do not rerun it over the archived results just to inspect them.
See [the analysis report](analysis/nearest_train_identity/report.md) and
[methodological details](analysis/nearest_train_identity/README.md).
Regression checks: `python code/test_nearest_train_identity.py`.

Final analysis results, search logs, and validation records are versioned.
Their parameters, source hashes, and software versions describe the recorded
runs. Historical usernames and local paths are intentionally retained.

## Sequence-similarity strata

Run `python code/make_similarity_bins.py` to create a new
`analysis/similarity_bins/` directory (existing outputs are never overwritten).
See [the binning decision and validation report](analysis/similarity_bins/report.md).
The common strata are <30%, 30-40%, 40-50%, 50-60%, 60-70%, and 70-100%,
with `no_hit` kept separate. Each numerical stratum has at least 100 proteins
per test group in every experiment. Regression checks:
`python code/test_make_similarity_bins.py`.
The strata were selected outcome-independently, before classifier performance
was examined.

## Embeddings and classifier training

Embedding extraction additionally requires PyTorch, Transformers, psutil,
filelock, and tqdm. Install the appropriate CPU or CUDA PyTorch build using the
[official PyTorch instructions](https://pytorch.org/get-started/locally/), then:

```text
python -m pip install transformers==4.57.6 psutil filelock tqdm
python code/generate_esm2_embeddings.py --pilot --output analysis/esm2_embeddings_pilot
```

The recorded full GPU command used:

```text
python code/generate_esm2_embeddings.py --output analysis/esm2_embeddings_cuda --device cuda --precision float16 --batch-size 1
```

Use `--resume` for an interrupted run. For a CPU run, use `--device cpu
--precision float32`. Inspect `progress.json` and `validation.json` before using
the generated embeddings. Numerical results can vary across hardware and
software environments.

After setting the classifier input paths described above, reproduce one
held-out-phylum model in a separate local output directory:

```text
python code/run_logreg.py --split splits_ec3/Actinomycetota/split_Actinomycetota.tsv --bins analysis/similarity_bins/Actinomycetota_similarity_bins.tsv --phylum Actinomycetota --C 0.001,0.01,0.1,1 --out .cache/reproduction/Actinomycetota
```

Repeat with Bacillota and Cyanobacteriota for the other experiments. Keep
retraining outputs separate from the frozen predictions.

## Bin-wise model evaluation

Run `python code/evaluate_similarity_bins.py` against the existing predictions
and frozen bins. It creates `analysis/binwise_evaluation/` without overwriting
existing outputs. The report compares seen/unseen performance within bins,
includes class-composition controls and the existing nearest-homolog baseline,
and records protein-level bootstrap intervals (2,000 replicates, seed 355).
No model is retrained and no upstream results are regenerated.
Regression checks: `python -m unittest discover -s code -p "test_*.py"`.

To recompute postprocessing from the included frozen inputs without replacing
the archived outputs, choose new direct subdirectories of `analysis/`:

```text
python code/make_similarity_bins.py --output analysis/similarity_bins_reproduced
python code/evaluate_similarity_bins.py --output analysis/binwise_evaluation_reproduced
```

Both commands refuse to overwrite an existing output directory. Evaluation
reads the included frozen predictions, not the separate retraining outputs.
The reports document the study's statistical and dataset limitations.

## 2-mer baseline

`code/run_kmer_logreg.py` trains the same logistic-regression pipeline on 400
overlapping amino-acid 2-mer frequencies (counts divided by L-1) instead of the
ESM-2 embeddings, using the same splits, training-only standardization, C grid
and validation-based selection. It needs only numpy and scikit-learn, and reads
`data/proteins_ec3.tsv` (run `git lfs pull` first). Frozen predictions and the
bin-wise evaluation are in `analysis/kmer_baseline/`; see
[its README](analysis/kmer_baseline/README.md) for the method, results, and how
the bin-wise evaluation was run. To retrain without replacing the frozen
predictions, use a new local output directory:

```text
python code/run_kmer_logreg.py --phylum Actinomycetota --out .cache/reproduction/kmer_Actinomycetota
```

Repeat with Bacillota and Cyanobacteriota for the other experiments.
