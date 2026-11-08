# k-mer baseline (proposal baseline 1)

Amino-acid 2-mer frequencies + multinomial logistic regression, for the same three held-out phyla
and evaluated exactly like the ESM-2 model. Only the input features differ from
`analysis/logistic_regression/`.

## Method

- **Features:** 400 overlapping 2-mer frequencies per protein, counts / (L-1) with L the full
  sequence length. 2-mers containing a letter outside the 20 standard amino acids are skipped,
  while the divisor stays L-1. This affects 52 of 146,677 proteins (51 with U, 1 with O).
- **Data and splits:** `data/proteins_ec3.tsv` (146,677 rows), `splits_ec3/<phylum>/split_<phylum>.tsv`
  and `analysis/similarity_bins/<phylum>_similarity_bins.tsv`, all unchanged. The group sizes and
  the test proteins, bins and labels are identical to the ESM-2 run, checked row by row.
- **Model:** `code/run_kmer_logreg.py` keeps `code/run_logreg.py`'s split loading, training-only
  standardization, `LogisticRegression(solver="lbfgs", max_iter=300)`, C in {0.001, 0.01, 0.1, 1}
  chosen by validation accuracy, no refit on train+validation, and the output format. There were
  no convergence warnings (at most 267 iterations).
- **Evaluation:** `binwise_evaluation/` is the output of `code/evaluate_similarity_bins.py`, run
  unchanged. It used 2,000 protein-level resamples, seed 355, and percentile 95% intervals. The
  script reads predictions from fixed paths, so it was run on a copy of this repository in which
  `analysis/logistic_regression/<phylum>/predictions.tsv` was replaced by the k-mer predictions.
  The input hashes in `binwise_evaluation/settings.json` and `validation.json` therefore carry
  the `logistic_regression` path names but refer to the k-mer files.
- **Sanity check:** the same procedure on the ESM-2 predictions reproduces
  `analysis/binwise_evaluation/` (for example, Actinomycetota overall gap -3.62 [-4.24, -3.00]).
  A separate script without pandas recomputed every point gap and matched exactly.

Reproduce from the repository root:

```
python code/run_kmer_logreg.py --phylum Actinomycetota
python code/run_kmer_logreg.py --phylum Bacillota
python code/run_kmer_logreg.py --phylum Cyanobacteriota
```

## Overall test accuracy (%)

| Held-out phylum | Group | k-mer + LR | ESM-2 + LR | Nearest homolog (all records) |
|---|---|---:|---:|---:|
| Actinomycetota | Seen | 44.99 | 95.29 | 94.10 |
|  | Unseen | 42.69 | 91.67 | 82.33 |
| Bacillota | Seen | 46.68 | 94.65 | 89.74 |
|  | Unseen | 54.34 | 94.76 | 91.07 |
| Cyanobacteriota | Seen | 42.44 | 92.12 | 89.98 |
|  | Unseen | 56.94 | 91.93 | 87.99 |

Always predicting the most common training class scores 0-7% in every group.

## Unseen-minus-seen accuracy gap by identity bin (pp, 95% bootstrap CI)

| Bin | Actinomycetota k-mer | ESM-2 | Bacillota k-mer | ESM-2 | Cyanobacteriota k-mer | ESM-2 |
|---|---:|---:|---:|---:|---:|---:|
| <30% | +0.11 [-4.35, +4.05] | -2.51 | +2.15 [-1.66, +5.37] | +0.95 | +18.39 [+12.19, +25.24] | +27.41 |
| 30-40% | +7.09 [+3.97, +10.16] | -5.76 | +7.69 [+4.87, +10.49] | +10.05 | +27.93 [+23.63, +32.18] | +4.41 |
| 40-50% | +4.90 [+2.92, +6.85] | -0.76 | +6.56 [+4.91, +8.27] | +1.61 | +12.21 [+9.80, +14.37] | -0.27 |
| 50-60% | +2.96 [+0.51, +5.44] | -0.38 | +6.47 [+4.58, +8.31] | +0.27 | +17.69 [+15.11, +20.17] | +0.17 |
| 60-70% | +8.26 [+4.25, +12.20] | -0.71 | +8.87 [+4.70, +13.20] | -0.39 | +11.32 [+7.29, +15.05] | 0.00 |
| 70-100% | +4.16 [-4.25, +12.49] | -0.76 | +15.49 [+10.26, +20.62] | -0.25 | +5.92 [-0.15, +11.70] | 0.00 |
| no_hit | -3.46 [-6.67, -0.54] | +26.37 | +2.60 [+0.17, +4.97] | -18.25 | +1.36 [-1.88, +4.38] | -6.43 |
| all test | -2.31 [-3.52, -1.06] | -3.62 | +7.66 [+6.62, +8.73] | +0.11 | +14.50 [+13.06, +15.95] | -0.19 |

Per-bin accuracies, macro-F1, balanced accuracy, the class-mix reweighting and the lab-tested
subsets are in `binwise_evaluation/`.

## What it shows

- The 2-mer model is much less accurate than ESM-2, and its accuracy rises steeply with identity
  to the nearest training protein: 3-21% below 30% identity, 69-84% at 70-100%.
- With ESM-2, Actinomycetota unseen accuracy was lower than seen in all 6 identity bins. With
  2-mers it is higher in all 6. The overall 2-mer gap is still negative, because more unseen
  proteins fall in `no_hit` and the low bins, where the model does badly.
- In Bacillota and Cyanobacteriota, 2-mer unseen accuracy is higher in every bin.
- Class-mix reweighting (`binwise_evaluation/reweighted_metrics.tsv`) explains part of this. For
  example, Actinomycetota 30-40% goes from +7.09 to -2.92 pp.

## Caveats

- Below 30% identity and in `no_hit`, 2-mer accuracy is near chance, so gaps there say little.
- For Actinomycetota the chosen C (1) is the largest in the grid. Validation accuracy rose only
  0.35 pp from C=0.1, and the grid was kept identical to ESM-2.
- Only k=2 was tried.
- Why unseen accuracy is higher with 2-mers is not established beyond the class-mix part.
- The intervals are protein-level bootstrap intervals with no cluster resampling and no
  multiple-comparison adjustment, as in `analysis/binwise_evaluation/`.
