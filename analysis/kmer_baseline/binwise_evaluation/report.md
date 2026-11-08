# Bin-wise logistic-regression evaluation

## 1. Objective

Assess whether frozen ESM-2/logistic-regression performance differs between `seen_test` and `unseen_test` within the already-defined nearest-training-identity strata. This is post-processing of existing predictions, not a causal test of taxonomy.

## 2. Inputs and frozen upstream decisions

Experiments: Actinomycetota, Bacillota, Cyanobacteriota. Target: EC3. The prediction rows, split assignments, and similarity bins are read from their frozen upstream files. The numerical bins remain `<30%`, `30-40%`, `40-50%`, `50-60%`, `60-70%`, and `70-100%`; `no_hit` remains separate. No model was retrained, and no embedding, MMseqs2 search, split, filter, or bin definition was regenerated.

## 3. Validation checks

All validation checks passed. Prediction accessions are unique within each experiment; all rows have a valid group/bin and nonmissing true/predicted EC3; stored correctness matches the labels; per-accession bin and group assignments match the frozen bin and split files; and every group-by-bin count reconciles with `final_bin_counts.tsv`. See `validation.json` for counts, baseline missingness, and source hashes.

## 4. Statistical methods

Accuracy, macro-F1, and balanced accuracy are computed from the existing test predictions. Macro-F1 averages over labels present in either the true or predicted labels for that group/bin; zero F1 is assigned to a label with no true or predicted support in a bootstrap replicate. Balanced accuracy averages per-class recall over classes represented in the true labels for that evaluation sample. Shared-class macro-F1 is limited to the classes shared by both groups. The sign convention is `gap = unseen - seen` for every metric.

Uncertainty uses 2,000 ordinary protein-level bootstrap replicates with seed 355, with independent resampling within seen/unseen groups and percentile 95% intervals. Cluster bootstrap was not used: split metadata provides cluster50 IDs for seen-test proteins, but the held-out unseen phylum was intentionally not clustered, so both comparison groups cannot be resampled on a clean, common cluster mapping. Intervals do not correct for multiple comparisons.

## 5. Bin-wise sample sizes

| held_out_phylum | similarity_bin | group | n_proteins | n_ec3_classes |
| --- | --- | --- | --- | --- |
| Actinomycetota | <30% | seen_test | 178 | 48 |
| Actinomycetota | <30% | unseen_test | 461 | 69 |
| Actinomycetota | 30-40% | seen_test | 1170 | 113 |
| Actinomycetota | 30-40% | unseen_test | 1513 | 104 |
| Actinomycetota | 40-50% | seen_test | 6199 | 127 |
| Actinomycetota | 40-50% | unseen_test | 3908 | 117 |
| Actinomycetota | 50-60% | seen_test | 4199 | 84 |
| Actinomycetota | 50-60% | unseen_test | 2481 | 107 |
| Actinomycetota | 60-70% | seen_test | 887 | 47 |
| Actinomycetota | 60-70% | unseen_test | 981 | 69 |
| Actinomycetota | 70-100% | seen_test | 375 | 42 |
| Actinomycetota | 70-100% | unseen_test | 132 | 26 |
| Actinomycetota | no_hit | seen_test | 658 | 86 |
| Actinomycetota | no_hit | unseen_test | 1586 | 117 |
| Bacillota | <30% | seen_test | 268 | 51 |
| Bacillota | <30% | unseen_test | 727 | 83 |
| Bacillota | 30-40% | seen_test | 1129 | 104 |
| Bacillota | 30-40% | unseen_test | 3510 | 124 |
| Bacillota | 40-50% | seen_test | 5220 | 122 |
| Bacillota | 40-50% | unseen_test | 8392 | 124 |
| Bacillota | 50-60% | seen_test | 3543 | 84 |
| Bacillota | 50-60% | unseen_test | 8078 | 110 |
| Bacillota | 60-70% | seen_test | 503 | 45 |
| Bacillota | 60-70% | unseen_test | 3352 | 80 |
| Bacillota | 70-100% | seen_test | 375 | 38 |
| Bacillota | 70-100% | unseen_test | 815 | 49 |
| Bacillota | no_hit | seen_test | 924 | 88 |
| Bacillota | no_hit | unseen_test | 2014 | 101 |
| Cyanobacteriota | <30% | seen_test | 794 | 56 |
| Cyanobacteriota | <30% | unseen_test | 155 | 27 |
| Cyanobacteriota | 30-40% | seen_test | 1398 | 100 |
| Cyanobacteriota | 30-40% | unseen_test | 715 | 69 |
| Cyanobacteriota | 40-50% | seen_test | 5750 | 116 |
| Cyanobacteriota | 40-50% | unseen_test | 2563 | 89 |
| Cyanobacteriota | 50-60% | seen_test | 3736 | 80 |
| Cyanobacteriota | 50-60% | unseen_test | 1978 | 82 |
| Cyanobacteriota | 60-70% | seen_test | 874 | 48 |
| Cyanobacteriota | 60-70% | unseen_test | 833 | 56 |
| Cyanobacteriota | 70-100% | seen_test | 607 | 46 |
| Cyanobacteriota | 70-100% | unseen_test | 220 | 23 |
| Cyanobacteriota | no_hit | seen_test | 749 | 80 |
| Cyanobacteriota | no_hit | unseen_test | 812 | 60 |

The original binning rule required at least 100 proteins per group in every numerical bin. All bin/group counts are checked against the frozen assignments. That aggregate floor does not ensure adequate per-class support; the smallest counts and shared-class supports should be considered when interpreting macro metrics.

## 6. Raw seen-vs-unseen results

| held_out_phylum | similarity_bin | seen_n | unseen_n | seen_accuracy | unseen_accuracy | gap_accuracy | seen_macro_f1 | unseen_macro_f1 | gap_macro_f1 | seen_balanced_accuracy | unseen_balanced_accuracy | gap_balanced_accuracy |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Actinomycetota | <30% | 178 | 461 | 0.062 | 0.063 | 0.001 | 0.039 | 0.031 | -0.008 | 0.050 | 0.061 | 0.012 |
| Actinomycetota | 30-40% | 1170 | 1513 | 0.183 | 0.254 | 0.071 | 0.158 | 0.141 | -0.018 | 0.227 | 0.163 | -0.064 |
| Actinomycetota | 40-50% | 6199 | 3908 | 0.404 | 0.453 | 0.049 | 0.332 | 0.272 | -0.060 | 0.411 | 0.332 | -0.079 |
| Actinomycetota | 50-60% | 4199 | 2481 | 0.581 | 0.611 | 0.030 | 0.387 | 0.424 | 0.038 | 0.562 | 0.522 | -0.040 |
| Actinomycetota | 60-70% | 887 | 981 | 0.691 | 0.774 | 0.083 | 0.400 | 0.590 | 0.190 | 0.653 | 0.744 | 0.091 |
| Actinomycetota | 70-100% | 375 | 132 | 0.739 | 0.780 | 0.042 | 0.560 | 0.487 | -0.074 | 0.764 | 0.756 | -0.008 |
| Actinomycetota | no_hit | 658 | 1586 | 0.137 | 0.102 | -0.035 | 0.066 | 0.055 | -0.011 | 0.115 | 0.085 | -0.029 |
| Actinomycetota | all_test | 13666 | 11062 | 0.450 | 0.427 | -0.023 | 0.317 | 0.304 | -0.014 | 0.351 | 0.314 | -0.038 |
| Bacillota | <30% | 268 | 727 | 0.060 | 0.081 | 0.021 | 0.080 | 0.050 | -0.030 | 0.108 | 0.084 | -0.024 |
| Bacillota | 30-40% | 1129 | 3510 | 0.228 | 0.305 | 0.077 | 0.158 | 0.146 | -0.012 | 0.183 | 0.179 | -0.004 |
| Bacillota | 40-50% | 5220 | 8392 | 0.464 | 0.529 | 0.066 | 0.312 | 0.351 | 0.039 | 0.366 | 0.388 | 0.022 |
| Bacillota | 50-60% | 3543 | 8078 | 0.621 | 0.685 | 0.065 | 0.435 | 0.550 | 0.115 | 0.554 | 0.622 | 0.068 |
| Bacillota | 60-70% | 503 | 3352 | 0.676 | 0.765 | 0.089 | 0.497 | 0.557 | 0.060 | 0.650 | 0.671 | 0.021 |
| Bacillota | 70-100% | 375 | 815 | 0.688 | 0.843 | 0.155 | 0.564 | 0.522 | -0.042 | 0.661 | 0.721 | 0.060 |
| Bacillota | no_hit | 924 | 2014 | 0.101 | 0.127 | 0.026 | 0.056 | 0.047 | -0.009 | 0.080 | 0.076 | -0.004 |
| Bacillota | all_test | 11962 | 26888 | 0.467 | 0.543 | 0.077 | 0.280 | 0.400 | 0.120 | 0.290 | 0.393 | 0.103 |
| Cyanobacteriota | <30% | 794 | 155 | 0.029 | 0.213 | 0.184 | 0.058 | 0.082 | 0.025 | 0.094 | 0.143 | 0.049 |
| Cyanobacteriota | 30-40% | 1398 | 715 | 0.175 | 0.455 | 0.279 | 0.111 | 0.189 | 0.078 | 0.161 | 0.268 | 0.107 |
| Cyanobacteriota | 40-50% | 5750 | 2563 | 0.426 | 0.548 | 0.122 | 0.333 | 0.348 | 0.015 | 0.415 | 0.450 | 0.035 |
| Cyanobacteriota | 50-60% | 3736 | 1978 | 0.541 | 0.718 | 0.177 | 0.394 | 0.536 | 0.143 | 0.617 | 0.695 | 0.078 |
| Cyanobacteriota | 60-70% | 874 | 833 | 0.714 | 0.827 | 0.113 | 0.492 | 0.524 | 0.032 | 0.759 | 0.693 | -0.066 |
| Cyanobacteriota | 70-100% | 607 | 220 | 0.773 | 0.832 | 0.059 | 0.547 | 0.539 | -0.009 | 0.730 | 0.721 | -0.008 |
| Cyanobacteriota | no_hit | 749 | 812 | 0.095 | 0.108 | 0.014 | 0.056 | 0.029 | -0.028 | 0.087 | 0.090 | 0.002 |
| Cyanobacteriota | all_test | 13908 | 7276 | 0.424 | 0.569 | 0.145 | 0.299 | 0.435 | 0.136 | 0.336 | 0.489 | 0.153 |

`all_test` is the overall seen/unseen reference across bins and is not the primary similarity-controlled comparison.

## 7. Bootstrap uncertainty

The main accuracy estimates and unseen-minus-seen accuracy gaps, including 95% intervals, are shown in `accuracy_by_similarity_bin.png` and `accuracy_gap_by_similarity_bin.png`. All group and gap intervals for accuracy, macro-F1, and balanced accuracy—including `all_test`—are in `bootstrap_intervals.tsv`.

Across 21 experiment/bin accuracy gaps, 1 point estimates are negative; 1 intervals lie entirely below zero and 15 entirely above zero. These counts summarize the pattern and are not multiplicity-adjusted hypothesis tests.

## 8. EC3 class-composition differences

The full true-EC3 class counts and frequencies for each group/bin, together with total variation distance (TV = one half the sum of absolute frequency differences), are in `class_distribution_tv.tsv`. Larger TV indicates more different class mixtures; it is not itself a performance measure.

| held_out_phylum | similarity_bin | seen_n | unseen_n | seen_ec3_classes | unseen_ec3_classes | shared_ec3_classes | total_variation_distance |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Actinomycetota | <30% | 178 | 461 | 48 | 69 | 34 | 0.678 |
| Actinomycetota | 30-40% | 1170 | 1513 | 113 | 104 | 83 | 0.478 |
| Actinomycetota | 40-50% | 6199 | 3908 | 127 | 117 | 97 | 0.411 |
| Actinomycetota | 50-60% | 4199 | 2481 | 84 | 107 | 72 | 0.535 |
| Actinomycetota | 60-70% | 887 | 981 | 47 | 69 | 31 | 0.694 |
| Actinomycetota | 70-100% | 375 | 132 | 42 | 26 | 11 | 0.898 |
| Actinomycetota | no_hit | 658 | 1586 | 86 | 117 | 73 | 0.585 |
| Bacillota | <30% | 268 | 727 | 51 | 83 | 38 | 0.667 |
| Bacillota | 30-40% | 1129 | 3510 | 104 | 124 | 92 | 0.459 |
| Bacillota | 40-50% | 5220 | 8392 | 122 | 124 | 99 | 0.357 |
| Bacillota | 50-60% | 3543 | 8078 | 84 | 110 | 77 | 0.410 |
| Bacillota | 60-70% | 503 | 3352 | 45 | 80 | 38 | 0.620 |
| Bacillota | 70-100% | 375 | 815 | 38 | 49 | 20 | 0.735 |
| Bacillota | no_hit | 924 | 2014 | 88 | 101 | 67 | 0.610 |
| Cyanobacteriota | <30% | 794 | 155 | 56 | 27 | 16 | 0.673 |
| Cyanobacteriota | 30-40% | 1398 | 715 | 100 | 69 | 57 | 0.643 |
| Cyanobacteriota | 40-50% | 5750 | 2563 | 116 | 89 | 79 | 0.380 |
| Cyanobacteriota | 50-60% | 3736 | 1978 | 80 | 82 | 58 | 0.424 |
| Cyanobacteriota | 60-70% | 874 | 833 | 48 | 56 | 28 | 0.776 |
| Cyanobacteriota | 70-100% | 607 | 220 | 46 | 23 | 11 | 0.932 |
| Cyanobacteriota | no_hit | 749 | 812 | 80 | 60 | 44 | 0.753 |

## 9. Shared-class results

Both groups are restricted to proteins whose true EC3 class occurs in both groups within that experiment/bin. The class list is held fixed for shared-class macro-F1; out-of-shared-set predictions count as incorrect for accuracy but do not add labels to shared-class macro-F1. Small retained samples or few shared classes are weak evidence.

| held_out_phylum | similarity_bin | n_shared_classes | seen_n | unseen_n | seen_accuracy | unseen_accuracy | gap_accuracy | seen_macro_f1 | unseen_macro_f1 | gap_macro_f1 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Actinomycetota | <30% | 34 | 149 | 265 | 0.060 | 0.091 | 0.030 | 0.043 | 0.073 | 0.030 |
| Actinomycetota | 30-40% | 83 | 968 | 1450 | 0.184 | 0.259 | 0.075 | 0.203 | 0.178 | -0.026 |
| Actinomycetota | 40-50% | 97 | 5799 | 3833 | 0.411 | 0.458 | 0.047 | 0.367 | 0.341 | -0.026 |
| Actinomycetota | 50-60% | 72 | 4131 | 2132 | 0.580 | 0.608 | 0.028 | 0.487 | 0.501 | 0.014 |
| Actinomycetota | 60-70% | 31 | 744 | 615 | 0.677 | 0.707 | 0.030 | 0.623 | 0.674 | 0.051 |
| Actinomycetota | 70-100% | 11 | 79 | 43 | 0.709 | 0.651 | -0.058 | 0.689 | 0.655 | -0.034 |
| Actinomycetota | no_hit | 73 | 622 | 1353 | 0.138 | 0.114 | -0.024 | 0.086 | 0.089 | 0.003 |
| Bacillota | <30% | 38 | 235 | 489 | 0.043 | 0.110 | 0.068 | 0.091 | 0.090 | -0.001 |
| Bacillota | 30-40% | 92 | 1038 | 3260 | 0.247 | 0.319 | 0.073 | 0.198 | 0.178 | -0.020 |
| Bacillota | 40-50% | 99 | 5059 | 8168 | 0.467 | 0.534 | 0.067 | 0.362 | 0.430 | 0.069 |
| Bacillota | 50-60% | 77 | 3459 | 7403 | 0.628 | 0.686 | 0.057 | 0.531 | 0.618 | 0.087 |
| Bacillota | 60-70% | 38 | 483 | 2266 | 0.675 | 0.725 | 0.050 | 0.610 | 0.663 | 0.052 |
| Bacillota | 70-100% | 20 | 257 | 384 | 0.646 | 0.854 | 0.208 | 0.506 | 0.701 | 0.195 |
| Bacillota | no_hit | 67 | 636 | 1838 | 0.124 | 0.137 | 0.012 | 0.082 | 0.062 | -0.020 |
| Cyanobacteriota | <30% | 16 | 597 | 92 | 0.015 | 0.141 | 0.126 | 0.130 | 0.121 | -0.009 |
| Cyanobacteriota | 30-40% | 57 | 944 | 673 | 0.224 | 0.478 | 0.255 | 0.180 | 0.275 | 0.095 |
| Cyanobacteriota | 40-50% | 79 | 5349 | 2495 | 0.426 | 0.554 | 0.128 | 0.423 | 0.447 | 0.025 |
| Cyanobacteriota | 50-60% | 58 | 3364 | 1748 | 0.526 | 0.716 | 0.190 | 0.564 | 0.720 | 0.156 |
| Cyanobacteriota | 60-70% | 28 | 639 | 443 | 0.684 | 0.806 | 0.122 | 0.695 | 0.669 | -0.026 |
| Cyanobacteriota | 70-100% | 11 | 56 | 90 | 0.375 | 0.711 | 0.336 | 0.465 | 0.601 | 0.136 |
| Cyanobacteriota | no_hit | 44 | 507 | 580 | 0.130 | 0.148 | 0.018 | 0.112 | 0.070 | -0.042 |

## 10. Class-reweighted results

Only classes present in both groups are retained. For class c, each seen observation receives per-observation weight `w_i = (n_unseen,c / N_unseen,common) / n_seen,c`. Thus all seen observations in class c carry total weight equal to the unseen class proportion, the weights sum to one, and the seen class mix matches unseen over common support. Unseen accuracy is computed on the same common-support subset. Classes absent from either group cannot be matched by weighting; their protein counts and class counts are reported explicitly.

| held_out_phylum | similarity_bin | common_support_class_count | seen_n | unseen_n | seen_excluded_n | unseen_excluded_n | class_mix_adjusted_seen_accuracy | unseen_common_support_accuracy | adjusted_gap_unseen_minus_seen |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Actinomycetota | <30% | 34 | 178 | 461 | 29 | 196 | 0.029 | 0.091 | 0.062 |
| Actinomycetota | 30-40% | 83 | 1170 | 1513 | 202 | 63 | 0.289 | 0.259 | -0.029 |
| Actinomycetota | 40-50% | 97 | 6199 | 3908 | 400 | 75 | 0.468 | 0.458 | -0.010 |
| Actinomycetota | 50-60% | 72 | 4199 | 2481 | 68 | 349 | 0.571 | 0.608 | 0.038 |
| Actinomycetota | 60-70% | 31 | 887 | 981 | 143 | 366 | 0.571 | 0.707 | 0.136 |
| Actinomycetota | 70-100% | 11 | 375 | 132 | 296 | 89 | 0.199 | 0.651 | 0.452 |
| Actinomycetota | no_hit | 73 | 658 | 1586 | 36 | 233 | 0.179 | 0.114 | -0.065 |
| Bacillota | <30% | 38 | 268 | 727 | 33 | 238 | 0.043 | 0.110 | 0.068 |
| Bacillota | 30-40% | 92 | 1129 | 3510 | 91 | 250 | 0.284 | 0.319 | 0.035 |
| Bacillota | 40-50% | 99 | 5220 | 8392 | 161 | 224 | 0.487 | 0.534 | 0.046 |
| Bacillota | 50-60% | 77 | 3543 | 8078 | 84 | 675 | 0.634 | 0.686 | 0.051 |
| Bacillota | 60-70% | 38 | 503 | 3352 | 20 | 1086 | 0.785 | 0.725 | -0.060 |
| Bacillota | 70-100% | 20 | 375 | 815 | 118 | 431 | 0.371 | 0.854 | 0.483 |
| Bacillota | no_hit | 67 | 924 | 2014 | 288 | 176 | 0.120 | 0.137 | 0.016 |
| Cyanobacteriota | <30% | 16 | 794 | 155 | 197 | 63 | 0.129 | 0.141 | 0.012 |
| Cyanobacteriota | 30-40% | 57 | 1398 | 715 | 454 | 42 | 0.348 | 0.478 | 0.130 |
| Cyanobacteriota | 40-50% | 79 | 5750 | 2563 | 401 | 68 | 0.489 | 0.554 | 0.066 |
| Cyanobacteriota | 50-60% | 58 | 3736 | 1978 | 372 | 230 | 0.526 | 0.716 | 0.189 |
| Cyanobacteriota | 60-70% | 28 | 874 | 833 | 235 | 390 | 0.735 | 0.806 | 0.070 |
| Cyanobacteriota | 70-100% | 11 | 607 | 220 | 551 | 130 | 0.313 | 0.711 | 0.398 |
| Cyanobacteriota | no_hit | 44 | 749 | 812 | 242 | 232 | 0.175 | 0.148 | -0.027 |

## 11. `no_hit` results

`no_hit` is kept as its own category and is not assigned identity zero. It means no MMseqs2 training alignment passed the qualifying search criteria; it does not mean zero sequence identity. Its counts and model performance are included separately in all bin-wise tables and plots.

## 12. Lab-tested sensitivity analysis

This repeats bin-wise accuracy for rows with `lab_tested=yes`; it is a sensitivity analysis, not the primary result. Wilson 95% intervals are descriptive. Counts below 100 are flagged as small-support cells and should not be overinterpreted.

| held_out_phylum | similarity_bin | group | n_proteins | accuracy | accuracy_ci_low | accuracy_ci_high |
| --- | --- | --- | --- | --- | --- | --- |
| Actinomycetota | <30% | seen_test | 41 | 0.073 | 0.025 | 0.194 |
| Actinomycetota | <30% | unseen_test | 179 | 0.034 | 0.015 | 0.071 |
| Actinomycetota | 30-40% | seen_test | 152 | 0.151 | 0.103 | 0.217 |
| Actinomycetota | 30-40% | unseen_test | 271 | 0.137 | 0.101 | 0.183 |
| Actinomycetota | 40-50% | seen_test | 146 | 0.192 | 0.136 | 0.263 |
| Actinomycetota | 40-50% | unseen_test | 181 | 0.271 | 0.211 | 0.340 |
| Actinomycetota | 50-60% | seen_test | 23 | 0.391 | 0.222 | 0.592 |
| Actinomycetota | 50-60% | unseen_test | 85 | 0.376 | 0.281 | 0.483 |
| Actinomycetota | 60-70% | seen_test | 5 | 0.600 | 0.231 | 0.882 |
| Actinomycetota | 60-70% | unseen_test | 34 | 0.500 | 0.341 | 0.659 |
| Actinomycetota | 70-100% | seen_test | 6 | 0.500 | 0.188 | 0.812 |
| Actinomycetota | 70-100% | unseen_test | 10 | 0.400 | 0.168 | 0.687 |
| Actinomycetota | no_hit | seen_test | 108 | 0.074 | 0.038 | 0.139 |
| Actinomycetota | no_hit | unseen_test | 433 | 0.090 | 0.067 | 0.121 |
| Actinomycetota | all_test | seen_test | 481 | 0.160 | 0.130 | 0.196 |
| Actinomycetota | all_test | unseen_test | 1193 | 0.154 | 0.135 | 0.176 |
| Bacillota | <30% | seen_test | 51 | 0.157 | 0.082 | 0.280 |
| Bacillota | <30% | unseen_test | 101 | 0.089 | 0.048 | 0.161 |
| Bacillota | 30-40% | seen_test | 139 | 0.158 | 0.107 | 0.228 |
| Bacillota | 30-40% | unseen_test | 237 | 0.143 | 0.105 | 0.194 |
| Bacillota | 40-50% | seen_test | 142 | 0.211 | 0.152 | 0.286 |
| Bacillota | 40-50% | unseen_test | 214 | 0.220 | 0.169 | 0.280 |
| Bacillota | 50-60% | seen_test | 25 | 0.520 | 0.335 | 0.700 |
| Bacillota | 50-60% | unseen_test | 119 | 0.395 | 0.312 | 0.485 |
| Bacillota | 60-70% | seen_test | 2 | 0.500 | 0.095 | 0.905 |
| Bacillota | 60-70% | unseen_test | 62 | 0.516 | 0.394 | 0.636 |
| Bacillota | 70-100% | seen_test | 5 | 0.400 | 0.118 | 0.769 |
| Bacillota | 70-100% | unseen_test | 21 | 0.619 | 0.409 | 0.792 |
| Bacillota | no_hit | seen_test | 127 | 0.031 | 0.012 | 0.078 |
| Bacillota | no_hit | unseen_test | 328 | 0.085 | 0.060 | 0.121 |
| Bacillota | all_test | seen_test | 491 | 0.163 | 0.133 | 0.198 |
| Bacillota | all_test | unseen_test | 1082 | 0.194 | 0.172 | 0.219 |
| Cyanobacteriota | <30% | seen_test | 48 | 0.062 | 0.021 | 0.168 |
| Cyanobacteriota | <30% | unseen_test | 13 | 0.000 | 0.000 | 0.228 |
| Cyanobacteriota | 30-40% | seen_test | 122 | 0.074 | 0.039 | 0.134 |
| Cyanobacteriota | 30-40% | unseen_test | 21 | 0.143 | 0.050 | 0.346 |
| Cyanobacteriota | 40-50% | seen_test | 161 | 0.248 | 0.188 | 0.321 |
| Cyanobacteriota | 40-50% | unseen_test | 22 | 0.273 | 0.132 | 0.482 |
| Cyanobacteriota | 50-60% | seen_test | 39 | 0.231 | 0.126 | 0.383 |
| Cyanobacteriota | 50-60% | unseen_test | 25 | 0.360 | 0.202 | 0.555 |
| Cyanobacteriota | 60-70% | seen_test | 12 | 0.667 | 0.391 | 0.862 |
| Cyanobacteriota | 60-70% | unseen_test | 4 | 0.750 | 0.301 | 0.954 |
| Cyanobacteriota | 70-100% | seen_test | 2 | 0.000 | 0.000 | 0.658 |
| Cyanobacteriota | 70-100% | unseen_test | 5 | 0.600 | 0.231 | 0.882 |
| Cyanobacteriota | no_hit | seen_test | 117 | 0.085 | 0.047 | 0.150 |
| Cyanobacteriota | no_hit | unseen_test | 56 | 0.036 | 0.010 | 0.121 |
| Cyanobacteriota | all_test | seen_test | 501 | 0.158 | 0.128 | 0.192 |
| Cyanobacteriota | all_test | unseen_test | 146 | 0.178 | 0.125 | 0.248 |

Cells with fewer than 100 lab-tested proteins: Actinomycetota/<30%/seen_test (n=41), Actinomycetota/50-60%/seen_test (n=23), Actinomycetota/50-60%/unseen_test (n=85), Actinomycetota/60-70%/seen_test (n=5), Actinomycetota/60-70%/unseen_test (n=34), Actinomycetota/70-100%/seen_test (n=6), Actinomycetota/70-100%/unseen_test (n=10), Bacillota/<30%/seen_test (n=51), Bacillota/50-60%/seen_test (n=25), Bacillota/60-70%/seen_test (n=2), Bacillota/60-70%/unseen_test (n=62), Bacillota/70-100%/seen_test (n=5), Bacillota/70-100%/unseen_test (n=21), Cyanobacteriota/<30%/seen_test (n=48), Cyanobacteriota/<30%/unseen_test (n=13), Cyanobacteriota/30-40%/unseen_test (n=21), Cyanobacteriota/40-50%/unseen_test (n=22), Cyanobacteriota/50-60%/seen_test (n=39), Cyanobacteriota/50-60%/unseen_test (n=25), Cyanobacteriota/60-70%/seen_test (n=12), Cyanobacteriota/60-70%/unseen_test (n=4), Cyanobacteriota/70-100%/seen_test (n=2), Cyanobacteriota/70-100%/unseen_test (n=5), Cyanobacteriota/no_hit/unseen_test (n=56)

## 13. Nearest-homolog baseline

Baseline 1 is the existing nearest-training-accession EC3 label transfer; no baseline model was created. `accuracy_all_test_missing_as_wrong` counts missing predictions as incorrect, consistent with the upstream logistic-regression README. Available-only accuracy is also shown. Missing values are counted by bin and their reason is in `validation.json`.

| held_out_phylum | similarity_bin | group | n_proteins | n_baseline_predictions | n_missing_baseline_predictions | accuracy_all_test_missing_as_wrong | accuracy_when_baseline_available |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Actinomycetota | all_test | seen_test | 13666 | 13008 | 658 | 0.941 | 0.989 |
| Actinomycetota | all_test | unseen_test | 11062 | 9476 | 1586 | 0.823 | 0.961 |
| Bacillota | all_test | seen_test | 11962 | 11038 | 924 | 0.897 | 0.973 |
| Bacillota | all_test | unseen_test | 26888 | 24874 | 2014 | 0.911 | 0.984 |
| Cyanobacteriota | all_test | seen_test | 13908 | 13159 | 749 | 0.900 | 0.951 |
| Cyanobacteriota | all_test | unseen_test | 7276 | 6464 | 812 | 0.880 | 0.990 |

Across numerical-bin/group cells, nearest-homolog accuracy is often near ceiling at high identity and exceeds the logistic model in many hit-bearing cells (unweighted cell means: 0.934 vs 0.494). The logistic model improves on nearest-neighbor transfer in some low-identity cells and on every `no_hit` cell; baseline accuracy is zero for all 6,743 `no_hit` proteins because no donor prediction exists. The available-only baseline accuracy omits these rows and must not be compared to all-test logistic accuracy without noting its reduced coverage.

## 14. Interpretation

Overall seen/unseen accuracy gaps are -0.023 for Actinomycetota, +0.077 for Bacillota, and +0.145 for Cyanobacteriota. Within strata, 1 of 21 accuracy gaps are negative; 1 intervals lie wholly below zero and 15 wholly above zero. The direction is not consistent across phyla or bins:

- **Actinomycetota:** 0/6 numerical-bin accuracy gaps are negative (range +0.001 to +0.083); the separate `no_hit` gap is -0.035.
- **Bacillota:** 0/6 numerical-bin accuracy gaps are negative (range +0.021 to +0.155); the separate `no_hit` gap is +0.026.
- **Cyanobacteriota:** 0/6 numerical-bin accuracy gaps are negative (range +0.059 to +0.279); the separate `no_hit` gap is +0.014.

Some 30-50% bins favor unseen, while other matched bins favor seen. The `no_hit` direction also reverses across phyla. Class reweighting shrinks several differences toward zero, but residual gaps remain in both directions and the `no_hit` results change materially. Shared-class and reweighted analyses therefore do not support a single consistent unseen penalty. Similarity bins are broad categorical controls, not exact matching, and can retain within-bin identity differences. EC3 composition differs substantially across many strata. Even a persistent gap is an association under this evaluation design, not proof that taxonomy causes it.

## 15. Limitations

- Protein-level resampling does not account for dependence among homologous sequences; comparable sequence-cluster IDs are unavailable for the unseen group.
- Similarity bins are coarse, including the broad 70-100% interval; within-bin similarity may differ.
- Macro-F1 and balanced accuracy are sensitive to class support and composition; do not use macro-F1 alone to infer a taxonomic effect.
- Bootstrap intervals are descriptive and not adjusted for multiple bins, metrics, or experiments. Shared-class/reweighted analyses can exclude substantial numbers of proteins.
- `no_hit` reflects a search qualification outcome, not absence of biological homology.
- The lab-tested flag is analyzed as supplied and does not guarantee independent experimental validation for every label.

## 16. Exact output files generated

- `report.md`
- `binwise_metrics.tsv`
- `binwise_gaps.tsv`
- `class_distribution_tv.tsv`
- `shared_class_metrics.tsv`
- `reweighted_metrics.tsv`
- `bootstrap_intervals.tsv`
- `lab_tested_metrics.tsv`
- `baseline1_binwise_metrics.tsv`
- `validation.json`
- `settings.json`
- `accuracy_by_similarity_bin.png`
- `accuracy_by_similarity_bin.pdf`
- `accuracy_gap_by_similarity_bin.png`
- `accuracy_gap_by_similarity_bin.pdf`
