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
| Actinomycetota | <30% | 178 | 461 | 0.787 | 0.761 | -0.025 | 0.611 | 0.514 | -0.097 | 0.716 | 0.629 | -0.087 |
| Actinomycetota | 30-40% | 1170 | 1513 | 0.919 | 0.861 | -0.058 | 0.735 | 0.660 | -0.076 | 0.802 | 0.734 | -0.069 |
| Actinomycetota | 40-50% | 6199 | 3908 | 0.981 | 0.974 | -0.008 | 0.901 | 0.870 | -0.031 | 0.935 | 0.907 | -0.028 |
| Actinomycetota | 50-60% | 4199 | 2481 | 0.993 | 0.990 | -0.004 | 0.970 | 0.942 | -0.028 | 0.979 | 0.968 | -0.011 |
| Actinomycetota | 60-70% | 887 | 981 | 1.000 | 0.993 | -0.007 | 1.000 | 0.927 | -0.073 | 1.000 | 0.942 | -0.058 |
| Actinomycetota | 70-100% | 375 | 132 | 1.000 | 0.992 | -0.008 | 1.000 | 0.961 | -0.039 | 1.000 | 0.962 | -0.038 |
| Actinomycetota | no_hit | 658 | 1586 | 0.444 | 0.707 | 0.264 | 0.479 | 0.412 | -0.067 | 0.571 | 0.466 | -0.105 |
| Actinomycetota | all_test | 13666 | 11062 | 0.953 | 0.917 | -0.036 | 0.794 | 0.739 | -0.055 | 0.816 | 0.759 | -0.057 |
| Bacillota | <30% | 268 | 727 | 0.791 | 0.801 | 0.010 | 0.666 | 0.480 | -0.185 | 0.802 | 0.584 | -0.218 |
| Bacillota | 30-40% | 1129 | 3510 | 0.808 | 0.908 | 0.100 | 0.738 | 0.697 | -0.041 | 0.803 | 0.762 | -0.041 |
| Bacillota | 40-50% | 5220 | 8392 | 0.971 | 0.987 | 0.016 | 0.873 | 0.917 | 0.043 | 0.896 | 0.942 | 0.046 |
| Bacillota | 50-60% | 3543 | 8078 | 0.994 | 0.996 | 0.003 | 0.971 | 0.934 | -0.038 | 0.973 | 0.959 | -0.013 |
| Bacillota | 60-70% | 503 | 3352 | 1.000 | 0.996 | -0.004 | 1.000 | 0.959 | -0.041 | 1.000 | 0.967 | -0.033 |
| Bacillota | 70-100% | 375 | 815 | 1.000 | 0.998 | -0.002 | 1.000 | 0.957 | -0.043 | 1.000 | 0.959 | -0.041 |
| Bacillota | no_hit | 924 | 2014 | 0.792 | 0.610 | -0.182 | 0.525 | 0.370 | -0.155 | 0.601 | 0.473 | -0.128 |
| Bacillota | all_test | 11962 | 26888 | 0.946 | 0.948 | 0.001 | 0.799 | 0.772 | -0.027 | 0.828 | 0.811 | -0.017 |
| Cyanobacteriota | <30% | 794 | 155 | 0.339 | 0.613 | 0.274 | 0.571 | 0.522 | -0.049 | 0.727 | 0.719 | -0.008 |
| Cyanobacteriota | 30-40% | 1398 | 715 | 0.935 | 0.979 | 0.044 | 0.794 | 0.857 | 0.063 | 0.868 | 0.922 | 0.053 |
| Cyanobacteriota | 40-50% | 5750 | 2563 | 0.984 | 0.981 | -0.003 | 0.909 | 0.975 | 0.066 | 0.936 | 0.987 | 0.051 |
| Cyanobacteriota | 50-60% | 3736 | 1978 | 0.997 | 0.999 | 0.002 | 0.972 | 0.984 | 0.012 | 0.973 | 0.984 | 0.010 |
| Cyanobacteriota | 60-70% | 874 | 833 | 1.000 | 1.000 | 0.000 | 1.000 | 1.000 | 0.000 | 1.000 | 1.000 | 0.000 |
| Cyanobacteriota | 70-100% | 607 | 220 | 1.000 | 1.000 | 0.000 | 1.000 | 1.000 | 0.000 | 1.000 | 1.000 | 0.000 |
| Cyanobacteriota | no_hit | 749 | 812 | 0.495 | 0.431 | -0.064 | 0.501 | 0.314 | -0.187 | 0.620 | 0.514 | -0.106 |
| Cyanobacteriota | all_test | 13908 | 7276 | 0.921 | 0.919 | -0.002 | 0.776 | 0.792 | 0.016 | 0.819 | 0.862 | 0.042 |

`all_test` is the overall seen/unseen reference across bins and is not the primary similarity-controlled comparison.

## 7. Bootstrap uncertainty

The main accuracy estimates and unseen-minus-seen accuracy gaps, including 95% intervals, are shown in `accuracy_by_similarity_bin.png` and `accuracy_gap_by_similarity_bin.png`. All group and gap intervals for accuracy, macro-F1, and balanced accuracy—including `all_test`—are in `bootstrap_intervals.tsv`.

Across 21 experiment/bin accuracy gaps, 11 point estimates are negative; 6 intervals lie entirely below zero and 5 entirely above zero. These counts summarize the pattern and are not multiplicity-adjusted hypothesis tests.

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
| Actinomycetota | <30% | 34 | 149 | 265 | 0.839 | 0.834 | -0.005 | 0.777 | 0.721 | -0.056 |
| Actinomycetota | 30-40% | 83 | 968 | 1450 | 0.916 | 0.874 | -0.043 | 0.805 | 0.760 | -0.045 |
| Actinomycetota | 40-50% | 97 | 5799 | 3833 | 0.981 | 0.976 | -0.004 | 0.933 | 0.942 | 0.009 |
| Actinomycetota | 50-60% | 72 | 4131 | 2132 | 0.993 | 0.991 | -0.002 | 0.965 | 0.970 | 0.005 |
| Actinomycetota | 60-70% | 31 | 744 | 615 | 1.000 | 0.997 | -0.003 | 1.000 | 0.964 | -0.036 |
| Actinomycetota | 70-100% | 11 | 79 | 43 | 1.000 | 1.000 | 0.000 | 1.000 | 1.000 | 0.000 |
| Actinomycetota | no_hit | 73 | 622 | 1353 | 0.434 | 0.758 | 0.323 | 0.563 | 0.548 | -0.015 |
| Bacillota | <30% | 38 | 235 | 489 | 0.796 | 0.883 | 0.088 | 0.850 | 0.691 | -0.158 |
| Bacillota | 30-40% | 92 | 1038 | 3260 | 0.798 | 0.933 | 0.135 | 0.786 | 0.786 | 0.001 |
| Bacillota | 40-50% | 99 | 5059 | 8168 | 0.971 | 0.988 | 0.017 | 0.915 | 0.964 | 0.049 |
| Bacillota | 50-60% | 77 | 3459 | 7403 | 0.993 | 0.999 | 0.005 | 0.969 | 0.992 | 0.023 |
| Bacillota | 60-70% | 38 | 483 | 2266 | 1.000 | 0.998 | -0.002 | 1.000 | 0.997 | -0.003 |
| Bacillota | 70-100% | 20 | 257 | 384 | 1.000 | 1.000 | 0.000 | 1.000 | 1.000 | 0.000 |
| Bacillota | no_hit | 67 | 636 | 1838 | 0.737 | 0.650 | -0.088 | 0.582 | 0.504 | -0.078 |
| Cyanobacteriota | <30% | 16 | 597 | 92 | 0.151 | 0.543 | 0.393 | 0.648 | 0.701 | 0.053 |
| Cyanobacteriota | 30-40% | 57 | 944 | 673 | 0.928 | 0.981 | 0.053 | 0.895 | 0.921 | 0.026 |
| Cyanobacteriota | 40-50% | 79 | 5349 | 2495 | 0.985 | 0.981 | -0.004 | 0.954 | 0.972 | 0.018 |
| Cyanobacteriota | 50-60% | 58 | 3364 | 1748 | 0.998 | 1.000 | 0.002 | 0.998 | 1.000 | 0.002 |
| Cyanobacteriota | 60-70% | 28 | 639 | 443 | 1.000 | 1.000 | 0.000 | 1.000 | 1.000 | 0.000 |
| Cyanobacteriota | 70-100% | 11 | 56 | 90 | 1.000 | 1.000 | 0.000 | 1.000 | 1.000 | 0.000 |
| Cyanobacteriota | no_hit | 44 | 507 | 580 | 0.473 | 0.579 | 0.106 | 0.692 | 0.501 | -0.192 |

## 10. Class-reweighted results

Only classes present in both groups are retained. For class c, each seen observation receives per-observation weight `w_i = (n_unseen,c / N_unseen,common) / n_seen,c`. Thus all seen observations in class c carry total weight equal to the unseen class proportion, the weights sum to one, and the seen class mix matches unseen over common support. Unseen accuracy is computed on the same common-support subset. Classes absent from either group cannot be matched by weighting; their protein counts and class counts are reported explicitly.

| held_out_phylum | similarity_bin | common_support_class_count | seen_n | unseen_n | seen_excluded_n | unseen_excluded_n | class_mix_adjusted_seen_accuracy | unseen_common_support_accuracy | adjusted_gap_unseen_minus_seen |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Actinomycetota | <30% | 34 | 178 | 461 | 29 | 196 | 0.914 | 0.834 | -0.080 |
| Actinomycetota | 30-40% | 83 | 1170 | 1513 | 202 | 63 | 0.901 | 0.874 | -0.027 |
| Actinomycetota | 40-50% | 97 | 6199 | 3908 | 400 | 75 | 0.976 | 0.976 | -0.000 |
| Actinomycetota | 50-60% | 72 | 4199 | 2481 | 68 | 349 | 0.994 | 0.991 | -0.003 |
| Actinomycetota | 60-70% | 31 | 887 | 981 | 143 | 366 | 1.000 | 0.997 | -0.003 |
| Actinomycetota | 70-100% | 11 | 375 | 132 | 296 | 89 | 1.000 | 1.000 | 0.000 |
| Actinomycetota | no_hit | 73 | 658 | 1586 | 36 | 233 | 0.730 | 0.758 | 0.027 |
| Bacillota | <30% | 38 | 268 | 727 | 33 | 238 | 0.924 | 0.883 | -0.040 |
| Bacillota | 30-40% | 92 | 1129 | 3510 | 91 | 250 | 0.893 | 0.933 | 0.040 |
| Bacillota | 40-50% | 99 | 5220 | 8392 | 161 | 224 | 0.977 | 0.988 | 0.011 |
| Bacillota | 50-60% | 77 | 3543 | 8078 | 84 | 675 | 0.991 | 0.999 | 0.008 |
| Bacillota | 60-70% | 38 | 503 | 3352 | 20 | 1086 | 1.000 | 0.998 | -0.002 |
| Bacillota | 70-100% | 20 | 375 | 815 | 118 | 431 | 1.000 | 1.000 | 0.000 |
| Bacillota | no_hit | 67 | 924 | 2014 | 288 | 176 | 0.771 | 0.650 | -0.122 |
| Cyanobacteriota | <30% | 16 | 794 | 155 | 197 | 63 | 0.561 | 0.543 | -0.017 |
| Cyanobacteriota | 30-40% | 57 | 1398 | 715 | 454 | 42 | 0.975 | 0.981 | 0.006 |
| Cyanobacteriota | 40-50% | 79 | 5750 | 2563 | 401 | 68 | 0.969 | 0.981 | 0.011 |
| Cyanobacteriota | 50-60% | 58 | 3736 | 1978 | 372 | 230 | 0.998 | 1.000 | 0.002 |
| Cyanobacteriota | 60-70% | 28 | 874 | 833 | 235 | 390 | 1.000 | 1.000 | 0.000 |
| Cyanobacteriota | 70-100% | 11 | 607 | 220 | 551 | 130 | 1.000 | 1.000 | 0.000 |
| Cyanobacteriota | no_hit | 44 | 749 | 812 | 242 | 232 | 0.814 | 0.579 | -0.235 |

## 11. `no_hit` results

`no_hit` is kept as its own category and is not assigned identity zero. It means no MMseqs2 training alignment passed the qualifying search criteria; it does not mean zero sequence identity. Its counts and model performance are included separately in all bin-wise tables and plots.

## 12. Lab-tested sensitivity analysis

This repeats bin-wise accuracy for rows with `lab_tested=yes`; it is a sensitivity analysis, not the primary result. Wilson 95% intervals are descriptive. Counts below 100 are flagged as small-support cells and should not be overinterpreted.

| held_out_phylum | similarity_bin | group | n_proteins | accuracy | accuracy_ci_low | accuracy_ci_high |
| --- | --- | --- | --- | --- | --- | --- |
| Actinomycetota | <30% | seen_test | 41 | 0.659 | 0.505 | 0.784 |
| Actinomycetota | <30% | unseen_test | 179 | 0.642 | 0.570 | 0.709 |
| Actinomycetota | 30-40% | seen_test | 152 | 0.822 | 0.754 | 0.875 |
| Actinomycetota | 30-40% | unseen_test | 271 | 0.764 | 0.710 | 0.810 |
| Actinomycetota | 40-50% | seen_test | 146 | 0.863 | 0.798 | 0.910 |
| Actinomycetota | 40-50% | unseen_test | 181 | 0.840 | 0.779 | 0.886 |
| Actinomycetota | 50-60% | seen_test | 23 | 1.000 | 0.857 | 1.000 |
| Actinomycetota | 50-60% | unseen_test | 85 | 0.906 | 0.825 | 0.952 |
| Actinomycetota | 60-70% | seen_test | 5 | 1.000 | 0.566 | 1.000 |
| Actinomycetota | 60-70% | unseen_test | 34 | 0.941 | 0.809 | 0.984 |
| Actinomycetota | 70-100% | seen_test | 6 | 1.000 | 0.610 | 1.000 |
| Actinomycetota | 70-100% | unseen_test | 10 | 0.900 | 0.596 | 0.982 |
| Actinomycetota | no_hit | seen_test | 108 | 0.574 | 0.480 | 0.663 |
| Actinomycetota | no_hit | unseen_test | 433 | 0.554 | 0.507 | 0.600 |
| Actinomycetota | all_test | seen_test | 481 | 0.778 | 0.738 | 0.812 |
| Actinomycetota | all_test | unseen_test | 1193 | 0.697 | 0.671 | 0.723 |
| Bacillota | <30% | seen_test | 51 | 0.804 | 0.675 | 0.890 |
| Bacillota | <30% | unseen_test | 101 | 0.723 | 0.629 | 0.801 |
| Bacillota | 30-40% | seen_test | 139 | 0.727 | 0.647 | 0.794 |
| Bacillota | 30-40% | unseen_test | 237 | 0.814 | 0.760 | 0.859 |
| Bacillota | 40-50% | seen_test | 142 | 0.915 | 0.858 | 0.951 |
| Bacillota | 40-50% | unseen_test | 214 | 0.921 | 0.876 | 0.950 |
| Bacillota | 50-60% | seen_test | 25 | 1.000 | 0.867 | 1.000 |
| Bacillota | 50-60% | unseen_test | 119 | 0.958 | 0.905 | 0.982 |
| Bacillota | 60-70% | seen_test | 2 | 1.000 | 0.342 | 1.000 |
| Bacillota | 60-70% | unseen_test | 62 | 0.935 | 0.846 | 0.975 |
| Bacillota | 70-100% | seen_test | 5 | 1.000 | 0.566 | 1.000 |
| Bacillota | 70-100% | unseen_test | 21 | 1.000 | 0.845 | 1.000 |
| Bacillota | no_hit | seen_test | 127 | 0.583 | 0.496 | 0.665 |
| Bacillota | no_hit | unseen_test | 328 | 0.518 | 0.464 | 0.572 |
| Bacillota | all_test | seen_test | 491 | 0.770 | 0.731 | 0.805 |
| Bacillota | all_test | unseen_test | 1082 | 0.763 | 0.737 | 0.788 |
| Cyanobacteriota | <30% | seen_test | 48 | 0.771 | 0.635 | 0.867 |
| Cyanobacteriota | <30% | unseen_test | 13 | 0.462 | 0.232 | 0.709 |
| Cyanobacteriota | 30-40% | seen_test | 122 | 0.861 | 0.788 | 0.911 |
| Cyanobacteriota | 30-40% | unseen_test | 21 | 0.810 | 0.600 | 0.923 |
| Cyanobacteriota | 40-50% | seen_test | 161 | 0.876 | 0.816 | 0.918 |
| Cyanobacteriota | 40-50% | unseen_test | 22 | 0.909 | 0.722 | 0.975 |
| Cyanobacteriota | 50-60% | seen_test | 39 | 0.923 | 0.797 | 0.973 |
| Cyanobacteriota | 50-60% | unseen_test | 25 | 1.000 | 0.867 | 1.000 |
| Cyanobacteriota | 60-70% | seen_test | 12 | 1.000 | 0.758 | 1.000 |
| Cyanobacteriota | 60-70% | unseen_test | 4 | 1.000 | 0.510 | 1.000 |
| Cyanobacteriota | 70-100% | seen_test | 2 | 1.000 | 0.342 | 1.000 |
| Cyanobacteriota | 70-100% | unseen_test | 5 | 1.000 | 0.566 | 1.000 |
| Cyanobacteriota | no_hit | seen_test | 117 | 0.590 | 0.499 | 0.675 |
| Cyanobacteriota | no_hit | unseen_test | 56 | 0.500 | 0.373 | 0.627 |
| Cyanobacteriota | all_test | seen_test | 501 | 0.802 | 0.765 | 0.835 |
| Cyanobacteriota | all_test | unseen_test | 146 | 0.719 | 0.641 | 0.786 |

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

Across numerical-bin/group cells, nearest-homolog accuracy is often near ceiling at high identity and exceeds the logistic model in many hit-bearing cells (unweighted cell means: 0.934 vs 0.926). The logistic model improves on nearest-neighbor transfer in some low-identity cells and on every `no_hit` cell; baseline accuracy is zero for all 6,743 `no_hit` proteins because no donor prediction exists. The available-only baseline accuracy omits these rows and must not be compared to all-test logistic accuracy without noting its reduced coverage.

## 14. Interpretation

Overall seen/unseen accuracy gaps are -0.036 for Actinomycetota, +0.001 for Bacillota, and -0.002 for Cyanobacteriota. Within strata, 11 of 21 accuracy gaps are negative; 6 intervals lie wholly below zero and 5 wholly above zero. The direction is not consistent across phyla or bins:

- **Actinomycetota:** 6/6 numerical-bin accuracy gaps are negative (range -0.058 to -0.004); the separate `no_hit` gap is +0.264.
- **Bacillota:** 2/6 numerical-bin accuracy gaps are negative (range -0.004 to +0.100); the separate `no_hit` gap is -0.182.
- **Cyanobacteriota:** 1/6 numerical-bin accuracy gaps are negative (range -0.003 to +0.274); the separate `no_hit` gap is -0.064.

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
