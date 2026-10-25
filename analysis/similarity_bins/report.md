# Final sequence-similarity strata

Decision made before any classifier or embedding step, using only saved identity, group membership and sample support. No prediction or performance outputs were read.

## Support rule fixed before selection

Require at least 100 hit-bearing seen-test AND 100 hit-bearing unseen-test proteins in each numerical bin, separately for ALL three experiments. This pragmatic floor gives a worst-case standard error of 5 percentage points for a simple independent proportion (roughly +/-9.8 points at 95% under a normal approximation). It does not guarantee precision for class-specific metrics, small differences between groups, or correlated proteins within sequence clusters. Future uncertainty assessment should account for these limits.

The threshold was stated before candidate counts were evaluated. A stricter 200-per-group rule fails three group cells in the selected scheme (see support_sensitivity.tsv), requiring further merging; 100 retains useful resolution without accepting very sparse groups. Counts of 100 are an aggregate floor, not class-level support.

## Candidates and decision

| scheme | boundaries_percent | numerical_bins | smallest_group_count | failing_experiment_bins | selected |
| --- | --- | --- | --- | --- | --- |
| ten_percentage_points_with_below20 | [0, 20, 30, 40, 50, 60, 70, 80, 90, 100] | 9 | 0 | 8 | False |
| twenty_percentage_points_with_below20 | [0, 20, 40, 60, 80, 100] | 5 | 0 | 5 | False |
| coarse_twenty_with_merged_tails | [0, 40, 60, 100] | 3 | 870 | 0 | False |
| central_tens_broad_tails | [0, 30, 40, 50, 60, 100] | 5 | 155 | 0 | False |
| selected_common_adjacent_merges | [0, 30, 40, 50, 60, 70, 100] | 6 | 132 | 0 | True |

The selected scheme is the finest supported adjacent coarsening of the fixed 10-percentage-point grid with a separate below-20 interval. All 256 adjacent coarsenings were evaluated with an explicit resolution-only objective; no criterion rewards making distributions or group sizes look similar. Coarser supported alternatives sacrifice central resolution without being needed for the stated support floor.

**Final labels:** `<30%`, `30-40%`, `40-50%`, `50-60%`, `60-70%`, `70-100%`, `no_hit`.

Numerical boundaries in fractions: [0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0]. Intervals include their lower bound and exclude their upper bound, except the last includes 1.000. Thus exactly 30% is in 30-40%, exactly 70% is in 70-100%, and exactly 100% is retained. The first bin is [0,30%); this includes qualifying low identities, not imputed zeros.

## Tail merges

- Below 20%: one qualifying unseen protein each in Actinomycetota and Bacillota, and no seen proteins; no such proteins in Cyanobacteriota. Merge with 20-30% into <30%; no record is discarded.
- The 90-100% tail has only 3, 12, and 0 unseen proteins, respectively. Merging 80-100% still gives only 22 Actinomycetota unseen and 32 Cyanobacteriota unseen proteins. Merge these tails with 70-80% across ALL experiments, giving 70-100% and at least 132 per group.
- Central 30-40%, 40-50%, 50-60%, and 60-70% intervals already meet support and remain separate. Final smallest group support is 132.

## Final counts

| held_out_phylum | similarity_bin | seen_test | unseen_test | total |
| --- | --- | --- | --- | --- |
| Actinomycetota | <30% | 178 | 461 | 639 |
| Actinomycetota | 30-40% | 1170 | 1513 | 2683 |
| Actinomycetota | 40-50% | 6199 | 3908 | 10107 |
| Actinomycetota | 50-60% | 4199 | 2481 | 6680 |
| Actinomycetota | 60-70% | 887 | 981 | 1868 |
| Actinomycetota | 70-100% | 375 | 132 | 507 |
| Actinomycetota | no_hit | 658 | 1586 | 2244 |
| Bacillota | <30% | 268 | 727 | 995 |
| Bacillota | 30-40% | 1129 | 3510 | 4639 |
| Bacillota | 40-50% | 5220 | 8392 | 13612 |
| Bacillota | 50-60% | 3543 | 8078 | 11621 |
| Bacillota | 60-70% | 503 | 3352 | 3855 |
| Bacillota | 70-100% | 375 | 815 | 1190 |
| Bacillota | no_hit | 924 | 2014 | 2938 |
| Cyanobacteriota | <30% | 794 | 155 | 949 |
| Cyanobacteriota | 30-40% | 1398 | 715 | 2113 |
| Cyanobacteriota | 40-50% | 5750 | 2563 | 8313 |
| Cyanobacteriota | 50-60% | 3736 | 1978 | 5714 |
| Cyanobacteriota | 60-70% | 874 | 833 | 1707 |
| Cyanobacteriota | 70-100% | 607 | 220 | 827 |
| Cyanobacteriota | no_hit | 749 | 812 | 1561 |

candidate_bin_counts.tsv gives every candidate/experiment/bin count, both group proportions (denominator: the entire test group including no_hit), totals, and support flags. final_bin_counts.tsv supplies final counts by experiment, group and stratum, with both whole-test-group and hit-bearing proportions. no_hit has no hit-bearing proportion and is excluded from numerical-bin support decisions.

## Validation and interpretation

All 84,762 experiment-specific test records remain represented: 78,019 in exactly one numerical stratum and 6,743 in no_hit. Uniqueness is checked within experiment, since the same accession can legitimately appear in different held-out experiments. Combined and per-group identity inputs agree, test accessions/groups match frozen split files, assignment totals match each group, every original source column and identity string survives an output round trip unchanged, and all input-file hashes remain unchanged.

No-hit identities remain missing, never zero. These proteins must be reported and evaluated separately; no qualifying hit does not imply absence of homology. The upstream identities are approximate maxima among returned MMseqs2 local alignments requiring >=80% coverage of both sequences and E<=0.001, not guaranteed global maxima. Broad bins only coarsely control similarity: within-bin identity distributions and EC3 composition can still differ, and a count floor alone does not establish exchangeability or independent observations. In particular, 70-100% pools a broad high-identity tail. Do not infer comparable identities from equal bin labels alone.

No classifier was inspected or trained. The upstream 2% histogram intervals were not used as analysis strata. Existing nearest-identity and split files were not modified. No Git commit or push is performed by this script.

## Reproduction and files

Run `python code/make_similarity_bins.py` from the project root (or use a new `--output analysis/<name>`). The script refuses to overwrite an existing output directory. Requires pandas, numpy and matplotlib; no MMseqs2 rerun is needed.

- all_test_similarity_bins.tsv and <phylum>_similarity_bins.tsv: original per-protein columns plus similarity_bin.
- candidate_schemes.tsv, candidate_bin_counts.tsv, adjacent_coarsening_evaluation.tsv: alternatives and support decisions.
- final_bin_counts.tsv, final_bin_counts_wide.tsv, below20_counts.tsv, support_sensitivity.tsv: final counts and diagnostics.
- similarity_bin_counts.png and .pdf: grouped counts with no_hit shown separately by hatching.
- settings.json, validation.json, report.md: decision rule, provenance and validation.
