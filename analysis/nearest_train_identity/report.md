# Nearest-training identity analysis

All 84,762 expected experiment-specific test records are represented. Existing split hashes are unchanged.
Sensitive MMseqs2 search: both coverages >=80%, E-value <=0.001, sensitivity 7.5. Values are approximate maxima over returned qualifying hits.

Identity statistics below are percentages, conditional on having a qualifying hit. Missing identities are never assigned zero.

| Held-out phylum | Group | Expected | With hit | No hit | Min | Median | Mean | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Actinomycetota | seen_test | 13,666 | 13,008 | 658 | 20.7 | 49.0 | 49.6 | 99.7 |
| Actinomycetota | unseen_test | 11,062 | 9,476 | 1,586 | 18.8 | 47.2 | 47.5 | 100.0 |
| Bacillota | seen_test | 11,962 | 11,038 | 924 | 21.0 | 48.9 | 48.9 | 99.7 |
| Bacillota | unseen_test | 26,888 | 24,874 | 2,014 | 18.8 | 49.7 | 49.9 | 99.6 |
| Cyanobacteriota | seen_test | 13,908 | 13,159 | 749 | 21.8 | 48.2 | 48.7 | 100.0 |
| Cyanobacteriota | unseen_test | 7,276 | 6,464 | 812 | 20.1 | 49.1 | 50.0 | 86.9 |

## Observed overlap

| Held-out phylum | Shared observed range (%) | Seen in range | Unseen in range | Maximum ECDF difference |
|---|---:|---:|---:|---:|
| Actinomycetota | 20.7-99.7 | 100.00% | 99.97% | 0.1557 |
| Bacillota | 21.0-99.6 | 99.98% | 99.98% | 0.1682 |
| Cyanobacteriota | 21.8-86.9 | 98.71% | 99.97% | 0.1044 |

Observed ranges overlap broadly. These comparisons describe hit-bearing proteins; different missing-hit rates can bias distribution comparisons. No conclusions about classifier performance are available.

## Validation and caveats

- All expected accessions are present exactly once per group; all targets belong to the corresponding final training set. Independent raw-hit maxima checks passed for all six groups.
- The first returned hit had lower reported identity than another returned hit for 34,519 test records, confirming the need to reduce all returned hits.
- Search is heuristic, uses default low-complexity masking, and optimizes alignment score. It does not guarantee the unrestricted mathematical maximum.
- An absence of qualifying hits can reflect coverage, significance, masking, or search sensitivity. Missing identities remain missing.
- Identities above 50% are possible across the existing clusters. A 100% local match need not be a full-sequence duplicate.
- MMseqs2 maps U/O to X internally, verified with the installed binary. Affected query/winning-target pairs are flagged; source sequences are unchanged.
- MMseqs2 emitted conservative worst-case prefilter disk-space warnings. Actual disk usage stayed within available capacity; all searches completed.
- Reported fident is maximized directly; reconstructed nident is used only for consistency checking. No extra precision is inferred.
- Four reducer regression tests passed. Preflight checks used eight test proteins per group against full training sets, plus an alignment-backtrace check.
- No similarity strata were selected and no classifier was trained. The 2% histogram intervals are fixed display intervals only.

## Files

- Reproducible script: ../../code/nearest_train_identity.py; regression tests: ../../code/test_nearest_train_identity.py.
- summary.tsv; distribution_overlap.tsv; identity_distributions.png and .pdf.
- Per phylum: train/seen_test/unseen_test FASTAs, MMseqs2 databases, raw hit tables, search logs, temporary files, per-group nearest-identity tables and no-hit lists, and combined test_nearest_train_identity.tsv.
- settings.json; commands_*.json; validation.json; independent_validation.json; postprocessing_revision.json; README.md; this report.md.
- Preflight artifacts are kept separately under pilot/ and ../nearest_train_identity_preflight/.
