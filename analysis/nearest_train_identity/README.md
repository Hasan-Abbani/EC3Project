# Nearest-training sequence identity

Reproduce from the project root:

```powershell
python code/nearest_train_identity.py
```

Dependencies: Python with pandas, numpy, matplotlib; WSL with
`/home/hsnab/mmseqs/bin/mmseqs`. To validate inputs without writing outputs:
`python code/nearest_train_identity.py --validate-only`.
`--pilot` runs eight test accessions per group against each complete training set
in a separate `pilot/` subdirectory. Pilot results are not full-dataset results.

## Inputs and frozen assignments

Each experiment reads its existing `proteins_split_<phylum>.tsv`, including its
final filtered `train` group and the original `seen_test` and `unseen_test`
groups. Compact `split_<phylum>.tsv` assignments must agree exactly. Sequences
are checked against `data/proteins_ec3.tsv`. No clustering, splitting, filtering
of retained proteins, embeddings, or model fitting is performed.
The pre-clustering `nonheldout.fasta` files are not final training FASTAs and
are deliberately not reused. All split TSVs and settings are hashed before
analysis and rechecked afterward.

## Operational definition and limits

Following the user's methodological choices, this is a **sensitive approximate
search**, requiring **at least 80% coverage of both query and target**. Results
also require E-value <= 0.001, the standard search significance cutoff used here.
Identity is `fident`, the reported fraction of identical alignment columns,
including gaps in alignment length, rather than percent units or bit score.
`max_train_identity` is the largest reported `fident` among returned qualifying
alignments to that experiment's final training set. Ties prefer higher bit
score, then lexicographically smaller training accession.

Search uses `-s 7.5 --alignment-mode 3 -a 1 --seq-id-mode 0 --min-seq-id 0
-c 0.8 --cov-mode 0 -e 0.001`. `--max-seqs`, `--max-accept`, and
`--max-rejected` are set to the experiment's entire training count, avoiding a
small top-hit cap. The k-mer prefilter and default low-complexity masking still
make the search heuristic. Alignment scoring uses MMseqs2 defaults (BLOSUM62,
gap-open 11, gap-extension 1, composition bias correction); full effective
parameters are in search logs. Alignments optimize score, not percent identity
over all possible alignments. This is not a proof of the mathematical maximum
over every training sequence or every possible alignment. First returned hits
are never presumed to maximize identity.

Whole-cluster splitting at a 50% threshold does not enforce a 50% upper bound
on every test/training pair. Values above 50% therefore remain in the analysis;
they do not trigger reassignment. Likewise 100% local identity need not imply
identical full sequences. For example, Actinomycetota unseen accession P0DUB9
aligns to training P0DUB8 at 100% identity over 319 columns, with query coverage
1.000 and target coverage 0.891. Exact full-sequence duplicate checks still pass.

MMseqs2 rounds serialized identity. The script maximizes the reported `fident`
directly and uses `nident/alnlen` only as a consistency check, not to manufacture
extra precision. Source references inspected:
[identity computation](https://github.com/soedinglab/MMseqs2/blob/master/src/alignment/Matcher.cpp)
and [conversion fields](https://github.com/soedinglab/MMseqs2/blob/master/src/util/convertalignments.cpp).
The installed version hash is saved in `settings.json`.

The cleaned dataset contains 52 proteins with rare amino acids U or O, which
the standard MMseqs2 amino-acid alphabet maps to X internally. These sequences
are retained as supplied. Per-protein outputs flag U/O in both query and winning
target, and summaries count affected queries and winning targets. Consequently,
reported MMseqs2 identity can differ from literal identity of original letters
at these positions. See the upstream
[amino-acid mapping](https://github.com/soedinglab/MMseqs2/blob/master/src/commons/SubstitutionMatrix.cpp).
An installed-binary probe confirmed this: a 351-residue query containing U
versus a synthetic copy with U replaced by O returned `fident=1.000` over all
351 columns. Probe FASTAs and output are stored under
`../nearest_train_identity_preflight/rare_residue_probe/` and are not study
proteins or split changes.

## Missing hits and distributions

Every expected test accession appears once in its result table. A protein with
no qualifying hit has an empty `max_train_identity` and
`hit_status=no_qualifying_hit`; it is never assigned zero. This may reflect
coverage, significance, masking, or heuristic search failure and does not
establish absence of homology. No extra fallback method mixes definitions.

Summary statistics, histograms, ECDFs, and overlap measures describe proteins
**with qualifying hits only**. Missing counts are reported separately, so
distribution comparisons must account for possible selective missingness.
`distribution_overlap.tsv` reports the intersection of observed ranges, the
fraction of each hit-bearing group in that interval, and the maximum absolute
ECDF difference. Range overlap does not imply equal distributions or establish
adequate support for a future controlled performance comparison.

Histogram intervals are fixed at 2 percentage points for display only. They
are not selected analysis strata or similarity bins for model comparisons.

## Outputs

- `summary.tsv`: expected, hit-bearing, missing counts and min/median/mean/max.
- `distribution_overlap.tsv`: descriptive overlap and ECDF distance.
- `identity_distributions.png` and `.pdf`: all three experiments, seen vs unseen.
- `settings.json`, `commands_*.json`, `validation.json`: configuration, executed
  commands, input hashes, target/accession checks, and split preservation.
- Each phylum directory contains final `train.fasta`, `seen_test.fasta`,
  `unseen_test.fasta`; MMseqs2 databases, alignments, temporary files, logs, raw
  `<group>_hits.tsv` (headerless, schema in settings), search completion markers;
  `<group>_nearest_train_identity.tsv`, `<group>_no_hits.tsv`, and combined
  `test_nearest_train_identity.tsv`. Per-protein tables preserve original
  accession, phylum, EC3 class, cluster, and group without duplicating sequences.

Reruns verify FASTA contents and hit checksums before reuse. A changed input,
configuration, or script hash requires a new output directory. Intermediate
pilot files are separate from full analysis outputs.
