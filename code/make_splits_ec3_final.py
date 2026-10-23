from pathlib import Path
import json
import random
import shutil
import subprocess

import pandas as pd


# ============================================================
# CONFIGURATION
# ============================================================

# This script is meant to live in:
#   Cleaned Data - EC Level 3/code/
#
# It will automatically read:
#   Cleaned Data - EC Level 3/data/proteins_ec3.tsv
#
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
OUTPUT_DIR = PROJECT_DIR / "splits_ec3"

DATA_TSV = DATA_DIR / "proteins_ec3.tsv"

# Edit this path if MMseqs2 is somewhere else on your computer.
MMSEQS_WSL = "/home/hsnab/mmseqs/bin/mmseqs"

HELD_OUT_PHYLA = [
    "Actinomycetota",
    "Bacillota",
    "Cyanobacteriota",
]

# MMseqs2 clustering settings
CLUSTER_MIN_IDENTITY = 0.50
CLUSTER_COVERAGE = 0.80
CLUSTER_COV_MODE = 0

# Split settings
SEED = 355
SEEN_TEST_FRACTION = 0.10

# "10% of what remains" after selecting the seen-test clusters
VALIDATION_FRACTION_OF_REMAINDER = 0.10

# Apply this AFTER the split, using TRAINING counts only.
MIN_TRAIN_PER_CLASS = 10

# If a held-out test protein has an exact sequence copy in training,
# remove the TRAINING copy for that run.
REMOVE_EXACT_TRAIN_UNSEEN_DUPLICATES = True


# ============================================================
# HELPERS
# ============================================================

def safe_name(text):
    return "".join(
        c if c.isalnum() or c in "-_" else "_"
        for c in str(text)
    )

def to_wsl_path(path):
    path = str(Path(path).resolve())
    drive = path[0].lower()
    rest = path[2:].replace("\\", "/")
    return f"/mnt/{drive}{rest}"

def write_fasta(df, path):
    """Write accession + sequence from a dataframe to FASTA."""
    with path.open("w", encoding="utf-8") as f:
        for accession, sequence in zip(df["accession"], df["sequence"]):
            f.write(f">{accession}\n")
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i + 80] + "\n")


def run_command(command):
    print("\nRunning:")
    print(" ".join(map(str, command)))
    subprocess.run(command, check=True)


# ============================================================
# LOAD CLEANED DATA
# ============================================================

def load_cleaned_data():
    if not DATA_TSV.exists():
        raise FileNotFoundError(
            f"Could not find:\n{DATA_TSV}\n\n"
            "Put this script inside the ZIP's 'code' folder."
        )

    df = pd.read_csv(DATA_TSV, sep="\t")

    required = {
        "accession",
        "phylum",
        "ec3",
        "lab_tested",
        "sequence",
    }

    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"The cleaned TSV is missing required columns: {sorted(missing)}"
        )

    if df["accession"].duplicated().any():
        n = int(df["accession"].duplicated().sum())
        raise ValueError(
            f"The cleaned TSV contains {n} duplicated accessions."
        )

    if df["sequence"].isna().any():
        raise ValueError("The cleaned TSV contains missing sequences.")

    print("Loaded cleaned dataset")
    print(f"  Proteins:    {len(df):,}")
    print(f"  EC3 classes: {df['ec3'].nunique():,}")
    print(f"  Phyla:       {df['phylum'].nunique():,}")

    print("\nHeld-out phylum sizes in the CURRENT TSV:")
    for phylum in HELD_OUT_PHYLA:
        n = int(df["phylum"].eq(phylum).sum())
        print(f"  {phylum}: {n:,}")

    return df


# ============================================================
# MMSEQS2 CLUSTERING
# ============================================================

def cluster_nonheldout(df, held_out_phylum, run_dir):
    """
    IMPORTANT DESIGN CHOICE:

    The held-out phylum is removed FIRST.
    Only the remaining (seen) phyla are clustered.

    This prevents the unseen-test proteins from influencing how
    train / validation / seen_test are constructed.
    """
    nonheld = df.loc[
        ~df["phylum"].eq(held_out_phylum)
    ].copy()

    fasta_path = run_dir / "nonheldout.fasta"
    write_fasta(nonheld, fasta_path)

    cluster_prefix = run_dir / "clu50"
    cluster_tsv = run_dir / "clu50_cluster.tsv"
    tmp_dir = run_dir / "mmseqs_tmp"

    if cluster_tsv.exists():
        print(
            f"\n[{held_out_phylum}] Found existing "
            f"{cluster_tsv.name}; skipping MMseqs2 clustering."
        )
    else:
        command = [
            "wsl",
            MMSEQS_WSL,
            "easy-cluster",
            to_wsl_path(fasta_path),
            to_wsl_path(cluster_prefix),
            to_wsl_path(tmp_dir),
            "--min-seq-id", str(CLUSTER_MIN_IDENTITY),
            "-c", str(CLUSTER_COVERAGE),
            "--cov-mode", str(CLUSTER_COV_MODE),
        ]

        run_command(command)

    if not cluster_tsv.exists():
        raise RuntimeError(
            f"MMseqs2 did not create expected file:\n{cluster_tsv}"
        )

    clusters = pd.read_csv(
        cluster_tsv,
        sep="\t",
        header=None,
        names=["cluster50", "accession"],
    )

    if clusters["accession"].duplicated().any():
        raise ValueError(
            f"[{held_out_phylum}] An accession appears in more than one "
            "MMseqs2 cluster."
        )

    missing = set(nonheld["accession"]) - set(clusters["accession"])
    if missing:
        raise ValueError(
            f"[{held_out_phylum}] {len(missing)} non-held-out proteins "
            "are missing from the MMseqs2 cluster file."
        )

    nonheld = nonheld.merge(
        clusters,
        on="accession",
        how="left",
        validate="one_to_one",
    )

    return nonheld

# ============================================================
# WHOLE-CLUSTER SPLITTING
# ============================================================

def choose_clusters_close_to_target(cluster_sizes, target_n, seed):
    """
    Select whole clusters until their total number of proteins is
    approximately target_n.

    Cluster IDs are sorted before shuffling so the result is reproducible.
    """
    items = sorted(
        cluster_sizes.items(),
        key=lambda x: str(x[0]),
    )

    rng = random.Random(seed)
    rng.shuffle(items)

    chosen = []
    total = 0

    for cluster_id, size in items:
        if total >= target_n:
            break

        chosen.append(cluster_id)
        total += int(size)

    return set(chosen), total


def split_seen_phyla(nonheld):
    """
    Split ONLY the non-held-out phyla.

    About 10% of proteins -> seen_test
    About 10% of what remains -> validation
    Rest -> train

    Whole MMseqs2 clusters stay together.
    """
    work = nonheld.copy()
    work["group"] = "unassigned"

    cluster_sizes = (
        work.groupby("cluster50")
        .size()
        .to_dict()
    )

    n_total = len(work)

    target_seen = round(
        SEEN_TEST_FRACTION * n_total
    )

    seen_clusters, n_seen = choose_clusters_close_to_target(
        cluster_sizes,
        target_seen,
        SEED,
    )

    work.loc[
        work["cluster50"].isin(seen_clusters),
        "group",
    ] = "seen_test"

    remaining = work.loc[
        work["group"].eq("unassigned")
    ].copy()

    remaining_cluster_sizes = (
        remaining.groupby("cluster50")
        .size()
        .to_dict()
    )

    target_validation = round(
        VALIDATION_FRACTION_OF_REMAINDER * len(remaining)
    )

    validation_clusters, n_validation = choose_clusters_close_to_target(
        remaining_cluster_sizes,
        target_validation,
        SEED + 1,
    )

    work.loc[
        work["cluster50"].isin(validation_clusters)
        & work["group"].eq("unassigned"),
        "group",
    ] = "validation"

    work.loc[
        work["group"].eq("unassigned"),
        "group",
    ] = "train"

    print("\nInitial seen-phylum split:")
    print(
        work["group"]
        .value_counts()
        .reindex(["train", "validation", "seen_test"])
        .fillna(0)
        .astype(int)
        .to_string()
    )

    return work


# ============================================================
# ADD UNSEEN TEST
# ============================================================

def add_unseen_test(df, held_out_phylum, seen_split):
    unseen = df.loc[
        df["phylum"].eq(held_out_phylum)
    ].copy()

    # The unseen phylum was not part of clustering.
    unseen["cluster50"] = pd.NA
    unseen["group"] = "unseen_test"

    combined = pd.concat(
        [seen_split, unseen],
        ignore_index=True,
        sort=False,
    )

    return combined


# ============================================================
# EXACT TRAIN / UNSEEN DUPLICATES
# ============================================================

def remove_exact_train_unseen_duplicates(work):
    """
    If the exact amino-acid sequence of an unseen-test protein occurs
    in TRAINING, remove the training copy.

    This is generic and does not hard-code the two known sequence pairs.
    """
    if not REMOVE_EXACT_TRAIN_UNSEEN_DUPLICATES:
        return work.copy(), pd.DataFrame()

    unseen_sequences = set(
        work.loc[
            work["group"].eq("unseen_test"),
            "sequence",
        ]
    )

    remove_mask = (
        work["group"].eq("train")
        & work["sequence"].isin(unseen_sequences)
    )

    removed = work.loc[
        remove_mask,
        [
            "accession",
            "phylum",
            "ec3",
            "sequence",
            "cluster50",
            "group",
        ],
    ].copy()

    cleaned = work.loc[~remove_mask].copy()

    return cleaned, removed


# ============================================================
# TRAIN-ONLY EC3 CLASS FILTER
# ============================================================

def filter_classes_from_training(work):
    """
    Determine allowed EC3 classes FROM TRAINING ONLY.

    A class is retained only if training contains at least
    MIN_TRAIN_PER_CLASS proteins from that class.

    Then the same allowed class set is applied to:
      train / validation / seen_test / unseen_test

    Therefore every retained test class is actually learnable.
    """
    train_counts = (
        work.loc[
            work["group"].eq("train"),
            "ec3",
        ]
        .value_counts()
    )

    allowed_classes = set(
        train_counts.loc[
            train_counts >= MIN_TRAIN_PER_CLASS
        ].index
    )

    remove_mask = ~work["ec3"].isin(allowed_classes)

    removed = work.loc[remove_mask].copy()
    kept = work.loc[~remove_mask].copy()

    return kept, removed, train_counts, allowed_classes


# ============================================================
# VALIDATION CHECKS
# ============================================================

def lab_tested_percent(series):
    if len(series) == 0:
        return float("nan")

    return 100.0 * (
        series.astype(str)
        .str.strip()
        .str.lower()
        .eq("yes")
        .mean()
    )


def validate_split(work, held_out_phylum):
    """
    Check the split before any model is trained.
    """
    order = [
        "train",
        "validation",
        "seen_test",
        "unseen_test",
    ]

    # --------------------------------------------------------
    # 1. Held-out phylum must occur only in unseen_test.
    # --------------------------------------------------------
    unseen = work.loc[
        work["group"].eq("unseen_test")
    ]

    if not unseen["phylum"].eq(held_out_phylum).all():
        raise AssertionError(
            "unseen_test contains proteins outside the held-out phylum."
        )

    elsewhere = work.loc[
        ~work["group"].eq("unseen_test")
    ]

    if elsewhere["phylum"].eq(held_out_phylum).any():
        raise AssertionError(
            "Held-out phylum appears outside unseen_test."
        )

    # --------------------------------------------------------
    # 2. No non-held-out cluster may cross train/validation/seen_test.
    # --------------------------------------------------------
    seen_groups = work.loc[
        work["group"].isin(
            ["train", "validation", "seen_test"]
        )
    ]

    cluster_group_counts = (
        seen_groups.groupby("cluster50")["group"]
        .nunique()
    )

    bad_clusters = cluster_group_counts.loc[
        cluster_group_counts > 1
    ]

    if len(bad_clusters):
        raise AssertionError(
            f"{len(bad_clusters)} MMseqs2 clusters cross "
            "train/validation/seen_test."
        )

    # --------------------------------------------------------
    # 3. No exact sequence may occur in both training and either test.
    # --------------------------------------------------------
    train_sequences = set(
        work.loc[
            work["group"].eq("train"),
            "sequence",
        ]
    )

    seen_test_sequences = set(
        work.loc[
            work["group"].eq("seen_test"),
            "sequence",
        ]
    )

    unseen_test_sequences = set(
        work.loc[
            work["group"].eq("unseen_test"),
            "sequence",
        ]
    )

    train_seen_overlap = train_sequences & seen_test_sequences
    train_unseen_overlap = train_sequences & unseen_test_sequences

    if train_seen_overlap:
        raise AssertionError(
            f"{len(train_seen_overlap)} exact sequences occur in both "
            "train and seen_test."
        )

    if train_unseen_overlap:
        raise AssertionError(
            f"{len(train_unseen_overlap)} exact sequences occur in both "
            "train and unseen_test."
        )

    # --------------------------------------------------------
    # 4. Every retained validation/test class exists in training.
    # --------------------------------------------------------
    train_classes = set(
        work.loc[
            work["group"].eq("train"),
            "ec3",
        ]
    )

    for group in [
        "validation",
        "seen_test",
        "unseen_test",
    ]:
        group_classes = set(
            work.loc[
                work["group"].eq(group),
                "ec3",
            ]
        )

        missing = group_classes - train_classes

        if missing:
            raise AssertionError(
                f"{group} contains EC3 classes absent from training: "
                f"{sorted(missing)}"
            )

    # --------------------------------------------------------
    # 5. All retained training classes meet the threshold.
    # --------------------------------------------------------
    final_train_counts = (
        work.loc[
            work["group"].eq("train"),
            "ec3",
        ]
        .value_counts()
    )

    if not (
        final_train_counts >= MIN_TRAIN_PER_CLASS
    ).all():
        raise AssertionError(
            "At least one retained EC3 class has fewer than "
            f"{MIN_TRAIN_PER_CLASS} training proteins."
        )

    # --------------------------------------------------------
    # Summary table
    # --------------------------------------------------------
    rows = []

    for group in order:
        g = work.loc[
            work["group"].eq(group)
        ]

        rows.append({
            "group": group,
            "proteins": len(g),
            "ec3_classes": g["ec3"].nunique(),
            "phyla": g["phylum"].nunique(),
            "lab_tested_percent": lab_tested_percent(
                g["lab_tested"]
            ),
        })

    summary = pd.DataFrame(rows)

    print("\nFinal split summary:")
    print(summary.to_string(index=False))

    print("\nAll validation checks passed.")

    return summary


# ============================================================
# SAVE OUTPUT
# ============================================================

def save_run_outputs(
    work,
    removed_exact,
    removed_classes,
    train_counts,
    allowed_classes,
    summary,
    held_out_phylum,
    run_dir,
):
    # Compact reproducibility file requested in the plan.
    compact = work[
        [
            "accession",
            "cluster50",
            "group",
        ]
    ].copy()

    compact.to_csv(
        run_dir / f"split_{held_out_phylum}.tsv",
        sep="\t",
        index=False,
    )

    # Full convenience file.
    work.to_csv(
        run_dir / f"proteins_split_{held_out_phylum}.tsv",
        sep="\t",
        index=False,
    )

    removed_exact.to_csv(
        run_dir / "removed_exact_train_unseen_duplicates.tsv",
        sep="\t",
        index=False,
    )

    removed_classes.to_csv(
        run_dir / "removed_low_training_class.tsv",
        sep="\t",
        index=False,
    )

    train_counts.rename("training_count").to_csv(
        run_dir / "training_ec3_counts_before_filter.tsv",
        sep="\t",
        header=True,
    )

    pd.DataFrame({
        "ec3": sorted(allowed_classes)
    }).to_csv(
        run_dir / "allowed_ec3_classes.tsv",
        sep="\t",
        index=False,
    )

    summary.to_csv(
        run_dir / "split_summary.tsv",
        sep="\t",
        index=False,
    )

    settings = {
        "held_out_phylum": held_out_phylum,
        "seed": SEED,
        "seen_test_fraction": SEEN_TEST_FRACTION,
        "validation_fraction_of_remainder":
            VALIDATION_FRACTION_OF_REMAINDER,
        "min_train_per_class": MIN_TRAIN_PER_CLASS,
        "cluster_min_identity": CLUSTER_MIN_IDENTITY,
        "cluster_coverage": CLUSTER_COVERAGE,
        "cluster_cov_mode": CLUSTER_COV_MODE,
        "remove_exact_train_unseen_duplicates":
            REMOVE_EXACT_TRAIN_UNSEEN_DUPLICATES,
        "final_proteins": int(len(work)),
        "final_train_classes": int(
            work.loc[
                work["group"].eq("train"),
                "ec3",
            ].nunique()
        ),
    }

    with (
        run_dir / "settings.json"
    ).open("w", encoding="utf-8") as f:
        json.dump(
            settings,
            f,
            indent=2,
        )


# ============================================================
# ONE HELD-OUT-PHYLUM EXPERIMENT
# ============================================================

def build_one_split(df, held_out_phylum):
    print("\n" + "=" * 80)
    print(f"HELD-OUT PHYLUM: {held_out_phylum}")
    print("=" * 80)

    run_dir = OUTPUT_DIR / safe_name(held_out_phylum)
    run_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    # 1. Remove held-out phylum FIRST, then cluster only seen phyla.
    nonheld = cluster_nonheldout(
        df,
        held_out_phylum,
        run_dir,
    )

    # 2. Whole-cluster train / validation / seen-test split.
    seen_split = split_seen_phyla(nonheld)

    # 3. Add held-out phylum as unseen_test.
    work = add_unseen_test(
        df,
        held_out_phylum,
        seen_split,
    )

    print("\nBefore exact-duplicate and class filtering:")
    print(
        work["group"]
        .value_counts()
        .reindex(
            [
                "train",
                "validation",
                "seen_test",
                "unseen_test",
            ]
        )
        .fillna(0)
        .astype(int)
        .to_string()
    )

    # 4. Remove exact TRAINING copies of unseen-test sequences.
    work, removed_exact = (
        remove_exact_train_unseen_duplicates(work)
    )

    print(
        "\nExact training copies of unseen-test sequences removed: "
        f"{len(removed_exact)}"
    )

    if len(removed_exact):
        print(
            removed_exact[
                [
                    "accession",
                    "phylum",
                    "ec3",
                ]
            ].to_string(index=False)
        )

    # 5. Apply 10-per-class minimum using TRAINING ONLY.
    (
        work,
        removed_classes,
        train_counts,
        allowed_classes,
    ) = filter_classes_from_training(work)

    print(
        f"\nEC3 classes with >= {MIN_TRAIN_PER_CLASS} "
        f"training proteins: {len(allowed_classes)}"
    )

    print(
        "Proteins removed because their EC3 class did not meet "
        f"the training threshold: {len(removed_classes):,}"
    )

    # 6. Validate.
    summary = validate_split(
        work,
        held_out_phylum,
    )

    # 7. Save.
    save_run_outputs(
        work=work,
        removed_exact=removed_exact,
        removed_classes=removed_classes,
        train_counts=train_counts,
        allowed_classes=allowed_classes,
        summary=summary,
        held_out_phylum=held_out_phylum,
        run_dir=run_dir,
    )

    print(f"\nSaved outputs to:\n{run_dir}")

    return summary


# ============================================================
# MAIN
# ============================================================

def main():
    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    df = load_cleaned_data()

    missing_heldouts = [
        phylum
        for phylum in HELD_OUT_PHYLA
        if phylum not in set(df["phylum"])
    ]

    if missing_heldouts:
        raise ValueError(
            "These requested held-out phyla are not in the cleaned TSV: "
            f"{missing_heldouts}"
        )

    combined_summaries = []

    for held_out_phylum in HELD_OUT_PHYLA:
        summary = build_one_split(
            df,
            held_out_phylum,
        )

        summary.insert(
            0,
            "held_out_phylum",
            held_out_phylum,
        )

        combined_summaries.append(summary)

    combined = pd.concat(
        combined_summaries,
        ignore_index=True,
    )

    combined.to_csv(
        OUTPUT_DIR / "all_split_summaries.tsv",
        sep="\t",
        index=False,
    )

    print("\n" + "=" * 80)
    print("ALL THREE SPLITS COMPLETED")
    print("=" * 80)

    print(
        f"\nCombined summary saved to:\n"
        f"{OUTPUT_DIR / 'all_split_summaries.tsv'}"
    )

    print(
        "\nNEXT STEP AFTER THIS SCRIPT:\n"
        "1. Generate frozen ESM-2 embeddings.\n"
        "2. Train the multinomial logistic regression on TRAIN only.\n"
        "3. Tune settings using VALIDATION only.\n"
        "4. Predict SEEN_TEST and UNSEEN_TEST.\n"
        "5. Use MMseqs2 easy-search against TRAIN to compute S(x), "
        "the nearest-training sequence identity, for every test protein.\n"
        "6. Compare seen vs unseen performance only inside identity "
        "ranges with enough proteins from both groups."
    )


if __name__ == "__main__":
    main()
