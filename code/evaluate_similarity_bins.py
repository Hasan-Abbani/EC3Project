"""Evaluate frozen logistic-regression predictions within frozen identity bins.

Run from any directory with:
    python code/evaluate_similarity_bins.py

This is a post-processing step. It does not retrain models, regenerate
embeddings, run MMseqs2, or modify any upstream files.
"""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
PHYLA = ["Actinomycetota", "Bacillota", "Cyanobacteriota"]
GROUPS = ["seen_test", "unseen_test"]
BINS = ["<30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-100%", "no_hit"]
ALL_BINS = [*BINS, "all_test"]
DEFAULT_OUTPUT = ROOT / "analysis" / "binwise_evaluation"
BOOTSTRAP_REPLICATES = 2000
BOOTSTRAP_SEED = 355
CI_LEVEL = 0.95
PREDICTION_COLUMNS = {
    "accession", "group", "lab_tested", "bin", "true_ec3", "pred_ec3",
    "correct", "nearest_train_accession", "baseline1_pred", "baseline1_correct",
}


def read_tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _encoded_metrics(true_labels, predicted_labels, macro_labels=None):
    """Return accuracy, macro-F1 and balanced accuracy from string labels."""
    truth = np.asarray(true_labels, dtype=str)
    predicted = np.asarray(predicted_labels, dtype=str)
    if len(truth) == 0:
        return {"accuracy": np.nan, "macro_f1": np.nan, "balanced_accuracy": np.nan}

    labels = sorted(set(truth.tolist()) | set(predicted.tolist()))
    label_index = {label: index for index, label in enumerate(labels)}
    true_codes = np.fromiter((label_index[value] for value in truth), dtype=np.int64)
    predicted_codes = np.fromiter((label_index[value] for value in predicted), dtype=np.int64)
    n_labels = len(labels)
    confusion = np.bincount(
        true_codes * n_labels + predicted_codes,
        minlength=n_labels * n_labels,
    ).reshape(n_labels, n_labels)
    true_support = confusion.sum(axis=1)
    predicted_support = confusion.sum(axis=0)
    true_positive = np.diag(confusion)
    f1_denominator = true_support + predicted_support
    f1 = np.divide(
        2 * true_positive,
        f1_denominator,
        out=np.zeros(n_labels, dtype=float),
        where=f1_denominator != 0,
    )
    if macro_labels is None:
        included = np.arange(n_labels)
    else:
        included = np.asarray([label_index[label] for label in macro_labels if label in label_index])
    recalls = np.divide(
        true_positive,
        true_support,
        out=np.zeros(n_labels, dtype=float),
        where=true_support != 0,
    )
    present_truth = np.flatnonzero(true_support > 0)
    return {
        "accuracy": float(true_positive.sum() / len(truth)),
        "macro_f1": float(f1[included].mean()) if len(included) else np.nan,
        "balanced_accuracy": float(recalls[present_truth].mean()) if len(present_truth) else np.nan,
    }


def calculate_metrics(frame, shared_classes=None):
    true_labels = frame["true_ec3"].to_numpy(dtype=str)
    predicted_labels = frame["pred_ec3"].to_numpy(dtype=str)
    labels = None if shared_classes is None else sorted(shared_classes)
    metrics = _encoded_metrics(true_labels, predicted_labels, labels)
    metrics["n_macro_f1_labels"] = (
        len(set(true_labels.tolist()) | set(predicted_labels.tolist()))
        if labels is None else len(labels)
    )
    metrics["n_true_classes"] = len(set(true_labels.tolist()))
    return metrics


def bootstrap_metrics(frame, replicates, rng):
    """Ordinary protein-level bootstrap with a fixed label set per group."""
    truth = frame["true_ec3"].to_numpy(dtype=str)
    predicted = frame["pred_ec3"].to_numpy(dtype=str)
    n = len(frame)
    if n == 0:
        return {name: np.full(replicates, np.nan) for name in
                ("accuracy", "macro_f1", "balanced_accuracy")}
    labels = sorted(set(truth.tolist()) | set(predicted.tolist()))
    label_index = {label: index for index, label in enumerate(labels)}
    true_codes = np.fromiter((label_index[value] for value in truth), dtype=np.int64)
    predicted_codes = np.fromiter((label_index[value] for value in predicted), dtype=np.int64)
    n_labels = len(labels)
    scores = np.empty((replicates, 3), dtype=float)
    for replicate in range(replicates):
        sample = rng.integers(0, n, size=n)
        confusion = np.bincount(
            true_codes[sample] * n_labels + predicted_codes[sample],
            minlength=n_labels * n_labels,
        ).reshape(n_labels, n_labels)
        true_support = confusion.sum(axis=1)
        predicted_support = confusion.sum(axis=0)
        true_positive = np.diag(confusion)
        f1_denominator = true_support + predicted_support
        f1 = np.divide(
            2 * true_positive,
            f1_denominator,
            out=np.zeros(n_labels, dtype=float),
            where=f1_denominator != 0,
        )
        recalls = np.divide(
            true_positive,
            true_support,
            out=np.zeros(n_labels, dtype=float),
            where=true_support != 0,
        )
        present_truth = true_support > 0
        scores[replicate] = (
            true_positive.sum() / n,
            f1.mean(),
            recalls[present_truth].mean(),
        )
    return dict(zip(("accuracy", "macro_f1", "balanced_accuracy"), scores.T))


def percentile_interval(values):
    values = np.asarray(values, dtype=float)
    return tuple(float(value) for value in np.nanpercentile(values, [2.5, 97.5]))


def wilson_interval(successes, total):
    if total == 0:
        return np.nan, np.nan
    z = 1.959963984540054
    proportion = successes / total
    denominator = 1 + z * z / total
    center = (proportion + z * z / (2 * total)) / denominator
    half_width = z * np.sqrt(
        proportion * (1 - proportion) / total + z * z / (4 * total * total)
    ) / denominator
    return float(center - half_width), float(center + half_width)


def _group_subset(data, group, similarity_bin):
    if similarity_bin == "all_test":
        return data.loc[data["group"].eq(group)]
    return data.loc[data["group"].eq(group) & data["bin"].eq(similarity_bin)]


def _metric_row(phylum, similarity_bin, group, subset, metrics):
    return {
        "held_out_phylum": phylum,
        "similarity_bin": similarity_bin,
        "group": group,
        "n_proteins": len(subset),
        "n_ec3_classes": metrics["n_true_classes"],
        "n_macro_f1_labels": metrics["n_macro_f1_labels"],
        "accuracy": metrics["accuracy"],
        "macro_f1": metrics["macro_f1"],
        "balanced_accuracy": metrics["balanced_accuracy"],
    }


def _validate_experiment(phylum, predictions, bins_path, split_path, count_table):
    required = PREDICTION_COLUMNS - set(predictions.columns)
    if required:
        raise ValueError(f"{phylum}: prediction columns missing: {sorted(required)}")
    if predictions["accession"].eq("").any() or predictions["accession"].duplicated().any():
        raise ValueError(f"{phylum}: prediction accessions are empty or duplicated")
    if not set(predictions["group"]).issubset(GROUPS) or set(predictions["group"]) != set(GROUPS):
        raise ValueError(f"{phylum}: unexpected or missing seen/unseen group")
    if not set(predictions["bin"]).issubset(BINS):
        invalid = sorted(set(predictions["bin"]) - set(BINS))
        raise ValueError(f"{phylum}: invalid frozen similarity bin(s): {invalid}")
    if predictions[["true_ec3", "pred_ec3"]].eq("").any().any():
        raise ValueError(f"{phylum}: evaluated model rows have missing EC3 labels")
    if not predictions["lab_tested"].isin(["yes", "no"]).all():
        raise ValueError(f"{phylum}: lab_tested has values other than yes/no")
    if not predictions["correct"].isin(["0", "1"]).all():
        raise ValueError(f"{phylum}: correct has values other than 0/1")
    calculated_correct = predictions["true_ec3"].eq(predictions["pred_ec3"])
    if not calculated_correct.equals(predictions["correct"].eq("1")):
        raise ValueError(f"{phylum}: stored model correctness disagrees with labels")

    frozen = read_tsv(bins_path)
    needed_bin_columns = {"accession", "group", "similarity_bin"}
    if not needed_bin_columns.issubset(frozen.columns):
        raise ValueError(f"{phylum}: frozen bin file lacks required columns")
    if frozen["accession"].duplicated().any():
        raise ValueError(f"{phylum}: frozen bin assignments have duplicate accessions")
    if not frozen["group"].isin(GROUPS).all() or not frozen["similarity_bin"].isin(BINS).all():
        raise ValueError(f"{phylum}: frozen assignments contain an invalid group or bin")

    actual_indexed = predictions.set_index("accession")[["group", "bin"]].sort_index()
    frozen_indexed = (frozen.set_index("accession")[["group", "similarity_bin"]]
                      .rename(columns={"similarity_bin": "bin"}).sort_index())
    if not actual_indexed.index.equals(frozen_indexed.index):
        missing_predictions = sorted(set(frozen_indexed.index) - set(actual_indexed.index))
        extra_predictions = sorted(set(actual_indexed.index) - set(frozen_indexed.index))
        raise ValueError(
            f"{phylum}: prediction/frozen-bin accessions disagree; "
            f"missing predictions={len(missing_predictions)}, extra predictions={len(extra_predictions)}"
        )
    if not actual_indexed.equals(frozen_indexed):
        raise ValueError(f"{phylum}: prediction groups/bins disagree with frozen assignments")

    split = read_tsv(split_path)
    if not {"accession", "group"}.issubset(split.columns) or split["accession"].duplicated().any():
        raise ValueError(f"{phylum}: invalid frozen split table")
    split_test = split.loc[split["group"].isin(GROUPS)].set_index("accession")["group"].sort_index()
    if not split_test.index.equals(frozen_indexed.index):
        raise ValueError(f"{phylum}: frozen bin accessions do not match frozen test split")
    if not split_test.equals(frozen_indexed["group"]):
        raise ValueError(f"{phylum}: frozen bin groups do not match frozen test split")

    final_counts = count_table.loc[count_table["held_out_phylum"].eq(phylum)]
    count_rows = []
    for group in GROUPS:
        for similarity_bin in BINS:
            observed = int((predictions["group"].eq(group) &
                            predictions["bin"].eq(similarity_bin)).sum())
            expected_rows = final_counts.loc[
                final_counts["group"].eq(group) &
                final_counts["similarity_bin"].eq(similarity_bin), "count"
            ]
            if len(expected_rows) != 1 or observed != int(expected_rows.iloc[0]):
                expected = None if len(expected_rows) != 1 else int(expected_rows.iloc[0])
                raise ValueError(
                    f"{phylum} {group} {similarity_bin}: prediction count {observed} "
                    f"does not match frozen count {expected}"
                )
            count_rows.append({
                "group": group,
                "similarity_bin": similarity_bin,
                "prediction_count": observed,
                "frozen_count": int(expected_rows.iloc[0]),
                "reconciled": True,
            })

    baseline_pred_missing = predictions["baseline1_pred"].eq("")
    baseline_correct_missing = predictions["baseline1_correct"].eq("")
    if not baseline_pred_missing.equals(baseline_correct_missing):
        raise ValueError(f"{phylum}: baseline prediction/correctness missingness differs")
    nonmissing_correct = predictions.loc[~baseline_correct_missing, "baseline1_correct"]
    if not nonmissing_correct.isin(["0", "1"]).all():
        raise ValueError(f"{phylum}: baseline1_correct has values other than 0/1")
    if not predictions.loc[~baseline_pred_missing, "baseline1_pred"].eq(
            predictions.loc[~baseline_correct_missing, "true_ec3"]).equals(
            nonmissing_correct.eq("1")):
        raise ValueError(f"{phylum}: stored baseline correctness disagrees with labels")
    unexplained_missing = int((baseline_pred_missing & predictions["bin"].ne("no_hit")).sum())

    return {
        "prediction_rows": len(predictions),
        "unique_accessions": bool(not predictions["accession"].duplicated().any()),
        "valid_bins": True,
        "valid_groups": True,
        "nonmissing_model_labels": True,
        "stored_correctness_matches_labels": True,
        "frozen_assignments_match_predictions": True,
        "frozen_assignments_match_split": True,
        "bin_counts": count_rows,
        "baseline_missing_predictions": int(baseline_pred_missing.sum()),
        "baseline_missing_on_no_hit": int((baseline_pred_missing & predictions["bin"].eq("no_hit")).sum()),
        "baseline_missing_on_hit": unexplained_missing,
        "baseline_missing_reason": (
            "No qualifying training alignment, so the nearest-homolog label-transfer "
            "baseline has no donor accession."
            if unexplained_missing == 0 else
            f"{unexplained_missing} missing baseline prediction(s) occur outside no_hit; "
            "their cause is not established by the available files."
        ),
    }


def _class_distribution_tv(phylum, similarity_bin, seen, unseen):
    seen_counts = seen["true_ec3"].value_counts()
    unseen_counts = unseen["true_ec3"].value_counts()
    classes = sorted(set(seen_counts.index) | set(unseen_counts.index))
    seen_total = len(seen)
    unseen_total = len(unseen)
    rows = []
    tv = 0.5 * sum(
        abs((seen_counts.get(ec3, 0) / seen_total if seen_total else 0) -
            (unseen_counts.get(ec3, 0) / unseen_total if unseen_total else 0))
        for ec3 in classes
    )
    for ec3 in classes:
        seen_count = int(seen_counts.get(ec3, 0))
        unseen_count = int(unseen_counts.get(ec3, 0))
        seen_frequency = seen_count / seen_total if seen_total else np.nan
        unseen_frequency = unseen_count / unseen_total if unseen_total else np.nan
        rows.append({
            "held_out_phylum": phylum,
            "similarity_bin": similarity_bin,
            "ec3": ec3,
            "seen_count": seen_count,
            "seen_frequency": seen_frequency,
            "unseen_count": unseen_count,
            "unseen_frequency": unseen_frequency,
            "absolute_frequency_difference": abs(seen_frequency - unseen_frequency),
            "total_variation_distance": tv,
        })
    return rows, float(tv), len(set(seen_counts.index) & set(unseen_counts.index))


def _shared_metrics(seen, unseen, shared_classes):
    seen_common = seen.loc[seen["true_ec3"].isin(shared_classes)]
    unseen_common = unseen.loc[unseen["true_ec3"].isin(shared_classes)]
    seen_metrics = calculate_metrics(seen_common, shared_classes)
    unseen_metrics = calculate_metrics(unseen_common, shared_classes)
    row = {
        "n_shared_classes": len(shared_classes),
        "seen_n": len(seen_common),
        "unseen_n": len(unseen_common),
    }
    for metric in ("accuracy", "macro_f1", "balanced_accuracy"):
        row[f"seen_{metric}"] = seen_metrics[metric]
        row[f"unseen_{metric}"] = unseen_metrics[metric]
        row[f"gap_{metric}"] = unseen_metrics[metric] - seen_metrics[metric]
    row["seen_macro_f1_label_count"] = seen_metrics["n_macro_f1_labels"]
    row["unseen_macro_f1_label_count"] = unseen_metrics["n_macro_f1_labels"]
    return row


def _reweighted_accuracy(seen, unseen):
    seen_counts = seen["true_ec3"].value_counts()
    unseen_counts = unseen["true_ec3"].value_counts()
    shared = sorted(set(seen_counts.index) & set(unseen_counts.index))
    seen_common = seen.loc[seen["true_ec3"].isin(shared)]
    unseen_common = unseen.loc[unseen["true_ec3"].isin(shared)]
    base = {
        "seen_n": len(seen),
        "unseen_n": len(unseen),
        "common_support_class_count": len(shared),
        "seen_common_support_n": len(seen_common),
        "unseen_common_support_n": len(unseen_common),
        "seen_excluded_n": len(seen) - len(seen_common),
        "unseen_excluded_n": len(unseen) - len(unseen_common),
        "seen_only_class_count": len(set(seen_counts.index) - set(unseen_counts.index)),
        "unseen_only_class_count": len(set(unseen_counts.index) - set(seen_counts.index)),
    }
    if not shared:
        base.update({
            "class_mix_adjusted_seen_accuracy": np.nan,
            "unseen_common_support_accuracy": np.nan,
            "adjusted_gap_unseen_minus_seen": np.nan,
            "weight_sum": 0.0,
        })
        return base

    seen_common_correct = seen_common["true_ec3"].eq(seen_common["pred_ec3"])
    weights = np.array([
        (unseen_counts[ec3] / len(unseen_common)) / seen_counts[ec3]
        for ec3 in seen_common["true_ec3"]
    ], dtype=float)
    adjusted_seen = float(np.average(seen_common_correct.to_numpy(dtype=float), weights=weights))
    unseen_accuracy = float(unseen_common["true_ec3"].eq(unseen_common["pred_ec3"]).mean())
    base.update({
        "class_mix_adjusted_seen_accuracy": adjusted_seen,
        "unseen_common_support_accuracy": unseen_accuracy,
        "adjusted_gap_unseen_minus_seen": unseen_accuracy - adjusted_seen,
        "weight_sum": float(weights.sum()),
    })
    return base


def _baseline_metrics(subset):
    missing = subset["baseline1_pred"].eq("")
    available = subset.loc[~missing]
    correct = available["baseline1_pred"].eq(available["true_ec3"])
    n = len(subset)
    return {
        "n_proteins": n,
        "n_baseline_predictions": len(available),
        "n_missing_baseline_predictions": int(missing.sum()),
        "n_missing_on_no_hit": int((missing & subset["bin"].eq("no_hit")).sum()),
        "n_missing_on_hit": int((missing & subset["bin"].ne("no_hit")).sum()),
        "accuracy_all_test_missing_as_wrong": float(correct.sum() / n) if n else np.nan,
        "accuracy_when_baseline_available": float(correct.mean()) if len(available) else np.nan,
    }


def _md_table(frame, columns, formats=None):
    formats = formats or {}
    lines = ["| " + " | ".join(columns) + " |",
             "| " + " | ".join("---" for _ in columns) + " |"]
    for row in frame[columns].itertuples(index=False, name=None):
        values = []
        for column, value in zip(columns, row):
            if pd.isna(value):
                values.append("NA")
            elif column in formats:
                values.append(formats[column].format(float(value)))
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def _make_plots(metrics, intervals, output):
    colors = {"seen_test": "#0072B2", "unseen_test": "#D55E00"}
    positions = np.arange(len(BINS), dtype=float)
    figure, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True, layout="constrained")
    for axis, phylum in zip(axes, PHYLA):
        for offset, group in [(-0.13, "seen_test"), (0.13, "unseen_test")]:
            group_metrics = (metrics.loc[
                metrics["held_out_phylum"].eq(phylum) &
                metrics["group"].eq(group) &
                metrics["similarity_bin"].isin(BINS)
            ].set_index("similarity_bin").reindex(BINS))
            group_intervals = (intervals.loc[
                intervals["held_out_phylum"].eq(phylum) &
                intervals["group"].eq(group) &
                intervals["metric"].eq("accuracy") &
                intervals["similarity_bin"].isin(BINS)
            ].set_index("similarity_bin").reindex(BINS))
            estimates = group_metrics["accuracy"].to_numpy(dtype=float)
            lower = group_intervals["ci_low"].to_numpy(dtype=float)
            upper = group_intervals["ci_high"].to_numpy(dtype=float)
            axis.errorbar(
                positions + offset, estimates,
                yerr=np.vstack((estimates - lower, upper - estimates)),
                fmt="o", capsize=3, color=colors[group], label=group,
            )
        axis.axvline(5.5, color="gray", linestyle="--", linewidth=0.8)
        axis.set(title=phylum, ylabel="Accuracy", ylim=(0, 1.02))
        axis.grid(axis="y", alpha=0.2)
        axis.legend(loc="lower right", fontsize=8)
    axes[-1].set_xticks(positions, BINS, rotation=20, ha="right")
    figure.suptitle("Logistic-regression accuracy by frozen similarity bin\n"
                    "Points and 95% protein-bootstrap percentile intervals")
    figure.savefig(output / "accuracy_by_similarity_bin.png", dpi=180)
    figure.savefig(output / "accuracy_by_similarity_bin.pdf")
    plt.close(figure)

    figure, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True, layout="constrained")
    for axis, phylum in zip(axes, PHYLA):
        selected = (intervals.loc[
            intervals["held_out_phylum"].eq(phylum) &
            intervals["group"].eq("gap") &
            intervals["metric"].eq("accuracy") &
            intervals["similarity_bin"].isin(BINS)
        ].set_index("similarity_bin").reindex(BINS))
        estimates = selected["estimate"].to_numpy(dtype=float)
        lower = selected["ci_low"].to_numpy(dtype=float)
        upper = selected["ci_high"].to_numpy(dtype=float)
        axis.errorbar(
            positions, estimates,
            yerr=np.vstack((estimates - lower, upper - estimates)),
            fmt="o", capsize=3, color="#332288",
        )
        axis.axhline(0, color="black", linewidth=0.9)
        axis.axvline(5.5, color="gray", linestyle="--", linewidth=0.8)
        axis.set(title=phylum, ylabel="Unseen - seen accuracy gap")
        axis.grid(axis="y", alpha=0.2)
    axes[-1].set_xticks(positions, BINS, rotation=20, ha="right")
    figure.suptitle("Seen/unseen accuracy gap by frozen similarity bin\n"
                    "95% protein-bootstrap percentile intervals; zero = no gap")
    figure.savefig(output / "accuracy_gap_by_similarity_bin.png", dpi=180)
    figure.savefig(output / "accuracy_gap_by_similarity_bin.pdf")
    plt.close(figure)


def _build_report(binwise, gaps, tv_summary, shared, reweighted, intervals, lab, baseline,
                  validation, settings):
    accuracy_gaps = gaps.loc[gaps["similarity_bin"].isin(BINS)]
    numerical_gaps = accuracy_gaps.loc[accuracy_gaps["similarity_bin"].ne("no_hit")]
    n_negative = int(accuracy_gaps["gap_accuracy"].lt(0).sum())
    n_total = len(accuracy_gaps)
    n_ci_below = int((
        intervals.loc[
            intervals["group"].eq("gap") & intervals["metric"].eq("accuracy") &
            intervals["similarity_bin"].isin(BINS), "ci_high"
        ] < 0
    ).sum())
    n_ci_above = int((
        intervals.loc[
            intervals["group"].eq("gap") & intervals["metric"].eq("accuracy") &
            intervals["similarity_bin"].isin(BINS), "ci_low"
        ] > 0
    ).sum())
    overall = gaps.loc[gaps["similarity_bin"].eq("all_test")]
    lab_small = lab.loc[lab["n_proteins"].lt(100)]
    baseline_overview = baseline.loc[baseline["similarity_bin"].eq("all_test")]
    baseline_hit = baseline.loc[baseline["similarity_bin"].isin(BINS[:-1])]
    model_hit = binwise.loc[binwise["similarity_bin"].isin(BINS[:-1])]
    no_hit_baseline = baseline.loc[baseline["similarity_bin"].eq("no_hit")]
    numerical_patterns = []
    for phylum in PHYLA:
        rows = numerical_gaps.loc[numerical_gaps["held_out_phylum"].eq(phylum)]
        no_hit_gap = float(gaps.loc[
            gaps["held_out_phylum"].eq(phylum) & gaps["similarity_bin"].eq("no_hit"),
            "gap_accuracy",
        ].iloc[0])
        numerical_patterns.append(
            f"- **{phylum}:** {int(rows['gap_accuracy'].lt(0).sum())}/6 numerical-bin "
            f"accuracy gaps are negative (range {rows['gap_accuracy'].min():+.3f} to "
            f"{rows['gap_accuracy'].max():+.3f}); the separate `no_hit` gap is "
            f"{no_hit_gap:+.3f}."
        )

    lines = [
        "# Bin-wise logistic-regression evaluation",
        "",
        "## 1. Objective",
        "",
        "Assess whether frozen ESM-2/logistic-regression performance differs between "
        "`seen_test` and `unseen_test` within the already-defined nearest-training-identity "
        "strata. This is post-processing of existing predictions, not a causal test of taxonomy.",
        "",
        "## 2. Inputs and frozen upstream decisions",
        "",
        f"Experiments: {', '.join(PHYLA)}. Target: EC3. The prediction rows, split assignments, "
        "and similarity bins are read from their frozen upstream files. The numerical bins "
        "remain `<30%`, `30-40%`, `40-50%`, `50-60%`, `60-70%`, and `70-100%`; `no_hit` "
        "remains separate. No model was retrained, and no embedding, MMseqs2 search, split, "
        "filter, or bin definition was regenerated.",
        "",
        "## 3. Validation checks",
        "",
        "All validation checks passed. Prediction accessions are unique within each experiment; "
        "all rows have a valid group/bin and nonmissing true/predicted EC3; stored correctness "
        "matches the labels; per-accession bin and group assignments match the frozen bin and "
        "split files; and every group-by-bin count reconciles with `final_bin_counts.tsv`. "
        "See `validation.json` for counts, baseline missingness, and source hashes.",
        "",
        "## 4. Statistical methods",
        "",
        f"Accuracy, macro-F1, and balanced accuracy are computed from the existing test predictions. "
        "Macro-F1 averages over labels present in either the true or predicted labels for that "
        "group/bin; zero F1 is assigned to a label with no true or predicted support in a bootstrap "
        "replicate. Balanced accuracy averages per-class recall over classes represented in the "
        "true labels for that evaluation sample. Shared-class macro-F1 is limited to the classes "
        "shared by both groups. The sign convention is `gap = unseen - seen` for every metric.",
        "",
        f"Uncertainty uses {settings['bootstrap_replicates']:,} ordinary protein-level bootstrap "
        f"replicates with seed {settings['bootstrap_seed']}, with independent resampling within "
        "seen/unseen groups and percentile 95% intervals. Cluster bootstrap was not used: split "
        "metadata provides cluster50 IDs for seen-test proteins, but the held-out unseen phylum "
        "was intentionally not clustered, so both comparison groups cannot be resampled on a "
        "clean, common cluster mapping. Intervals do not correct for multiple comparisons.",
        "",
        "## 5. Bin-wise sample sizes",
        "",
        _md_table(
            binwise.loc[binwise["similarity_bin"].isin(BINS)],
            ["held_out_phylum", "similarity_bin", "group", "n_proteins", "n_ec3_classes"],
        ),
        "",
        "The original binning rule required at least 100 proteins per group in every numerical "
        "bin. All bin/group counts are checked against the frozen assignments. That aggregate "
        "floor does not ensure adequate per-class support; the smallest counts and shared-class "
        "supports should be considered when interpreting macro metrics.",
        "",
        "## 6. Raw seen-vs-unseen results",
        "",
        _md_table(
            gaps.loc[gaps["similarity_bin"].isin([*BINS, "all_test"])],
            ["held_out_phylum", "similarity_bin", "seen_n", "unseen_n",
             "seen_accuracy", "unseen_accuracy", "gap_accuracy",
             "seen_macro_f1", "unseen_macro_f1", "gap_macro_f1",
             "seen_balanced_accuracy", "unseen_balanced_accuracy", "gap_balanced_accuracy"],
            {key: "{:.3f}" for key in [
                "seen_accuracy", "unseen_accuracy", "gap_accuracy", "seen_macro_f1",
                "unseen_macro_f1", "gap_macro_f1", "seen_balanced_accuracy",
                "unseen_balanced_accuracy", "gap_balanced_accuracy",
            ]},
        ),
        "",
        "`all_test` is the overall seen/unseen reference across bins and is not the primary "
        "similarity-controlled comparison.",
        "",
        "## 7. Bootstrap uncertainty",
        "",
        "The main accuracy estimates and unseen-minus-seen accuracy gaps, including 95% intervals, "
        "are shown in `accuracy_by_similarity_bin.png` and `accuracy_gap_by_similarity_bin.png`. "
        "All group and gap intervals for accuracy, macro-F1, and balanced accuracy—including "
        "`all_test`—are in `bootstrap_intervals.tsv`.",
        "",
        f"Across {n_total} experiment/bin accuracy gaps, {n_negative} point estimates are negative; "
        f"{n_ci_below} intervals lie entirely below zero and {n_ci_above} entirely above zero. "
        "These counts summarize the pattern and are not multiplicity-adjusted hypothesis tests.",
        "",
        "## 8. EC3 class-composition differences",
        "",
        "The full true-EC3 class counts and frequencies for each group/bin, together with total "
        "variation distance (TV = one half the sum of absolute frequency differences), are in "
        "`class_distribution_tv.tsv`. Larger TV indicates more different class mixtures; it is "
        "not itself a performance measure.",
        "",
        _md_table(
            tv_summary.loc[tv_summary["similarity_bin"].isin(BINS)],
            ["held_out_phylum", "similarity_bin", "seen_n", "unseen_n",
             "seen_ec3_classes", "unseen_ec3_classes", "shared_ec3_classes", "total_variation_distance"],
            {"total_variation_distance": "{:.3f}"},
        ),
        "",
        "## 9. Shared-class results",
        "",
        "Both groups are restricted to proteins whose true EC3 class occurs in both groups "
        "within that experiment/bin. The class list is held fixed for shared-class macro-F1; "
        "out-of-shared-set predictions count as incorrect for accuracy but do not add labels "
        "to shared-class macro-F1. Small retained samples or few shared classes are weak evidence.",
        "",
        _md_table(
            shared.loc[shared["similarity_bin"].isin(BINS)],
            ["held_out_phylum", "similarity_bin", "n_shared_classes", "seen_n", "unseen_n",
             "seen_accuracy", "unseen_accuracy", "gap_accuracy",
             "seen_macro_f1", "unseen_macro_f1", "gap_macro_f1"],
            {key: "{:.3f}" for key in [
                "seen_accuracy", "unseen_accuracy", "gap_accuracy",
                "seen_macro_f1", "unseen_macro_f1", "gap_macro_f1",
            ]},
        ),
        "",
        "## 10. Class-reweighted results",
        "",
        "Only classes present in both groups are retained. For class c, each seen observation "
        "receives per-observation weight `w_i = (n_unseen,c / N_unseen,common) / n_seen,c`. "
        "Thus all seen observations in class c carry total weight equal to the unseen class "
        "proportion, the weights sum to one, and the seen class mix matches unseen over common "
        "support. Unseen accuracy is computed on the same "
        "common-support subset. Classes absent from either group cannot be matched by weighting; "
        "their protein counts and class counts are reported explicitly.",
        "",
        _md_table(
            reweighted.loc[reweighted["similarity_bin"].isin(BINS)],
            ["held_out_phylum", "similarity_bin", "common_support_class_count",
             "seen_n", "unseen_n", "seen_excluded_n", "unseen_excluded_n",
             "class_mix_adjusted_seen_accuracy", "unseen_common_support_accuracy",
             "adjusted_gap_unseen_minus_seen"],
            {key: "{:.3f}" for key in [
                "class_mix_adjusted_seen_accuracy", "unseen_common_support_accuracy",
                "adjusted_gap_unseen_minus_seen",
            ]},
        ),
        "",
        "## 11. `no_hit` results",
        "",
        "`no_hit` is kept as its own category and is not assigned identity zero. It means no "
        "MMseqs2 training alignment passed the qualifying search criteria; it does not mean zero "
        "sequence identity. Its counts and model performance are included separately in all "
        "bin-wise tables and plots.",
        "",
        "## 12. Lab-tested sensitivity analysis",
        "",
        "This repeats bin-wise accuracy for rows with `lab_tested=yes`; it is a sensitivity "
        "analysis, not the primary result. Wilson 95% intervals are descriptive. Counts below "
        "100 are flagged as small-support cells and should not be overinterpreted.",
        "",
        _md_table(
            lab,
            ["held_out_phylum", "similarity_bin", "group", "n_proteins", "accuracy",
             "accuracy_ci_low", "accuracy_ci_high"],
            {key: "{:.3f}" for key in ["accuracy", "accuracy_ci_low", "accuracy_ci_high"]},
        ),
        "",
        ("Cells with fewer than 100 lab-tested proteins: " +
         (", ".join(
             f"{r.held_out_phylum}/{r.similarity_bin}/{r.group} (n={r.n_proteins})"
             for r in lab_small.itertuples(index=False)
         ) if len(lab_small) else "none.")),
        "",
        "## 13. Nearest-homolog baseline",
        "",
        "Baseline 1 is the existing nearest-training-accession EC3 label transfer; no baseline "
        "model was created. `accuracy_all_test_missing_as_wrong` counts missing predictions as "
        "incorrect, consistent with the upstream logistic-regression README. Available-only "
        "accuracy is also shown. Missing values are counted by bin and their reason is in "
        "`validation.json`.",
        "",
        _md_table(
            baseline_overview,
            ["held_out_phylum", "similarity_bin", "group", "n_proteins",
             "n_baseline_predictions", "n_missing_baseline_predictions",
             "accuracy_all_test_missing_as_wrong", "accuracy_when_baseline_available"],
            {key: "{:.3f}" for key in [
                "accuracy_all_test_missing_as_wrong", "accuracy_when_baseline_available",
            ]},
        ),
        "",
        f"Across numerical-bin/group cells, nearest-homolog accuracy is often near ceiling at "
        f"high identity and exceeds the logistic model in many hit-bearing cells (unweighted "
        f"cell means: {baseline_hit['accuracy_all_test_missing_as_wrong'].mean():.3f} vs "
        f"{model_hit['accuracy'].mean():.3f}). The logistic model improves on nearest-neighbor "
        "transfer in some low-identity cells and on every `no_hit` cell; baseline accuracy is "
        f"zero for all {int(no_hit_baseline['n_proteins'].sum()):,} `no_hit` proteins because no "
        "donor prediction exists. The available-only baseline accuracy omits these rows and "
        "must not be compared to all-test logistic accuracy without noting its reduced coverage.",
        "",
        "## 14. Interpretation",
        "",
        f"Overall seen/unseen accuracy gaps are "
        f"{overall.set_index('held_out_phylum').loc['Actinomycetota', 'gap_accuracy']:+.3f} "
        f"for Actinomycetota, "
        f"{overall.set_index('held_out_phylum').loc['Bacillota', 'gap_accuracy']:+.3f} "
        f"for Bacillota, and "
        f"{overall.set_index('held_out_phylum').loc['Cyanobacteriota', 'gap_accuracy']:+.3f} "
        f"for Cyanobacteriota. Within strata, {n_negative} of {n_total} accuracy gaps are "
        f"negative; {n_ci_below} intervals lie wholly below zero and {n_ci_above} wholly above "
        "zero. The direction is not consistent across phyla or bins:",
        "",
        *numerical_patterns,
        "",
        "Some 30-50% bins favor unseen, while other matched bins favor seen. The `no_hit` "
        "direction also reverses across phyla. Class reweighting shrinks several differences "
        "toward zero, but residual gaps remain in both directions and the `no_hit` results "
        "change materially. Shared-class and reweighted analyses therefore do not support a "
        "single consistent unseen penalty. Similarity bins are broad categorical controls, not "
        "exact matching, and can retain within-bin identity differences. EC3 composition differs "
        "substantially across many strata. Even a persistent gap is an association under this "
        "evaluation design, not proof that taxonomy causes it.",
        "",
        "## 15. Limitations",
        "",
        "- Protein-level resampling does not account for dependence among homologous sequences; "
        "comparable sequence-cluster IDs are unavailable for the unseen group.",
        "- Similarity bins are coarse, including the broad 70-100% interval; within-bin similarity "
        "may differ.",
        "- Macro-F1 and balanced accuracy are sensitive to class support and composition; do not "
        "use macro-F1 alone to infer a taxonomic effect.",
        "- Bootstrap intervals are descriptive and not adjusted for multiple bins, metrics, or "
        "experiments. Shared-class/reweighted analyses can exclude substantial numbers of proteins.",
        "- `no_hit` reflects a search qualification outcome, not absence of biological homology.",
        "- The lab-tested flag is analyzed as supplied and does not guarantee independent "
        "experimental validation for every label.",
        "",
        "## 16. Exact output files generated",
        "",
    ]
    for filename in [
        "report.md", "binwise_metrics.tsv", "binwise_gaps.tsv",
        "class_distribution_tv.tsv", "shared_class_metrics.tsv",
        "reweighted_metrics.tsv", "bootstrap_intervals.tsv", "lab_tested_metrics.tsv",
        "baseline1_binwise_metrics.tsv", "validation.json", "settings.json",
        "accuracy_by_similarity_bin.png", "accuracy_by_similarity_bin.pdf",
        "accuracy_gap_by_similarity_bin.png", "accuracy_gap_by_similarity_bin.pdf",
    ]:
        lines.append(f"- `{filename}`")
    lines.append("")
    return "\n".join(lines)


def run_analysis(output, replicates, seed):
    if output.exists() or output.parent != (ROOT / "analysis").resolve():
        raise ValueError("Use a new, direct subdirectory of the project's analysis directory")
    if replicates < 2000:
        raise ValueError("Use at least 2,000 bootstrap replicates")

    count_table = read_tsv(ROOT / "analysis" / "similarity_bins" / "final_bin_counts.tsv")
    required_count_cols = {"held_out_phylum", "group", "similarity_bin", "count"}
    if not required_count_cols.issubset(count_table.columns):
        raise ValueError("Frozen final-bin count table lacks required columns")
    count_table["count"] = pd.to_numeric(count_table["count"], errors="raise")

    inputs = {}
    experiment_data = {}
    validation = {"experiments": {}, "cluster_bootstrap_used": False}
    for phylum in PHYLA:
        prediction_path = ROOT / "analysis" / "logistic_regression" / phylum / "predictions.tsv"
        bins_path = ROOT / "analysis" / "similarity_bins" / f"{phylum}_similarity_bins.tsv"
        split_path = ROOT / "splits_ec3" / phylum / f"split_{phylum}.tsv"
        cluster_metadata_path = ROOT / "splits_ec3" / phylum / f"proteins_split_{phylum}.tsv"
        for path in [prediction_path, bins_path, split_path, cluster_metadata_path]:
            if not path.is_file():
                raise FileNotFoundError(path)
            inputs[str(path.relative_to(ROOT))] = sha256(path)
        predictions = read_tsv(prediction_path)
        validation["experiments"][phylum] = _validate_experiment(
            phylum, predictions, bins_path, split_path, count_table
        )
        cluster_metadata = read_tsv(cluster_metadata_path)
        if not {"accession", "group", "cluster50"}.issubset(cluster_metadata.columns):
            raise ValueError(f"{phylum}: cluster metadata lacks expected columns")
        cluster_test = cluster_metadata.loc[cluster_metadata["group"].isin(GROUPS)]
        if cluster_test["accession"].duplicated().any():
            raise ValueError(f"{phylum}: cluster metadata has duplicate test accessions")
        test_cluster_map = predictions[["accession", "group"]].merge(
            cluster_test[["accession", "group", "cluster50"]],
            on=["accession", "group"], how="left", validate="one_to_one",
        )
        if test_cluster_map["cluster50"].isna().any() or test_cluster_map["cluster50"].eq("").any():
            missing_by_group = {
                group: int((test_cluster_map["group"].eq(group) &
                            test_cluster_map["cluster50"].eq("")).sum())
                for group in GROUPS
            }
        else:
            missing_by_group = {group: 0 for group in GROUPS}
        unseen_cluster_values = test_cluster_map.loc[
            test_cluster_map["group"].eq("unseen_test"), "cluster50"
        ]
        validation["experiments"][phylum]["cluster_ids_missing_by_group"] = missing_by_group
        validation["experiments"][phylum]["unseen_cluster_ids_unavailable"] = bool(
            unseen_cluster_values.eq("").any()
        )
        experiment_data[phylum] = predictions

    validation["all_upstream_inputs_unchanged"] = None
    validation["source_hashes_before"] = inputs
    validation["cluster_bootstrap_unavailable_reason"] = (
        "The frozen split construction does not cluster the held-out phylum; unseen_test "
        "proteins have no cluster50 assignment. Protein-level bootstrap used for both groups."
    )

    rng = np.random.default_rng(seed)
    binwise_rows = []
    gap_rows = []
    tv_rows = []
    tv_summary_rows = []
    shared_rows = []
    reweighted_rows = []
    interval_rows = []
    lab_rows = []
    baseline_rows = []

    for phylum in PHYLA:
        data = experiment_data[phylum]
        for similarity_bin in ALL_BINS:
            subsets = {group: _group_subset(data, group, similarity_bin) for group in GROUPS}
            metric_results = {}
            bootstrap_results = {}
            for group in GROUPS:
                subset = subsets[group]
                metric_results[group] = calculate_metrics(subset)
                binwise_rows.append(_metric_row(
                    phylum, similarity_bin, group, subset, metric_results[group]
                ))
                bootstrap_results[group] = bootstrap_metrics(subset, replicates, rng)

                for metric in ("accuracy", "macro_f1", "balanced_accuracy"):
                    low, high = percentile_interval(bootstrap_results[group][metric])
                    interval_rows.append({
                        "held_out_phylum": phylum,
                        "similarity_bin": similarity_bin,
                        "group": group,
                        "metric": metric,
                        "estimate": metric_results[group][metric],
                        "ci_low": low,
                        "ci_high": high,
                        "confidence_level": CI_LEVEL,
                        "replicates": replicates,
                        "seed": seed,
                    })

            tv_detail, tv, shared_class_count = _class_distribution_tv(
                phylum, similarity_bin, subsets["seen_test"], subsets["unseen_test"]
            )
            tv_rows.extend(tv_detail)
            tv_summary_rows.append({
                "held_out_phylum": phylum,
                "similarity_bin": similarity_bin,
                "seen_n": len(subsets["seen_test"]),
                "unseen_n": len(subsets["unseen_test"]),
                "seen_ec3_classes": metric_results["seen_test"]["n_true_classes"],
                "unseen_ec3_classes": metric_results["unseen_test"]["n_true_classes"],
                "shared_ec3_classes": shared_class_count,
                "total_variation_distance": tv,
            })
            shared_classes = (
                set(subsets["seen_test"]["true_ec3"]) &
                set(subsets["unseen_test"]["true_ec3"])
            )
            shared_result = _shared_metrics(
                subsets["seen_test"], subsets["unseen_test"], shared_classes
            )
            shared_rows.append({
                "held_out_phylum": phylum,
                "similarity_bin": similarity_bin,
                **shared_result,
            })
            reweighted_rows.append({
                "held_out_phylum": phylum,
                "similarity_bin": similarity_bin,
                **_reweighted_accuracy(subsets["seen_test"], subsets["unseen_test"]),
            })

            gap_row = {
                "held_out_phylum": phylum,
                "similarity_bin": similarity_bin,
                "seen_n": len(subsets["seen_test"]),
                "unseen_n": len(subsets["unseen_test"]),
                "seen_ec3_classes": metric_results["seen_test"]["n_true_classes"],
                "unseen_ec3_classes": metric_results["unseen_test"]["n_true_classes"],
                "shared_ec3_classes": shared_class_count,
                "total_variation_distance": tv,
            }
            for metric in ("accuracy", "macro_f1", "balanced_accuracy"):
                seen_value = metric_results["seen_test"][metric]
                unseen_value = metric_results["unseen_test"][metric]
                gap = unseen_value - seen_value
                seen_boot = bootstrap_results["seen_test"][metric]
                unseen_boot = bootstrap_results["unseen_test"][metric]
                gap_low, gap_high = percentile_interval(unseen_boot - seen_boot)
                gap_row[f"seen_{metric}"] = seen_value
                gap_row[f"unseen_{metric}"] = unseen_value
                gap_row[f"gap_{metric}"] = gap
                interval_rows.append({
                    "held_out_phylum": phylum,
                    "similarity_bin": similarity_bin,
                    "group": "gap",
                    "metric": metric,
                    "estimate": gap,
                    "ci_low": gap_low,
                    "ci_high": gap_high,
                    "confidence_level": CI_LEVEL,
                    "replicates": replicates,
                    "seed": seed,
                })
            gap_rows.append(gap_row)

            for group in GROUPS:
                lab_subset = subsets[group].loc[subsets[group]["lab_tested"].eq("yes")]
                n_lab = len(lab_subset)
                n_correct = int(lab_subset["correct"].eq("1").sum())
                ci_low, ci_high = wilson_interval(n_correct, n_lab)
                lab_rows.append({
                    "held_out_phylum": phylum,
                    "similarity_bin": similarity_bin,
                    "group": group,
                    "n_proteins": n_lab,
                    "accuracy": n_correct / n_lab if n_lab else np.nan,
                    "accuracy_ci_low": ci_low,
                    "accuracy_ci_high": ci_high,
                    "small_support_below_100": n_lab < 100,
                })
                baseline = _baseline_metrics(subsets[group])
                baseline_rows.append({
                    "held_out_phylum": phylum,
                    "similarity_bin": similarity_bin,
                    "group": group,
                    **baseline,
                })

    binwise = pd.DataFrame(binwise_rows)
    gaps = pd.DataFrame(gap_rows)
    tv_detail = pd.DataFrame(tv_rows)
    tv_summary = pd.DataFrame(tv_summary_rows)
    shared = pd.DataFrame(shared_rows)
    reweighted = pd.DataFrame(reweighted_rows)
    intervals = pd.DataFrame(interval_rows)
    lab = pd.DataFrame(lab_rows)
    baseline = pd.DataFrame(baseline_rows)

    output.mkdir(parents=False)
    binwise.to_csv(output / "binwise_metrics.tsv", sep="\t", index=False)
    gaps.to_csv(output / "binwise_gaps.tsv", sep="\t", index=False)
    tv_detail.to_csv(output / "class_distribution_tv.tsv", sep="\t", index=False)
    shared.to_csv(output / "shared_class_metrics.tsv", sep="\t", index=False)
    reweighted.to_csv(output / "reweighted_metrics.tsv", sep="\t", index=False)
    intervals.to_csv(output / "bootstrap_intervals.tsv", sep="\t", index=False)
    lab.to_csv(output / "lab_tested_metrics.tsv", sep="\t", index=False)
    baseline.to_csv(output / "baseline1_binwise_metrics.tsv", sep="\t", index=False)
    _make_plots(binwise, intervals, output)

    after_hashes = {path: sha256(ROOT / path) for path in inputs}
    unchanged = after_hashes == inputs
    validation["all_upstream_inputs_unchanged"] = unchanged
    validation["source_hashes_after"] = after_hashes
    validation["all_passed"] = unchanged and all(
        all(row["reconciled"] for row in exp["bin_counts"])
        for exp in validation["experiments"].values()
    )
    if not validation["all_passed"]:
        raise RuntimeError("Validation failed; see computed source hashes/counts")

    settings = {
        "experiments": PHYLA,
        "groups": GROUPS,
        "frozen_bins": BINS,
        "bin_order": BINS,
        "gap_sign_convention": "unseen_test minus seen_test",
        "bootstrap_method": "ordinary protein-level nonparametric bootstrap",
        "bootstrap_replicates": replicates,
        "bootstrap_seed": seed,
        "confidence_level": CI_LEVEL,
        "bootstrap_interval": "percentile",
        "cluster_bootstrap_used": False,
        "cluster_bootstrap_limitation": validation["cluster_bootstrap_unavailable_reason"],
        "macro_f1_labels": "union of true and predicted labels within group/bin; fixed within each bootstrap",
        "balanced_accuracy_labels": "classes observed in true labels for each evaluation sample",
        "class_mix_weight_per_seen_observation": "(unseen class count / unseen common-support total) / seen class count",
        "class_mix_unseen_comparison": "unseen accuracy restricted to the same common-support classes",
        "input_sha256": inputs,
        "script_sha256": sha256(Path(__file__)),
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "matplotlib_version": matplotlib.__version__,
    }
    (output / "settings.json").write_text(json.dumps(settings, indent=2) + "\n", encoding="utf-8")
    (output / "validation.json").write_text(
        json.dumps(validation, indent=2) + "\n", encoding="utf-8"
    )
    report = _build_report(
        binwise, gaps, tv_summary, shared, reweighted, intervals, lab, baseline,
        validation, settings,
    )
    (output / "report.md").write_text(report, encoding="utf-8")
    print(f"Analysis complete: {output}")
    print(f"Bootstrap: {replicates} replicates; seed {seed}")
    print(f"Upstream inputs unchanged: {unchanged}")
    print(f"Accuracy gaps below zero: {int(gaps.loc[gaps['similarity_bin'].isin(BINS), 'gap_accuracy'].lt(0).sum())}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--bootstrap-replicates", type=int, default=BOOTSTRAP_REPLICATES)
    parser.add_argument("--seed", type=int, default=BOOTSTRAP_SEED)
    args = parser.parse_args()
    run_analysis(args.output.resolve(), args.bootstrap_replicates, args.seed)


if __name__ == "__main__":
    main()
