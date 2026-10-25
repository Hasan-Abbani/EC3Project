"""Assign common, outcome-independent similarity strata to frozen test records.

Run from any directory: python code/make_similarity_bins.py
Inputs are the completed nearest-training analysis and existing final splits.
Only identities, test-group membership, and sample support select boundaries.
No classifiers, embeddings, or performance outputs are read.
"""
import argparse
from bisect import bisect_right
from decimal import Decimal
import hashlib
import itertools
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT / "analysis" / "nearest_train_identity"
PHYLA = ["Actinomycetota", "Bacillota", "Cyanobacteriota"]
GROUPS = ["seen_test", "unseen_test"]
MIN_SUPPORT = 100  # Fixed before candidate evaluation; per group AND experiment.
BASE = [0, 20, 30, 40, 50, 60, 70, 80, 90, 100]


def read(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def labels(edges):
    return [f"<{b}%" if a == 0 else f"{a}-{b}%" for a, b in zip(edges, edges[1:])]


def assign(values, statuses, edges):
    """Use decimal comparisons to avoid ambiguity at exact boundary values."""
    boundaries = [Decimal(str(e)) / 100 for e in edges]
    names = labels(edges)
    result = []
    for value, status in zip(values, statuses):
        if status == "no_qualifying_hit":
            if value != "":
                raise ValueError("No-hit record has a numerical identity")
            result.append("no_hit")
            continue
        if status != "hit" or value == "":
            raise ValueError("Invalid hit status or missing identity for hit")
        identity = Decimal(value)
        if not identity.is_finite() or not Decimal(0) <= identity <= Decimal(1):
            raise ValueError("Identity must be finite and in [0,1]")
        index = min(bisect_right(boundaries, identity) - 1, len(names) - 1)
        result.append(names[index])
    return pd.Series(result, index=values.index, name="similarity_bin", dtype=str)


def evaluate(frame, name, edges):
    bins = assign(frame.max_train_identity, frame.hit_status, edges)
    counts = (frame.assign(similarity_bin=bins)
              .groupby(["held_out_phylum", "similarity_bin", "group"]).size())
    rows = []
    for p in PHYLA:
        for i, label in enumerate([*labels(edges), "no_hit"]):
            n = [int(counts.get((p, label, g), 0)) for g in GROUPS]
            denominators = [int((frame.held_out_phylum.eq(p) & frame.group.eq(g)).sum()) for g in GROUPS]
            rows.append(dict(scheme=name, held_out_phylum=p, similarity_bin=label,
                             lower_percent=edges[i] if label != "no_hit" else None,
                             upper_percent=edges[i+1] if label != "no_hit" else None,
                             seen_test=n[0], unseen_test=n[1], total=sum(n),
                             proportion_seen_test=n[0]/denominators[0],
                             proportion_unseen_test=n[1]/denominators[1],
                             numerical=label != "no_hit",
                             support_rule_passed=min(n) >= MIN_SUPPORT if label != "no_hit" else None))
    return pd.DataFrame(rows)


def select_common_edges(base_counts):
    """Find the finest supported adjacent coarsening of the declared grid.

    Optimize resolution, never similarity of the group distributions. Ties
    prefer smaller sum of squared widths, then lexicographic boundaries.
    """
    cells = np.array([[int(base_counts.loc[
        base_counts.held_out_phylum.eq(p) & base_counts.similarity_bin.eq(b), g].iloc[0])
        for p in PHYLA for g in GROUPS] for b in labels(BASE)])
    choices = []
    audit = []
    for mask in itertools.product([False, True], repeat=len(BASE)-2):
        indices = [0, *[i+1 for i, keep in enumerate(mask) if keep], len(BASE)-1]
        edges = [BASE[i] for i in indices]
        grouped = np.array([cells[a:b].sum(axis=0) for a, b in zip(indices, indices[1:])])
        minimum = int(grouped.min())
        passed = minimum >= MIN_SUPPORT
        audit.append(dict(boundaries_percent=json.dumps(edges), numerical_bins=len(edges)-1,
                          smallest_group_count=minimum, all_groups_supported=passed))
        if passed:
            choices.append(edges)
    if not choices:
        raise ValueError("Even the full identity range fails common support")
    selected = min(choices, key=lambda e: (-(len(e)-1), sum((b-a)**2 for a,b in zip(e,e[1:])), e))
    return selected, pd.DataFrame(audit)


def markdown_table(frame, columns):
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"]*len(columns)) + " |"]
    for row in frame[columns].itertuples(index=False, name=None):
        lines.append("| " + " | ".join(str(v) for v in row) + " |")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "analysis" / "similarity_bins")
    args = parser.parse_args()
    out = args.output.resolve()
    # Fail before writing if output overlaps any source or already exists.
    if out.exists() or out.parent != (ROOT / "analysis").resolve():
        raise ValueError("Use a new direct subdirectory of the project's analysis directory")
    paths = [SOURCE / n for n in ["summary.tsv", "distribution_overlap.tsv", "report.md",
             "settings.json", "validation.json", "independent_validation.json", "identity_distributions.pdf"]]
    for p in PHYLA:
        paths.extend([SOURCE/p/"test_nearest_train_identity.tsv",
                      *[SOURCE/p/f"{g}_nearest_train_identity.tsv" for g in GROUPS],
                      ROOT/"splits_ec3"/p/f"split_{p}.tsv"])
    hashes = {str(p.relative_to(ROOT)): sha256(p) for p in paths}
    summary = read(SOURCE / "summary.tsv").set_index(["held_out_phylum", "group"])
    frames = []
    for p in PHYLA:
        data = read(SOURCE/p/"test_nearest_train_identity.tsv")
        if not data.held_out_phylum.eq(p).all() or set(data.group) != set(GROUPS):
            raise ValueError("Unexpected experiment or test group")
        if data.accession.duplicated().any():
            raise ValueError("Duplicate accession within experiment")
        compact = read(ROOT/"splits_ec3"/p/f"split_{p}.tsv")
        expected = compact[compact.group.isin(GROUPS)].set_index("accession").group.sort_index()
        pd.testing.assert_series_equal(data.set_index("accession").group.sort_index(), expected)
        for g in GROUPS:
            sub = data[data.group.eq(g)].reset_index(drop=True)
            pd.testing.assert_frame_equal(sub, read(SOURCE/p/f"{g}_nearest_train_identity.tsv"))
            if len(sub) != int(summary.loc[(p, g), "expected"]):
                raise ValueError("Summary total disagrees")
            if sub.hit_status.eq("no_qualifying_hit").sum() != int(summary.loc[(p, g), "no_hit"]):
                raise ValueError("Summary missing count disagrees")
        frames.append(data)
    data = pd.concat(frames, ignore_index=True)
    if len(data) != 84762:
        raise ValueError("Expected 84,762 experiment-specific test records")
    # Selection reads only identities, groups, and experiment identifiers.
    base_counts = evaluate(data, "ten_percentage_points_with_below20", BASE)
    chosen, coarsenings = select_common_edges(base_counts)
    schemes = {
        "ten_percentage_points_with_below20": BASE,
        "twenty_percentage_points_with_below20": [0,20,40,60,80,100],
        "coarse_twenty_with_merged_tails": [0,40,60,100],
        "central_tens_broad_tails": [0,30,40,50,60,100],
        "selected_common_adjacent_merges": chosen,
    }
    candidates = pd.concat([evaluate(data, name, edges) for name, edges in schemes.items()], ignore_index=True)
    decisions = []
    for name, edges in schemes.items():
        numerical = candidates[candidates.scheme.eq(name) & candidates.numerical]
        decisions.append(dict(scheme=name, boundaries_percent=json.dumps(edges), numerical_bins=len(edges)-1,
                              failing_experiment_bins=int(numerical.support_rule_passed.eq(False).sum()),
                              smallest_group_count=int(numerical[["seen_test", "unseen_test"]].min().min()),
                              all_groups_supported=bool(numerical.support_rule_passed.all()),
                              selected=name == "selected_common_adjacent_merges"))
    decisions = pd.DataFrame(decisions)
    final = candidates[candidates.scheme.eq("selected_common_adjacent_merges")].copy()
    if not final.loc[final.numerical, "support_rule_passed"].all():
        raise ValueError("Final numerical bins fail support")
    assigned = data.assign(similarity_bin=assign(data.max_train_identity, data.hit_status, chosen))
    pd.testing.assert_frame_equal(assigned[data.columns], data)
    if assigned.duplicated(["held_out_phylum", "accession"]).any() or assigned.similarity_bin.isna().any():
        raise ValueError("Lost or duplicated bin assignments")
    out.mkdir(parents=False)
    candidates.to_csv(out/"candidate_bin_counts.tsv", sep="\t", index=False)
    decisions.to_csv(out/"candidate_schemes.tsv", sep="\t", index=False)
    coarsenings.to_csv(out/"adjacent_coarsening_evaluation.tsv", sep="\t", index=False)
    final.to_csv(out/"final_bin_counts_wide.tsv", sep="\t", index=False)
    validation, long_counts, low_counts = [], [], []
    for p in PHYLA:
        sub = assigned[assigned.held_out_phylum.eq(p)].reset_index(drop=True)
        path = out / f"{p}_similarity_bins.tsv"
        sub.to_csv(path, sep="\t", index=False)
        pd.testing.assert_frame_equal(read(path)[data.columns], data[data.held_out_phylum.eq(p)].reset_index(drop=True))
        for g in GROUPS:
            group = sub[sub.group.eq(g)]
            hits = group.hit_status.eq("hit")
            if not group.similarity_bin.eq("no_hit").equals(~hits):
                raise ValueError("No-hit records mixed into numerical strata")
            counts = group.similarity_bin.value_counts()
            if counts.sum() != len(group):
                raise ValueError("Assignment totals disagree")
            for b in [*labels(chosen), "no_hit"]:
                n = int(counts.get(b, 0))
                long_counts.append(dict(held_out_phylum=p, group=g, similarity_bin=b,
                                        count=n, proportion_of_test_group=n/len(group),
                                        proportion_of_hit_bearing_group=n/int(hits.sum()) if b != "no_hit" else None))
            values = pd.to_numeric(group.max_train_identity, errors="coerce")
            low_counts.append(dict(held_out_phylum=p, group=g, qualifying_below20=int(values.lt(.2).sum())))
            validation.append(dict(held_out_phylum=p, group=g, expected=len(group),
                                   numerical=int(hits.sum()), no_hit=int((~hits).sum()),
                                   assignment_total=int(counts.sum()), all_original_columns_preserved=True))
    assigned.to_csv(out/"all_test_similarity_bins.tsv", sep="\t", index=False)
    pd.testing.assert_frame_equal(read(out/"all_test_similarity_bins.tsv")[data.columns], data)
    pd.DataFrame(long_counts).to_csv(out/"final_bin_counts.tsv", sep="\t", index=False)
    pd.DataFrame(low_counts).to_csv(out/"below20_counts.tsv", sep="\t", index=False)
    sensitivity = []
    numerical = final[final.numerical]
    for threshold in [50,100,200]:
        sensitivity.append(dict(threshold_per_group=threshold,
            failing_experiment_bins=int(numerical[["seen_test","unseen_test"]].lt(threshold).any(axis=1).sum()),
            failing_group_cells=int(numerical[["seen_test","unseen_test"]].lt(threshold).sum().sum())))
    pd.DataFrame(sensitivity).to_csv(out/"support_sensitivity.tsv", sep="\t", index=False)
    fig, axes = plt.subplots(3,1,figsize=(10,10),layout="constrained")
    names = [*labels(chosen), "no_hit"]
    for ax, p in zip(axes, PHYLA):
        table = final[final.held_out_phylum.eq(p)].set_index("similarity_bin").loc[names]
        pos = np.arange(len(names))
        for offset, g, color in [(-.2,"seen_test","#0072B2"),(.2,"unseen_test","#D55E00")]:
            bars = ax.bar(pos+offset,table[g],width=.38,label=g,color=color)
            bars[-1].set_hatch("///")
            ax.bar_label(bars,fontsize=7,padding=2,fmt="%d")
        ax.set(title=p, ylabel="Test proteins", xticks=pos, xticklabels=names)
        ax.set_ylim(0,table[GROUPS].to_numpy().max()*1.2)
        ax.axvline(len(names)-1.5,color="gray",linestyle="--",linewidth=.8)
        ax.legend(loc="upper right",fontsize=8)
    fig.suptitle("Common similarity strata; no_hit is a separate category",fontsize=12)
    fig.savefig(out/"similarity_bin_counts.png",dpi=180)
    fig.savefig(out/"similarity_bin_counts.pdf")
    plt.close(fig)
    unchanged = all(sha256(ROOT/p) == h for p,h in hashes.items())
    if not unchanged:
        raise RuntimeError("An existing input changed during binning")
    settings = dict(minimum_per_group_per_experiment=MIN_SUPPORT,
        minimum_support_rule_fixed_before_selection=True,
        support_rationale="100 gives worst-case SE 0.05 for a simple independent proportion; a pragmatic floor, not a guarantee for class-specific or clustered estimates.",
        base_grid_percent=BASE, final_edges_percent=chosen, final_labels=labels(chosen),
        closure="Left-closed/right-open; last numerical bin includes 100%; no_hit is separate",
        selection="Finest common supported adjacent coarsening of the fixed grid; ties by squared widths then lexicographic edges",
        candidates=schemes, selection_variables=["max_train_identity","group","held_out_phylum"],
        input_sha256=hashes, script_sha256=sha256(Path(__file__)),
        pandas_version=pd.__version__, numpy_version=np.__version__, matplotlib_version=matplotlib.__version__)
    (out/"settings.json").write_text(json.dumps(settings,indent=2)+"\n",encoding="utf-8")
    (out/"validation.json").write_text(json.dumps(dict(all_passed=True,total_test_records=len(assigned),
        hit_bearing=int(assigned.hit_status.eq("hit").sum()),
        no_hit=int(assigned.hit_status.eq("no_qualifying_hit").sum()),unique_within_each_experiment=True,
        every_hit_assigned_once=True,no_hits_separate=True,original_input_hashes_unchanged=unchanged,
        all_final_numerical_bins_meet_support=True,groups=validation),indent=2)+"\n",encoding="utf-8")
    report = ["# Final sequence-similarity strata", "",
        "Decision made before any classifier or embedding step, using only saved identity, group membership and sample support. No prediction or performance outputs were read.", "",
        "## Support rule fixed before selection", "",
        "Require at least 100 hit-bearing seen-test AND 100 hit-bearing unseen-test proteins in each numerical bin, separately for ALL three experiments. This pragmatic floor gives a worst-case standard error of 5 percentage points for a simple independent proportion (roughly +/-9.8 points at 95% under a normal approximation). It does not guarantee precision for class-specific metrics, small differences between groups, or correlated proteins within sequence clusters. Future uncertainty assessment should account for these limits.", "",
        "The threshold was stated before candidate counts were evaluated. A stricter 200-per-group rule fails three group cells in the selected scheme (see support_sensitivity.tsv), requiring further merging; 100 retains useful resolution without accepting very sparse groups. Counts of 100 are an aggregate floor, not class-level support.", "",
        "## Candidates and decision", "", markdown_table(decisions,["scheme","boundaries_percent","numerical_bins","smallest_group_count","failing_experiment_bins","selected"]), "",
        "The selected scheme is the finest supported adjacent coarsening of the fixed 10-percentage-point grid with a separate below-20 interval. All 256 adjacent coarsenings were evaluated with an explicit resolution-only objective; no criterion rewards making distributions or group sizes look similar. Coarser supported alternatives sacrifice central resolution without being needed for the stated support floor.", "",
        "**Final labels:** " + ", ".join(f"`{n}`" for n in names) + ".", "",
        "Numerical boundaries in fractions: " + json.dumps([e/100 for e in chosen]) + ". Intervals include their lower bound and exclude their upper bound, except the last includes 1.000. Thus exactly 30% is in 30-40%, exactly 70% is in 70-100%, and exactly 100% is retained. The first bin is [0,30%); this includes qualifying low identities, not imputed zeros.", "",
        "## Tail merges", "",
        "- Below 20%: one qualifying unseen protein each in Actinomycetota and Bacillota, and no seen proteins; no such proteins in Cyanobacteriota. Merge with 20-30% into <30%; no record is discarded.",
        "- The 90-100% tail has only 3, 12, and 0 unseen proteins, respectively. Merging 80-100% still gives only 22 Actinomycetota unseen and 32 Cyanobacteriota unseen proteins. Merge these tails with 70-80% across ALL experiments, giving 70-100% and at least 132 per group.",
        "- Central 30-40%, 40-50%, 50-60%, and 60-70% intervals already meet support and remain separate. Final smallest group support is 132.", "",
        "## Final counts", "", markdown_table(final,["held_out_phylum","similarity_bin","seen_test","unseen_test","total"]), "",
        "candidate_bin_counts.tsv gives every candidate/experiment/bin count, both group proportions (denominator: the entire test group including no_hit), totals, and support flags. final_bin_counts.tsv supplies final counts by experiment, group and stratum, with both whole-test-group and hit-bearing proportions. no_hit has no hit-bearing proportion and is excluded from numerical-bin support decisions.", "",
        "## Validation and interpretation", "",
        "All 84,762 experiment-specific test records remain represented: 78,019 in exactly one numerical stratum and 6,743 in no_hit. Uniqueness is checked within experiment, since the same accession can legitimately appear in different held-out experiments. Combined and per-group identity inputs agree, test accessions/groups match frozen split files, assignment totals match each group, every original source column and identity string survives an output round trip unchanged, and all input-file hashes remain unchanged.", "",
        "No-hit identities remain missing, never zero. These proteins must be reported and evaluated separately; no qualifying hit does not imply absence of homology. The upstream identities are approximate maxima among returned MMseqs2 local alignments requiring >=80% coverage of both sequences and E<=0.001, not guaranteed global maxima. Broad bins only coarsely control similarity: within-bin identity distributions and EC3 composition can still differ, and a count floor alone does not establish exchangeability or independent observations. In particular, 70-100% pools a broad high-identity tail. Do not infer comparable identities from equal bin labels alone.", "",
        "No classifier was inspected or trained. The upstream 2% histogram intervals were not used as analysis strata. Existing nearest-identity and split files were not modified. No Git commit or push is performed by this script.", "",
        "## Reproduction and files", "",
        "Run `python code/make_similarity_bins.py` from the project root (or use a new `--output analysis/<name>`). The script refuses to overwrite an existing output directory. Requires pandas, numpy and matplotlib; no MMseqs2 rerun is needed.", "",
        "- all_test_similarity_bins.tsv and <phylum>_similarity_bins.tsv: original per-protein columns plus similarity_bin.",
        "- candidate_schemes.tsv, candidate_bin_counts.tsv, adjacent_coarsening_evaluation.tsv: alternatives and support decisions.",
        "- final_bin_counts.tsv, final_bin_counts_wide.tsv, below20_counts.tsv, support_sensitivity.tsv: final counts and diagnostics.",
        "- similarity_bin_counts.png and .pdf: grouped counts with no_hit shown separately by hatching.",
        "- settings.json, validation.json, report.md: decision rule, provenance and validation."]
    (out/"report.md").write_text("\n".join(report)+"\n",encoding="utf-8")
    print(decisions.to_string(index=False))
    print(final[["held_out_phylum","similarity_bin","seen_test","unseen_test"]].to_string(index=False))
    print("All 84,762 records validated; input files unchanged. Outputs:",out)


if __name__ == "__main__":
    main()
