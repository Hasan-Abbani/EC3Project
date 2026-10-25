"""Nearest-training identity for frozen EC3 splits; never regenerates splits.

Run: python code/nearest_train_identity.py
Requires pandas, numpy, matplotlib and the existing WSL MMseqs2 installation.
Identity is MMseqs2's reported fraction of identical alignment columns in the local
alignment, with >=80% coverage of BOTH sequences and E-value <=1e-3.
The sensitive k-mer search is approximate: maxima are over returned hits.
"""
import argparse
import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PHYLA = ["Actinomycetota", "Bacillota", "Cyanobacteriota"]
GROUPS = ["seen_test", "unseen_test"]
MMSEQS = "/home/hsnab/mmseqs/bin/mmseqs"
FORMAT = "query,target,fident,nident,alnlen,qcov,tcov,evalue,bits"
FIELDS = FORMAT.split(",")


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def wsl_path(path):
    p = Path(path).resolve()
    return f"/mnt/{p.drive[0].lower()}{p.as_posix()[2:]}"


def read(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def command(args, log, commands):
    cmd = ["wsl", MMSEQS, *map(str, args)]
    commands.append(cmd)
    print("Running:", " ".join(cmd), flush=True)
    with log.open("w", encoding="utf-8") as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=True)


def fasta(frame, path):
    # A rerun reuses only byte-identical FASTAs, never a stale representation.
    content = "".join(f">{a}\n{s}\n" for a, s in zip(frame.accession, frame.sequence))
    if path.exists():
        if path.read_text(encoding="utf-8") != content:
            raise ValueError(f"Existing FASTA differs from final split: {path}")
    else:
        path.write_text(content, encoding="utf-8", newline="\n")


def load_split(phylum):
    folder = ROOT / "splits_ec3" / phylum
    full = folder / f"proteins_split_{phylum}.tsv"
    compact = folder / f"split_{phylum}.tsv"
    x = read(full)
    pd.testing.assert_frame_equal(x[["accession", "cluster50", "group"]], read(compact))
    if x.accession.duplicated().any() or x.sequence.eq("").any():
        raise ValueError("Duplicate accession or missing sequence")
    if set(x.group) != {"train", "validation", *GROUPS}:
        raise ValueError("Unexpected split groups")
    if not x.sequence.str.len().eq(x.length.astype(int)).all():
        raise ValueError("Sequence length mismatch")
    if not x.phylum.eq(phylum).equals(x.group.eq("unseen_test")):
        raise ValueError("Held-out phylum assignment mismatch")
    seen = x[x.group.ne("unseen_test")]
    if seen.groupby("cluster50").group.nunique().max() != 1:
        raise ValueError("Cluster crosses split groups")
    train = x[x.group.eq("train")]
    if not train.ec3.value_counts().ge(10).all():
        raise ValueError("Training class threshold violated")
    if set(x.ec3) != set(train.ec3):
        raise ValueError("A retained class is absent from training")
    if set(read(folder / "allowed_ec3_classes.tsv").ec3) != set(train.ec3):
        raise ValueError("Allowed-class file disagrees")
    if set(train.sequence) & set(x[x.group.isin(GROUPS)].sequence):
        raise ValueError("Exact sequence occurs in train and test")
    summary = read(folder / "split_summary.tsv").set_index("group")
    for g, sub in x.groupby("group"):
        if len(sub) != int(summary.loc[g, "proteins"]):
            raise ValueError("Split summary count mismatch")
    return x


def reduce_hits(path, expected, train_ids):
    best = pd.DataFrame(columns=FIELDS)
    counts = pd.Series(0, index=expected, dtype="int64")
    first = {}
    with pd.read_csv(path, sep="\t", names=FIELDS, chunksize=200000,
                     dtype={"query": str, "target": str}) as chunks:
        for chunk in chunks:
            if chunk.empty:
                continue
            if not set(chunk["query"]).issubset(set(expected)):
                raise ValueError("Unexpected query accession in hits")
            if not set(chunk.target).issubset(train_ids):
                raise ValueError("Hit target is not in final training set")
            numeric = chunk[FIELDS[2:]]
            if not np.isfinite(numeric.to_numpy()).all():
                raise ValueError("Nonfinite alignment output")
            if not (chunk.alnlen.gt(0) & chunk.nident.ge(0)
                    & chunk.nident.le(chunk.alnlen) & chunk.fident.between(0, 1)
                    & chunk.qcov.ge(0.799) & chunk.tcov.ge(0.799)
                    & chunk.evalue.le(0.001001)).all():
                raise ValueError("Invalid identity/coverage/E-value in returned hit")
            if not np.allclose(chunk.fident, chunk.nident / chunk.alnlen, atol=0.0011):
                raise ValueError("Reported identity disagrees with identical-residue count")
            # Maximize reported fident, without inferring extra precision from nident.
            for q, v in chunk.drop_duplicates("query").set_index("query").fident.items():
                first.setdefault(q, v)
            counts = counts.add(chunk.groupby("query").size(), fill_value=0).astype("int64")
            candidates = pd.concat([best, chunk], ignore_index=True)
            best = (candidates.sort_values(["fident", "bits", "target"],
                                          ascending=[False, False, True])
                    .drop_duplicates("query"))
    best = best.rename(columns={"query": "accession", "target": "nearest_train_accession",
                                "fident": "max_train_identity"})
    lower_first = sum(first[a] < v - 1e-12
                      for a, v in zip(best.accession, best.max_train_identity))
    return best, counts, int(lower_first)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "analysis" / "nearest_train_identity")
    parser.add_argument("--threads", type=int, default=12)
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--pilot", action="store_true", help="Eight test proteins per group, all final training proteins")
    args = parser.parse_args()
    out = args.output.resolve()
    if args.pilot:
        out = out / "pilot"
    if out == ROOT or any(out == ROOT / n or ROOT / n in out.parents
                          for n in ["splits_ec3", "data", "raw", "code"]):
        raise ValueError("Output must be a separate analysis directory")
    # Inspect and validate every final split before creating analysis outputs.
    inputs = sorted((ROOT / "splits_ec3").rglob("*.tsv"))
    inputs += sorted((ROOT / "splits_ec3").rglob("settings.json"))
    hashes = {str(p.relative_to(ROOT)): digest(p) for p in inputs}
    splits = {p: load_split(p) for p in PHYLA}
    cleaned = read(ROOT / "data" / "proteins_ec3.tsv").set_index("accession")
    for p, x in splits.items():
        if not x.set_index("accession").sequence.eq(cleaned.loc[x.accession, "sequence"]).all():
            raise ValueError("Final split sequence disagrees with cleaned dataset")
        print(p, x.group.value_counts().to_dict(), flush=True)
    del cleaned
    if args.validate_only:
        print("All three frozen splits validated; no files changed.")
        return
    out.mkdir(parents=True, exist_ok=True)
    version = subprocess.check_output(["wsl", MMSEQS, "version"], text=True).strip()
    settings = dict(mmseqs=MMSEQS, version=version, coverage=0.8, cov_mode=0,
                    evalue=0.001, sensitivity=7.5, alignment_mode=3, seq_id_mode=0,
                    identity_units="fraction [0,1]", approximation=True,
                    definition="Maximum reported fident over returned local alignments; both coverages >=0.8",
                    format_output=FORMAT, threads=args.threads, pilot=args.pilot, input_sha256=hashes,
                    script_sha256=digest(Path(__file__)), python=platform.python_version())
    config = out / "settings.json"
    if config.exists() and json.loads(config.read_text()) != settings:
        raise ValueError("Existing analysis settings differ; use a fresh --output directory")
    write_json(config, settings)
    commands = []
    summaries, overlaps, validations = [], [], []
    fig, axes = plt.subplots(3, 2, figsize=(12, 12))
    for row, (p, x) in enumerate(splits.items()):
        if args.pilot:
            x = pd.concat([x[x.group.eq("train")],
                           *[x[x.group.eq(g)].head(8) for g in GROUPS]], ignore_index=True)
        folder = out / p
        folder.mkdir(exist_ok=True)
        train = x[x.group.eq("train")]
        train_ids = set(train.accession)
        for g in ["train", *GROUPS]:
            sub = x[x.group.eq(g)]
            fasta(sub, folder / f"{g}.fasta")
            db = folder / f"{g}_db"
            if not db.with_suffix(".dbtype").exists():
                command(["createdb", wsl_path(folder / f"{g}.fasta"), wsl_path(db),
                         "--dbtype", 1, "--shuffle", 0], folder / f"{g}_createdb.log", commands)
        results = []
        for g in GROUPS:
            test = x[x.group.eq(g)]
            hitfile = folder / f"{g}_hits.tsv"
            done = folder / f"{g}_search_complete.json"
            if not done.exists():
                result = folder / f"{g}_alignments"
                command(["search", wsl_path(folder / f"{g}_db"), wsl_path(folder / "train_db"),
                         wsl_path(result), wsl_path(folder / f"{g}_tmp"),
                         "-s", 7.5, "--max-seqs", len(train), "--max-accept", len(train),
                         "--max-rejected", len(train), "--alignment-mode", 3, "-a", 1,
                         "--seq-id-mode", 0, "--min-seq-id", 0, "-c", 0.8, "--cov-mode", 0,
                         "-e", 0.001, "--threads", args.threads, "--split-memory-limit", "2G"],
                        folder / f"{g}_search.log", commands)
                command(["convertalis", wsl_path(folder / f"{g}_db"), wsl_path(folder / "train_db"),
                         wsl_path(result), wsl_path(hitfile), "--format-output", FORMAT,
                         "--threads", args.threads], folder / f"{g}_convertalis.log", commands)
                write_json(done, {"hits_sha256": digest(hitfile)})
            elif json.loads(done.read_text())["hits_sha256"] != digest(hitfile):
                raise ValueError("Cached hit table checksum differs")
            best, counts, lower_first = reduce_hits(hitfile, list(test.accession), train_ids)
            cols = ["accession", "phylum", "ec3", "cluster50", "group"]
            data = test[cols].merge(best, on="accession", how="left", validate="one_to_one")
            data.insert(0, "held_out_phylum", p)
            data["returned_hit_count"] = data.accession.map(counts)
            data["hit_status"] = np.where(data.max_train_identity.isna(), "no_qualifying_hit", "hit")
            data["identity_is_approximate"] = True
            data["query_contains_U_or_O"] = data.accession.map(
                test.set_index("accession").sequence.str.contains("[UO]"))
            data["nearest_train_contains_U_or_O"] = data.nearest_train_accession.map(
                train.set_index("accession").sequence.str.contains("[UO]"))
            if len(data) != len(test) or set(data.accession) != set(test.accession):
                raise ValueError("Expected test accessions not preserved")
            data.to_csv(folder / f"{g}_nearest_train_identity.tsv", sep="\t", index=False)
            missing = data[data.hit_status.eq("no_qualifying_hit")]
            missing[cols + ["hit_status"]].to_csv(folder / f"{g}_no_hits.tsv", sep="\t", index=False)
            v = data.max_train_identity.dropna().astype(float)
            stats = {"held_out_phylum": p, "group": g, "expected": len(test), "with_hit": len(v),
                     "no_hit": len(missing),
                     "queries_with_U_or_O": int(data.query_contains_U_or_O.sum()),
                     "winning_targets_with_U_or_O": int(data.nearest_train_contains_U_or_O.fillna(False).sum()),
                     "min": float(v.min()) if len(v) else None,
                     "median": float(v.median()) if len(v) else None,
                     "mean": float(v.mean()) if len(v) else None, "max": float(v.max()) if len(v) else None}
            summaries.append(stats)
            validations.append(dict(held_out_phylum=p, group=g, all_expected_accessions_present=True,
                                    targets_in_final_train=True, first_hit_below_max=lower_first,
                                    no_qualifying_hit=len(missing), hit_table_sha256=digest(hitfile)))
            print(stats, flush=True)
            axes[row, 0].hist(v, bins=np.linspace(0, 1, 51), weights=np.ones(len(v))/len(v) if len(v) else None,
                              histtype="step", label=f"{g}: {len(v):,} hits; {len(missing):,} missing")
            if len(v):
                axes[row, 1].step(np.sort(v), np.arange(1, len(v)+1)/len(v), label=g)
            results.append(data)
        combined = pd.concat(results, ignore_index=True)
        combined.to_csv(folder / "test_nearest_train_identity.tsv", sep="\t", index=False)
        a, b = [d.max_train_identity.dropna().astype(float) for d in results]
        if len(a) and len(b):
            lo, hi = max(a.min(), b.min()), min(a.max(), b.max())
            grid = np.sort(np.unique(np.concatenate([a, b])))
            ks = np.max(np.abs(np.searchsorted(np.sort(a), grid, side="right")/len(a)
                               - np.searchsorted(np.sort(b), grid, side="right")/len(b)))
            overlaps.append(dict(held_out_phylum=p, common_range_min=float(lo), common_range_max=float(hi),
                                 seen_fraction_in_common_range=float(a.between(lo, hi).mean()),
                                 unseen_fraction_in_common_range=float(b.between(lo, hi).mean()),
                                 ecdf_max_absolute_difference=float(ks)))
        for j, title in enumerate(["Identity distribution (fixed 2% display intervals)", "Empirical CDF"]):
            axes[row, j].set(title=f"{p}: {title}", xlabel="Maximum returned training identity (fraction)",
                             ylabel="Fraction of proteins with hits" if j == 0 else "Cumulative fraction", xlim=(0, 1))
            axes[row, j].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "identity_distributions.png", dpi=180)
    fig.savefig(out / "identity_distributions.pdf")
    plt.close(fig)
    pd.DataFrame(summaries).to_csv(out / "summary.tsv", sep="\t", index=False)
    pd.DataFrame(overlaps).to_csv(out / "distribution_overlap.tsv", sep="\t", index=False)
    unchanged = all(digest(ROOT / p) == h for p, h in hashes.items())
    write_json(out / "validation.json", dict(split_files_unchanged=unchanged, groups=validations,
                                              finished_utc=datetime.now(timezone.utc).isoformat()))
    write_json(out / f"commands_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')}.json", commands)
    if not unchanged:
        raise RuntimeError("An existing split input changed during analysis")
    print("Analysis complete; frozen split hashes unchanged.", flush=True)


if __name__ == "__main__":
    main()
