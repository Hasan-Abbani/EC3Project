"""Proposal baseline 1: amino-acid 2-mer frequencies + multinomial logistic regression.

Identical to Phase E's run_logreg.py (same splits, same train-only standardization,
same LogisticRegression settings, same C grid chosen by validation accuracy, no refit
on train+validation, same predictions.tsv / summary.json format). The ONLY change is
the features: 400 overlapping 2-mer frequencies instead of the 1,280 ESM-2 numbers.

Features
  For a protein of length L there are L-1 overlapping 2-mers. Each of the 400 standard
  pairs (20 x 20 amino acids) is counted and divided by (L-1).
  2-mers that contain a non-standard letter (U, O, X, B, Z, ...) are skipped (not counted),
  but the divisor stays L-1, so such proteins' 400 frequencies sum to slightly below 1.

Usage (from the repository root, after `git lfs pull` so data/proteins_ec3.tsv is present)
  python code/run_kmer_logreg.py --phylum Actinomycetota
  Defaults: --split splits_ec3/<phylum>/split_<phylum>.tsv,
            --bins analysis/similarity_bins/<phylum>_similarity_bins.tsv,
            --out analysis/kmer_baseline/<phylum>; C grid 0.001,0.01,0.1,1 as in the ESM-2 runs.
  The feature cache (analysis/kmer_baseline/features/) and fitted models (<out>/models/) are
  not committed; they are rebuilt on the first run.

Each fitted model is saved under <out>/models, so a power cut only loses the C value
that was running; rerunning the same command resumes.
"""
import argparse
import csv
import json
import os
import sys
import time
import warnings

import numpy as np
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LABELS = os.path.join(ROOT, "data", "proteins_ec3.tsv")
FEAT_DIR = os.path.join(ROOT, "analysis", "kmer_baseline", "features")
AA = "ACDEFGHIKLMNPQRSTVWY"
EXPECTED_PROTEINS = 146677
csv.field_size_limit(10**9)


def log(msg):
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


# ---------------------------------------------------------------- features
def kmer2(seq, lut):
    """400 overlapping 2-mer frequencies, counts / (L-1); pairs with a non-standard letter skipped."""
    s = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    idx = lut[s]                                   # 0..19 standard, -1 otherwise
    a, b = idx[:-1], idx[1:]
    ok = (a >= 0) & (b >= 0)
    counts = np.bincount(a[ok] * 20 + b[ok], minlength=400).astype(np.float32)
    return counts / max(len(s) - 1, 1)


def build_features():
    """Compute once for all proteins; cached as features/kmer2.npy + features/accessions.tsv."""
    npy = os.path.join(FEAT_DIR, "kmer2.npy")
    acc_path = os.path.join(FEAT_DIR, "accessions.tsv")
    info_path = os.path.join(FEAT_DIR, "feature_notes.json")
    if os.path.exists(npy) and os.path.exists(acc_path):
        return
    os.makedirs(FEAT_DIR, exist_ok=True)
    rows = read_tsv(LABELS)
    if len(rows) != EXPECTED_PROTEINS:
        sys.exit("proteins_ec3.tsv has %d rows, expected %d" % (len(rows), EXPECTED_PROTEINS))
    lut = np.full(256, -1, dtype=np.int64)
    for i, c in enumerate(AA):
        lut[ord(c)] = i
    X = np.empty((len(rows), 400), dtype=np.float32)
    nonstd = {}
    n_with_nonstd = 0
    for i, r in enumerate(rows):
        seq = r["sequence"].strip().upper()
        X[i] = kmer2(seq, lut)
        bad = set(seq) - set(AA)
        if bad:
            n_with_nonstd += 1
            for c in bad:
                nonstd[c] = nonstd.get(c, 0) + 1
    np.save(npy, X)
    with open(acc_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["accession", "row_index"])
        for i, r in enumerate(rows):
            w.writerow([r["accession"], i])
    with open(info_path, "w") as f:
        json.dump({"proteins": len(rows), "features": 400, "k": 2,
                   "normalisation": "count / (L - 1), L = full sequence length",
                   "non_standard_rule": "2-mers containing any letter outside the 20 standard amino acids are skipped",
                   "proteins_with_non_standard_letters": n_with_nonstd,
                   "non_standard_letters_count_of_proteins": nonstd}, f, indent=2)
    log("features built: %d x 400; %d proteins contain non-standard letters %s" % (len(rows), n_with_nonstd, nonstd))


def load_features():
    X = np.load(os.path.join(FEAT_DIR, "kmer2.npy"), mmap_mode="r")
    row_of = {r["accession"]: int(r["row_index"]) for r in read_tsv(os.path.join(FEAT_DIR, "accessions.tsv"))}
    return X, row_of


# ---------------------------------------------------------------- copied from run_logreg.py (unchanged)
def load_labels():
    info = {}
    for r in read_tsv(LABELS):
        info[r["accession"]] = {"ec3": r["ec3"], "phylum": r["phylum"], "lab_tested": r["lab_tested"]}
    return info


def normalize_group(value):
    v = value.strip().lower()
    if "unseen" in v:
        return "unseen_test"
    if "seen" in v:
        return "seen_test"
    if v.startswith("val"):
        return "validation"
    if "train" in v:
        return "train"
    return None


def load_split(path, phylum):
    rows = read_tsv(path)
    cols = list(rows[0].keys())
    acc_col = next((c for c in cols if "accession" in c.lower()), None)
    grp_col = next((c for c in cols if c.lower() in ("split", "set", "partition", "subset", "group", "role")), None)
    pyl_col = next((c for c in cols if c.lower() in ("held_out", "heldout", "held_out_phylum", "experiment")), None)
    if not acc_col or not grp_col:
        sys.exit("Could not find the accession / split columns. Columns are: %s" % cols)
    groups = {}
    for r in rows:
        if pyl_col and phylum and r[pyl_col] != phylum:
            continue
        g = normalize_group(r[grp_col])
        if g:
            groups[r[acc_col]] = g
    return groups


def load_bins(path):
    out = {}
    for r in read_tsv(path):
        acc = r.get("accession") or r.get(next(c for c in r if "accession" in c.lower()))
        out[acc] = {"bin": r.get("bin") or r.get("similarity_bin") or "",
                    "nearest": r.get("nearest_train_accession") or ""}
    return out


def softmax_predict(X, W, b, classes):
    scores = X @ W.T + b
    scores -= scores.max(axis=1, keepdims=True)
    p = np.exp(scores)
    p /= p.sum(axis=1, keepdims=True)
    best = p.argmax(axis=1)
    return classes[best], p[np.arange(len(best)), best]


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--split")
    ap.add_argument("--bins")
    ap.add_argument("--phylum", required=True)
    ap.add_argument("--out")
    ap.add_argument("--C", default="0.001,0.01,0.1,1")
    ap.add_argument("--max_iter", type=int, default=300)
    a = ap.parse_args()
    a.split = a.split or os.path.join(ROOT, "splits_ec3", a.phylum, "split_%s.tsv" % a.phylum)
    a.bins = a.bins or os.path.join(ROOT, "analysis", "similarity_bins", "%s_similarity_bins.tsv" % a.phylum)
    a.out = a.out or os.path.join(ROOT, "analysis", "kmer_baseline", a.phylum)
    os.makedirs(os.path.join(a.out, "models"), exist_ok=True)

    build_features()
    log("loading features and labels")
    feats, row_of = load_features()
    info = load_labels()
    groups = load_split(a.split, a.phylum)
    bins = load_bins(a.bins)

    missing = [x for x in groups if x not in row_of]
    log("split accessions: %d, missing from proteins_ec3.tsv: %d" % (len(groups), len(missing)))
    if missing:
        sys.exit("split has accessions not in proteins_ec3.tsv, e.g. %s" % missing[:5])

    by_group = {g: sorted(acc for acc, gg in groups.items() if gg == g and acc in row_of)
                for g in ("train", "validation", "seen_test", "unseen_test")}
    for g, accs in by_group.items():
        log("  %-12s %7d proteins" % (g, len(accs)))

    def matrix(accs):
        idx = np.array([row_of[x] for x in accs])
        X = np.asarray(feats[idx], dtype=np.float32)
        y = np.array([info[x]["ec3"] for x in accs])
        return X, y

    Xtr, ytr = matrix(by_group["train"])
    Xva, yva = matrix(by_group["validation"])

    mu = Xtr.mean(axis=0)
    sd = Xtr.std(axis=0) + 1e-6
    Xtr = (Xtr - mu) / sd
    Xva = (Xva - mu) / sd

    tried = []
    for C in [float(c) for c in a.C.split(",")]:
        path = os.path.join(a.out, "models", "C_%g.npz" % C)
        if os.path.exists(path):
            m = np.load(path, allow_pickle=True)
            W, b, classes, n_iter, secs = m["W"], m["b"], m["classes"], int(m["n_iter"]), float(m["secs"])
            conv = [str(x) for x in m["warnings"]] if "warnings" in m else []
            log("C=%g: already trained, loaded from disk" % C)
        else:
            log("C=%g: training on %d proteins, %d classes ..." % (C, len(ytr), len(set(ytr))))
            t = time.time()
            model = LogisticRegression(C=C, max_iter=a.max_iter, solver="lbfgs")
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                model.fit(Xtr, ytr)
            conv = sorted({str(w.message).split("\n")[0] for w in caught if issubclass(w.category, ConvergenceWarning)})
            secs = time.time() - t
            W, b, classes, n_iter = model.coef_, model.intercept_, model.classes_, int(model.n_iter_[0])
            np.savez(path, W=W, b=b, classes=classes, n_iter=n_iter, secs=secs, warnings=np.array(conv, dtype=str))
        pred, _ = softmax_predict(Xva, W, b, classes)
        acc = float((pred == yva).mean())
        tried.append({"C": C, "val_accuracy": acc, "iterations": n_iter, "seconds": round(secs),
                      "hit_max_iter": n_iter >= a.max_iter, "convergence_warnings": conv})
        log("C=%g: validation accuracy %.4f (%d iterations, %.0f s)%s" % (
            C, acc, n_iter, secs, "  CONVERGENCE WARNING" if conv else ""))

    best = max(tried, key=lambda r: r["val_accuracy"])
    log("best C = %g" % best["C"])
    m = np.load(os.path.join(a.out, "models", "C_%g.npz" % best["C"]), allow_pickle=True)
    W, b, classes = m["W"], m["b"], m["classes"]
    train_classes = set(classes.tolist())

    rows, summary = [], {}
    for g in ("seen_test", "unseen_test"):
        X, y = matrix(by_group[g])
        pred, prob = softmax_predict((X - mu) / sd, W, b, classes)
        ok = pred == y
        summary[g] = {"proteins": len(y), "accuracy": round(float(ok.mean()), 4),
                      "class_not_in_train": int(sum(c not in train_classes for c in y))}
        for acc, t, p, pr, k in zip(by_group[g], y, pred, prob, ok):
            bn = bins.get(acc, {})
            near = bn.get("nearest", "")
            b1 = info[near]["ec3"] if near in info else ""
            rows.append([acc, g, info[acc]["lab_tested"], bn.get("bin", ""), t, p, "%.4f" % pr, int(k),
                         int(t in train_classes), near, b1, int(b1 == t) if b1 else ""])
        b1 = [r[11] for r in rows if r[1] == g]
        summary[g]["baseline1_accuracy_all"] = round(sum(x == 1 for x in b1) / len(b1), 4)
        log("%s: accuracy %.4f on %d proteins" % (g, summary[g]["accuracy"], len(y)))

    with open(os.path.join(a.out, "predictions.tsv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["accession", "group", "lab_tested", "bin", "true_ec3", "pred_ec3", "prob_pred", "correct",
                    "class_in_train", "nearest_train_accession", "baseline1_pred", "baseline1_correct"])
        w.writerows(rows)
    with open(os.path.join(a.out, "summary.json"), "w") as f:
        json.dump({"mode": "real", "features": "kmer2 (400 overlapping 2-mer frequencies)", "phylum": a.phylum,
                   "sizes": {g: len(v) for g, v in by_group.items()},
                   "tuning": tried, "best_C": best["C"], "test": summary}, f, indent=2)
    log("done: predictions.tsv and summary.json written to %s" % a.out)


if __name__ == "__main__":
    main()
