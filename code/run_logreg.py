"""Phase E: logistic regression on frozen ESM-2 embeddings (predict EC-3 class).

Two ways to run it:

  REAL RUN (Hasan's split + bins):
    python -I run_logreg.py --split <split.tsv> --bins <phylum_similarity_bins.tsv> --phylum Actinomycetota --out <folder>

  PRACTICE RUN (our own quick split, small training sample, one C value):
    python -I run_logreg.py --practice --out <folder>

What it does, in order:
  1. Load the 146,677 x 1280 embedding table and its row order (accessions.tsv).
  2. Load the true EC-3 label of every protein (proteins_ec3.tsv).
  3. Load the split: which proteins are train / validation / seen test / unseen test.
  4. Standardize the 1,280 numbers (using train only).
  5. For each C value: train on train, score on validation. Each fitted model is
     saved, so a power cut only loses the C value that was running.
  6. Keep the C with the best validation accuracy.
  7. Predict the seen and unseen test proteins; save one row per protein.
"""
import argparse
import csv
import json
import os
import random
import sys
import time

import numpy as np
from sklearn.linear_model import LogisticRegression

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))  # .../CMPS 355
EMB_DIR = os.path.join(PROJECT, "Me", "Other Claude", "Claude", "EC3_ESM2_Results", "esm2_embeddings_cuda")
LABELS = os.path.join(HERE, "..", "Data", "Cleaning - EC Level 3", "2 Cleaned Data (send after yes)",
                      "data", "proteins_ec3.tsv")
SEED = 355
csv.field_size_limit(10**9)


def log(msg):
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


# ---------------------------------------------------------------- 1-2. data
def load_embeddings():
    # mmap: the 750 MB file is not read into memory until we ask for rows.
    emb = np.load(os.path.join(EMB_DIR, "embeddings.npy"), mmap_mode="r")
    row_of = {r["accession"]: int(r["row_index"]) for r in read_tsv(os.path.join(EMB_DIR, "accessions.tsv"))}
    return emb, row_of


def load_labels():
    info = {}
    for r in read_tsv(LABELS):
        info[r["accession"]] = {"ec3": r["ec3"], "phylum": r["phylum"], "lab_tested": r["lab_tested"]}
    return info


# ---------------------------------------------------------------- 3. split
def normalize_group(value):
    v = value.strip().lower()
    if "unseen" in v:                       # check "unseen" before "seen"
        return "unseen_test"
    if "seen" in v:
        return "seen_test"
    if v.startswith("val"):
        return "validation"
    if "train" in v:
        return "train"
    return None


def load_split(path, phylum):
    """Read Hasan's split table. Column names are found, not assumed."""
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


def practice_split(info):
    """Our own split, like Hasan's: hold out Actinomycetota entirely (unseen test);
    from the other phyla take 10% seen test, 10% validation, rest train."""
    rng = random.Random(SEED)
    groups = {}
    others = []
    for acc, d in info.items():
        if d["phylum"] == "Actinomycetota":
            groups[acc] = "unseen_test"
        else:
            others.append(acc)
    others.sort()
    rng.shuffle(others)
    n = len(others) // 10
    for a in others[:n]:
        groups[a] = "seen_test"
    for a in others[n:2 * n]:
        groups[a] = "validation"
    for a in others[2 * n:]:
        groups[a] = "train"
    return groups


def load_bins(path):
    out = {}
    for r in read_tsv(path):
        acc = r.get("accession") or r.get(next(c for c in r if "accession" in c.lower()))
        out[acc] = {"bin": r.get("bin") or r.get("similarity_bin") or "",
                    "nearest": r.get("nearest_train_accession") or ""}
    return out


# ---------------------------------------------------------------- 4-7. model
def softmax_predict(X, W, b, classes):
    """This is the whole model: a score per class = weights . x + bias,
    then softmax turns the scores into probabilities that add up to 1."""
    scores = X @ W.T + b                                  # (n proteins, n classes)
    scores -= scores.max(axis=1, keepdims=True)           # for numerical safety
    p = np.exp(scores)
    p /= p.sum(axis=1, keepdims=True)
    best = p.argmax(axis=1)
    return classes[best], p[np.arange(len(best)), best]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--split")
    ap.add_argument("--bins")
    ap.add_argument("--phylum", default="Actinomycetota")
    ap.add_argument("--out", required=True)
    ap.add_argument("--practice", action="store_true")
    ap.add_argument("--C", default=None, help="C values to try, comma-separated (default 0.01,0.1,1)")
    ap.add_argument("--train_sample", type=int, default=0, help="use only this many training proteins (0 = all)")
    ap.add_argument("--max_iter", type=int, default=300)
    a = ap.parse_args()
    if a.practice:
        a.C = a.C or "1"
        a.train_sample = a.train_sample or 20000
    a.C = a.C or "0.01,0.1,1"
    os.makedirs(os.path.join(a.out, "models"), exist_ok=True)

    log("loading embeddings and labels")
    emb, row_of = load_embeddings()
    info = load_labels()
    groups = practice_split(info) if a.practice else load_split(a.split, a.phylum)
    bins = load_bins(a.bins) if a.bins else {}

    by_group = {g: sorted(acc for acc, gg in groups.items() if gg == g and acc in row_of)
                for g in ("train", "validation", "seen_test", "unseen_test")}
    if a.train_sample and a.train_sample < len(by_group["train"]):
        by_group["train"] = sorted(random.Random(SEED).sample(by_group["train"], a.train_sample))
    for g, accs in by_group.items():
        log("  %-12s %7d proteins" % (g, len(accs)))

    def matrix(accs):
        idx = np.array([row_of[x] for x in accs])
        order = np.argsort(idx)                    # read the file in row order (faster)
        X = np.empty((len(idx), emb.shape[1]), dtype=np.float32)
        X[order] = emb[idx[order]]
        y = np.array([info[x]["ec3"] for x in accs])
        return X, y

    Xtr, ytr = matrix(by_group["train"])
    Xva, yva = matrix(by_group["validation"])

    # 4. standardize: each of the 1,280 numbers gets mean 0 and spread 1, using TRAIN only
    mu = Xtr.mean(axis=0)
    sd = Xtr.std(axis=0) + 1e-6
    Xtr = (Xtr - mu) / sd
    Xva = (Xva - mu) / sd

    # 5-6. try each C, score on validation, keep the best
    tried = []
    for C in [float(c) for c in a.C.split(",")]:
        path = os.path.join(a.out, "models", "C_%g.npz" % C)
        if os.path.exists(path):
            m = np.load(path, allow_pickle=True)
            W, b, classes, n_iter, secs = m["W"], m["b"], m["classes"], int(m["n_iter"]), float(m["secs"])
            log("C=%g: already trained, loaded from disk" % C)
        else:
            log("C=%g: training on %d proteins, %d classes ..." % (C, len(ytr), len(set(ytr))))
            t = time.time()
            model = LogisticRegression(C=C, max_iter=a.max_iter, solver="lbfgs")
            model.fit(Xtr, ytr)
            secs = time.time() - t
            W, b, classes, n_iter = model.coef_, model.intercept_, model.classes_, int(model.n_iter_[0])
            np.savez(path, W=W, b=b, classes=classes, n_iter=n_iter, secs=secs)
        pred, _ = softmax_predict(Xva, W, b, classes)
        acc = float((pred == yva).mean())
        tried.append({"C": C, "val_accuracy": acc, "iterations": n_iter, "seconds": round(secs)})
        log("C=%g: validation accuracy %.4f (%d iterations, %.0f s)" % (C, acc, n_iter, secs))

    best = max(tried, key=lambda r: r["val_accuracy"])
    log("best C = %g" % best["C"])
    m = np.load(os.path.join(a.out, "models", "C_%g.npz" % best["C"]), allow_pickle=True)
    W, b, classes = m["W"], m["b"], m["classes"]
    train_classes = set(classes.tolist())

    # 7. predict both test sets, one row per protein
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
        if bins:
            b1 = [r[11] for r in rows if r[1] == g]
            summary[g]["baseline1_accuracy_all"] = round(sum(x == 1 for x in b1) / len(b1), 4)
        log("%s: accuracy %.4f on %d proteins" % (g, summary[g]["accuracy"], len(y)))

    with open(os.path.join(a.out, "predictions.tsv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["accession", "group", "lab_tested", "bin", "true_ec3", "pred_ec3", "prob_pred", "correct",
                    "class_in_train", "nearest_train_accession", "baseline1_pred", "baseline1_correct"])
        w.writerows(rows)
    with open(os.path.join(a.out, "summary.json"), "w") as f:
        json.dump({"mode": "practice" if a.practice else "real", "phylum": a.phylum,
                   "sizes": {g: len(v) for g, v in by_group.items()},
                   "tuning": tried, "best_C": best["C"], "test": summary}, f, indent=2)
    log("done: predictions.tsv and summary.json written to %s" % a.out)


if __name__ == "__main__":
    main()
