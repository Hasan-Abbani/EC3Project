"""
clean_ec3.py - Clean the UniProt download for our project (predicting 3 EC levels).

WHAT IT DOES
    Reads the original UniProt file (never changed) and applies the nine agreed
    cleaning rules, one after the other. Every protein that is removed is written
    to a log file together with the rule that removed it, so nothing disappears
    silently.

HOW TO RUN (from any folder)
    python -I clean_ec3.py  <raw .tsv.gz>  <enzclass.txt>  <output folder>

    raw .tsv.gz   = Data/raw/uniprot_reviewed_bacteria_ec_2026_03.tsv.gz
    enzclass.txt  = Data/raw/enzclass.txt  (EC class names, ExPASy ENZYME release 02-Sep-2026)
    output folder = where the cleaned files go (it is created if missing)

WHAT IT WRITES (into the output folder)
    proteins_ec3.tsv      one row per kept protein; the main file (tab separated)
    proteins_ec3.fasta    the same proteins as FASTA (input for ESM-2 and MMseqs2)
    removed_proteins.tsv  every removed protein and the rule that removed it
    cleaning_steps.csv    proteins left after each rule (the cleaning funnel)
    summary.json          counts used by the Excel file and the Data Explorer

THE RULES, IN ORDER (numbers match "Cleaning Decisions - To Confirm.pdf")
    1  drop fragments                       (partial sequences)
    2  drop proteins with no phylum         (cannot be called seen / unseen)
    3  keep proteins whose EC numbers are ALL known to at least 3 levels
       (3.1.-.- is too vague; one vague number drops the whole protein)
    4  keep proteins whose EC numbers all give the SAME 3-level label
       (3.2.2.6; 3.2.2.7 -> 3.2.2, kept.  3.5.4.19; 3.6.1.31 -> two labels, dropped)
    5  drop sequences labelled with different EC numbers in different entries
    6  drop sequences longer than 1,022 amino acids  (chosen analysis length cap)
    7  drop sequences with X, B or Z  (unknown amino acid)
    8  keep one copy of each identical sequence per phylum (prefer a lab-tested copy)
    9  keep every evidence level, and add a yes/no "lab_tested" column

    Not applied here, on purpose: the "at least 10 proteins per EC class" rule.
    It must hold in the TRAINING set, so it belongs to the splitting step. The
    column ec3_class_size is provided so it can be applied in one line later.

Only the Python standard library is used. Same input -> same output, every time.
"""
import csv
import gzip
import json
import os
import re
import sys
from collections import Counter, defaultdict

# --------------------------------------------------------------------------
# Settings
# --------------------------------------------------------------------------
MAX_LENGTH = 1022            # Chosen analysis cap for embedding without truncation
UNKNOWN_LETTERS = set("XBZ")  # X = unknown, B = D or N, Z = E or Q
LAB_CODE = "ECO:0000269"     # evidence code "experimental evidence used in manual assertion"
MIN_CLASS = 10               # used only for the summary numbers, not for removing proteins
TOP_PHYLA = 8                # how many phyla get their own entry in the summary charts

raw_path, enzclass_path, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(out_dir, exist_ok=True)
csv.field_size_limit(10**9)   # some UniProt fields are very long


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------
def rank_from_lineage(lineage, rank):
    """Take one rank out of UniProt's lineage text.

    The lineage looks like "Bacteria (domain), Pseudomonadota (phylum), ...",
    so we look for the name written just before "(phylum)" (or "(class)", ...).
    Returns "" when the rank is missing.
    """
    m = re.search(r"(?:^|, )([^,]+?) \(%s\)" % rank, lineage)
    return m.group(1) if m else ""


def known_levels(ec):
    """How many levels of an EC number are known, counting from the left.

    "1.1.1.1"  -> 4      "3.2.2.-"  -> 3
    "3.5.1.n3" -> 3      "3.1.-.-"  -> 2   (n3 = provisional number, not a real level)
    """
    depth = 0
    for part in ec.split("."):
        if not part.isdigit():
            break
        depth += 1
    return depth


def read_ec_names(path):
    """Read ExPASy's enzclass.txt into {"1.1.1": "With NAD(+) or NADP(+) as acceptor", ...}.

    Lines look like "1. 1. 1.-    With NAD(+) or NADP(+) as acceptor."
    The code part has spaces for alignment; we remove them and drop the "-" levels.
    """
    names = {}
    for line in open(path, encoding="utf-8", errors="replace"):
        m = re.match(r"^(\d)\.\s*([\d-]+)\.\s*([\d-]+)\.-\s+(.+?)\.?\s*$", line)
        if m:
            code = ".".join(p for p in m.groups()[:3] if p != "-")
            names[code] = m.group(4)
    return names


# --------------------------------------------------------------------------
# Read the original download
# --------------------------------------------------------------------------
with gzip.open(raw_path, "rt", encoding="utf-8", newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE))

# Work out the extra fields once, so each rule below is a short test.
for r in rows:
    r["ecs"] = [e.strip() for e in r["EC number"].split(";") if e.strip()]
    r["phylum"] = rank_from_lineage(r["Taxonomic lineage"], "phylum")
    r["class"] = rank_from_lineage(r["Taxonomic lineage"], "class")
    r["genus"] = rank_from_lineage(r["Taxonomic lineage"], "genus")
    r["length"] = int(r["Length"])
    r["lab"] = LAB_CODE in r["Catalytic activity"]
    r["copies"] = 1   # how many identical rows this row stands for (filled in by rule 8)
downloaded_per_phylum = Counter(r["phylum"] for r in rows)  # kept for the "fair to every phylum" table

steps = [("0 downloaded", len(rows), 0)]   # the funnel: (rule, proteins left, removed)
removed = []                                # (accession, phylum, EC, rule) for the log


def apply_rule(name, keep_if):
    """Keep the rows for which keep_if(row) is True; log the others under this rule's name."""
    global rows
    kept = []
    for r in rows:
        if keep_if(r):
            kept.append(r)
        else:
            removed.append((r["Entry"], r["phylum"], r["EC number"], name))
    steps.append((name, len(kept), len(rows) - len(kept)))
    rows = kept


# --------------------------------------------------------------------------
# The rules
# --------------------------------------------------------------------------
# Rule 1. A fragment is only a piece of a protein; its embedding would describe a piece.
apply_rule("1 fragment", lambda r: not r["Fragment"].strip())

# Rule 2. Our question is "seen phylum or unseen phylum"; without a phylum a protein fits neither.
apply_rule("2 no phylum", lambda r: r["phylum"] != "")

# Rule 3. We predict 3 levels, so every EC number must know at least 3 levels.
#         "3.2.2.-" passes (its 3-level label is 3.2.2); "3.1.-.-" does not.
#         If a protein has several EC numbers and ANY of them is vague, the whole
#         protein goes: "2.7.7.-; 3.1.-.-" may be two different jobs, and dropping
#         only the vague number would hide the second one.
apply_rule("3 EC known to fewer than 3 levels", lambda r: min(known_levels(e) for e in r["ecs"]) >= 3)

# Every remaining EC number can now be cut to its 3-level label.
for r in rows:
    r["labels3"] = {".".join(e.split(".")[:3]) for e in r["ecs"]}

# Rule 4. Multinomial logistic regression gives one class per protein. If the EC
#         numbers give two different 3-level labels (3.5.4.19; 3.6.1.31), the
#         protein has two right answers, so we leave it out. If they give the same
#         one (3.2.2.6; 3.2.2.7 -> 3.2.2), there is one right answer: kept.
apply_rule("4 several 3-level labels", lambda r: len(r["labels3"]) == 1)

# From here on every protein has exactly one 3-level label: split it into its levels.
for r in rows:
    label = next(iter(r["labels3"]))
    r["ec"] = "; ".join(sorted(r["ecs"]))   # all EC numbers, in a fixed order
    for lvl in (1, 2, 3):
        r["ec%d" % lvl] = ".".join(label.split(".")[:lvl])

# Rule 5. If the same sequence carries different EC numbers in different entries,
#         we cannot tell which is right. Judged on the full EC numbers (all of
#         them, for proteins with several), as in the cleaning report.
labels_per_sequence = defaultdict(set)
for r in rows:
    labels_per_sequence[r["Sequence"]].add(r["ec"])
apply_rule("5 contradicting labels", lambda r: len(labels_per_sequence[r["Sequence"]]) == 1)

# Rule 6. Apply the chosen 1,022-residue analysis cap so retained proteins are
#         embedded without truncation; this is not a universal ESM-2 limit.
apply_rule("6 longer than 1,022", lambda r: r["length"] <= MAX_LENGTH)

# Rule 7. X, B, Z mean "unknown or ambiguous amino acid". U and O are real, rare ones: kept.
apply_rule("7 unknown amino acid (X/B/Z)", lambda r: not (set(r["Sequence"]) & UNKNOWN_LETTERS))

# Rule 8. The same sequence is stored once per strain. Identical copies could end up
#         in both training and test, which would inflate accuracy. Keep one copy per
#         (sequence, phylum). Sorting first means the kept copy is the lab-tested one
#         when there is one, and otherwise the one with the smallest accession (stable).
rows.sort(key=lambda r: (not r["lab"], r["Entry"]))
first_copy = {}
kept = []
for r in rows:
    key = (r["Sequence"], r["phylum"])
    if key in first_copy:
        first_copy[key]["copies"] += 1
        removed.append((r["Entry"], r["phylum"], r["EC number"], "8 extra copy of " + first_copy[key]["Entry"]))
    else:
        first_copy[key] = r
        kept.append(r)
steps.append(("8 extra identical copy", len(kept), len(rows) - len(kept)))
rows = sorted(kept, key=lambda r: r["Entry"])

# Rule 9 removes nothing: the "lab_tested" column is written below.

# --------------------------------------------------------------------------
# Write the cleaned files
# --------------------------------------------------------------------------
ec_names = read_ec_names(enzclass_path)
class_size = Counter(r["ec3"] for r in rows)

COLUMNS = [  # (column name in the output, how to get it from a row)
    ("accession", lambda r: r["Entry"]),
    ("entry_name", lambda r: r["Entry Name"]),
    ("protein_name", lambda r: r["Protein names"]),
    ("organism", lambda r: r["Organism"]),
    ("taxon_id", lambda r: r["Organism (ID)"]),
    ("phylum", lambda r: r["phylum"]),
    ("class", lambda r: r["class"]),
    ("genus", lambda r: r["genus"]),
    ("ec_full", lambda r: r["ec"]),
    ("ec_count", lambda r: len(r["ecs"])),
    ("ec1", lambda r: r["ec1"]),
    ("ec2", lambda r: r["ec2"]),
    ("ec3", lambda r: r["ec3"]),
    ("ec3_name", lambda r: ec_names.get(r["ec3"], "")),
    ("ec3_class_size", lambda r: class_size[r["ec3"]]),
    ("lab_tested", lambda r: "yes" if r["lab"] else "no"),
    ("identical_copies_merged", lambda r: r["copies"]),
    ("length", lambda r: r["length"]),
    ("sequence", lambda r: r["Sequence"]),
]

with open(os.path.join(out_dir, "proteins_ec3.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow([name for name, _ in COLUMNS])
    for r in rows:
        w.writerow([get(r) for _, get in COLUMNS])

with open(os.path.join(out_dir, "proteins_ec3.fasta"), "w", encoding="utf-8", newline="\n") as f:
    for r in rows:
        f.write(">%s %s %s\n%s\n" % (r["Entry"], r["phylum"], r["ec3"], r["Sequence"]))

with open(os.path.join(out_dir, "removed_proteins.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["accession", "phylum", "ec_as_downloaded", "removed_by_rule"])
    w.writerows(sorted(removed))

with open(os.path.join(out_dir, "cleaning_steps.csv"), "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["rule", "proteins_left", "removed"])
    w.writerows(steps)

# --------------------------------------------------------------------------
# Summary numbers (read by make_excel.py and make_explorer.py)
# --------------------------------------------------------------------------
EC1_NAMES = {"1": "Oxidoreductases", "2": "Transferases", "3": "Hydrolases", "4": "Lyases",
             "5": "Isomerases", "6": "Ligases", "7": "Translocases"}

phylum_count = Counter(r["phylum"] for r in rows)
top = [p for p, _ in phylum_count.most_common(TOP_PHYLA)]


def group(p):
    """Small phyla are grouped as "Other phyla" in the charts."""
    return p if p in top else "Other phyla"


groups = top + ["Other phyla"]
BIN = 50  # length histogram bin width


def describe(sub):
    """Counts for one set of proteins (all of them, or one phylum)."""
    return {
        "proteins": len(sub),
        "lab_tested": sum(1 for r in sub if r["lab"]),
        "several_ec": sum(1 for r in sub if len(r["ecs"]) > 1),
        "ec3_classes": len({r["ec3"] for r in sub}),
        "ec1": {k: sum(1 for r in sub if r["ec1"] == k) for k in EC1_NAMES},
        "length_bins": [sum(1 for r in sub if b <= r["length"] < b + BIN) for b in range(0, MAX_LENGTH + 1, BIN)],
        "median_length": sorted(r["length"] for r in sub)[len(sub) // 2] if sub else 0,
    }


summary = {
    "source": "UniProtKB/Swiss-Prot, reviewed bacterial proteins with an EC number, release 2026_03",
    "steps": steps,
    "groups": groups,
    "ec1_names": EC1_NAMES,
    "length_bin": BIN,
    "all": describe(rows),
    "by_group": {g: describe([r for r in rows if group(r["phylum"]) == g]) for g in groups},
    "phyla": [{"phylum": p, "proteins": n, "lab_tested": sum(1 for r in rows if r["phylum"] == p and r["lab"]),
               "ec3_classes": len({r["ec3"] for r in rows if r["phylum"] == p})}
              for p, n in phylum_count.most_common()],
}

# One line per EC class (3 levels): size, name, how many per phylum group.
summary["ec3"] = []
for c, n in class_size.most_common():
    sub = [r for r in rows if r["ec3"] == c]
    summary["ec3"].append({
        "ec3": c, "name": ec_names.get(c, ""), "ec1": c.split(".")[0], "proteins": n,
        "lab_tested": sum(1 for r in sub if r["lab"]),
        "by_group": {g: sum(1 for r in sub if group(r["phylum"]) == g) for g in groups},
    })
summary["min_class"] = {k: {"classes": sum(1 for v in class_size.values() if v >= k),
                            "proteins": sum(v for v in class_size.values() if v >= k)}
                        for k in (1, 5, 10, 20, 50, 100)}

# Removed proteins per rule and phylum group: does a rule hit one phylum harder?
hit = defaultdict(Counter)
for acc, phylum, ec, rule in removed:
    hit[rule.split(" ")[0]][group(phylum) if phylum else "no phylum"] += 1
summary["removed_by_rule_and_group"] = {k: dict(v) for k, v in sorted(hit.items())}
before = Counter()
for p, n in downloaded_per_phylum.items():
    before[group(p) if p else "no phylum"] += n
summary["downloaded_by_group"] = dict(before)

# If one phylum is held out: how many of its proteins have a class with at least
# MIN_CLASS proteins in the other phyla (so the classifier could have learned it)?
summary["holdout"] = []
for p in top:
    train = Counter(r["ec3"] for r in rows if r["phylum"] != p)
    learnable = {c for c, m in train.items() if m >= MIN_CLASS}
    own = [r for r in rows if r["phylum"] == p]
    summary["holdout"].append({
        "phylum": p, "proteins": len(own), "training_left": len(rows) - len(own),
        "usable": sum(1 for r in own if r["ec3"] in learnable),
        "own_classes": len({r["ec3"] for r in own}),
        "classes_not_in_other_phyla": len({r["ec3"] for r in own} - set(train)),
    })

with open(os.path.join(out_dir, "summary.json"), "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=1)

# --------------------------------------------------------------------------
# Print the funnel so the run can be checked at a glance
# --------------------------------------------------------------------------
for name, left, gone in steps:
    print("%-36s left %7d   removed %6d" % (name, left, gone))
print("EC classes (3 levels):", len(class_size), "| lab-tested:", summary["all"]["lab_tested"])
print("classes with >= %d proteins: %d holding %d proteins"
      % (MIN_CLASS, summary["min_class"][MIN_CLASS]["classes"], summary["min_class"][MIN_CLASS]["proteins"]))
