LOGISTIC REGRESSION ON FROZEN ESM-2 EMBEDDINGS
==============================================

What was done
-------------
For each held-out phylum (your splits and similarity bins), a multinomial
logistic regression was trained on the 1,280-number ESM-2 embeddings to
predict the EC-3 class. Embeddings were standardized using the training
proteins only. C was tried at 0.001, 0.01, 0.1 and 1 and picked by
validation accuracy; the test sets were not used for the choice.
Each run takes about 4 minutes.

Baseline 1 = copy the EC-3 of nearest_train_accession; proteins with no
hit ("no_hit") are counted as wrong.

Folder layout
-------------
README.txt                 this file
results_overview.tsv       one row per held-out phylum (sizes, best C, accuracies)
code/run_logreg.py         the script that produced everything below
<phylum>/predictions.tsv   one row per test protein (seen and unseen test)
<phylum>/summary.json      sizes, validation accuracy for every C, best C, test accuracy
<phylum>/models/C_<C>.npz  the fitted model for each C value tried

predictions.tsv columns
-----------------------
accession, group (seen_test / unseen_test), lab_tested (yes/no), bin
(similarity bin from your bins file, or no_hit), true_ec3, pred_ec3,
prob_pred (model's probability for its prediction), correct (1/0),
class_in_train (1 if the true class appears in train; always 1 here),
nearest_train_accession, baseline1_pred, baseline1_correct

Model files (.npz, load with numpy.load)
----------------------------------------
W (classes x 1280 weights), b (one intercept per class), classes (EC-3 label
of each row), n_iter, secs. Inputs must be standardized first with the mean
and standard deviation of the TRAINING embeddings of that split (not stored;
the script recomputes them, and they are deterministic).

Re-running
----------
Needs numpy and scikit-learn (1.9.1 used). The paths to embeddings.npy /
accessions.tsv (EMB_DIR) and proteins_ec3.tsv (LABELS) are set at the top
of run_logreg.py for our folder layout; edit those two lines to your paths.
  python run_logreg.py --split split_<phylum>.tsv --bins <phylum>_similarity_bins.tsv --phylum <phylum> --C 0.001,0.01,0.1,1 --out <folder>
