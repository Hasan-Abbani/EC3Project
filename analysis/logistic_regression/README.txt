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
Paths below are relative to analysis/logistic_regression/.

README.txt                 this file
results_overview.tsv       one row per held-out phylum (sizes, best C, accuracies)
../../code/run_logreg.py   the script in the repository's top-level code/ directory
<phylum>/predictions.tsv   one row per test protein (seen and unseen test)
<phylum>/summary.json      sizes, validation accuracy for every C, best C, test accuracy

Generated locally and excluded from Git:
<phylum>/models/C_<C>.npz  the fitted model for each C value tried

The model files described below belong to the local training archive. Cloning
the repository does not download those fitted arrays; retraining regenerates them.

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
The normal default C grid now matches the four values used in the recorded
runs. The separate --practice mode retains its single-value default.

From the repository root, using a new local output directory:
  python code/run_logreg.py --split splits_ec3/<phylum>/split_<phylum>.tsv --bins analysis/similarity_bins/<phylum>_similarity_bins.tsv --phylum <phylum> --C 0.001,0.01,0.1,1 --out .cache/reproduction/<phylum>
