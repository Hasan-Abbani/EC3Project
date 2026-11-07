import unittest

import numpy as np
import pandas as pd

from evaluate_similarity_bins import (
    _class_distribution_tv,
    _encoded_metrics,
    _reweighted_accuracy,
    bootstrap_metrics,
)


class BinwiseEvaluationTests(unittest.TestCase):
    def test_metrics_for_known_predictions(self):
        result = _encoded_metrics(["a", "a", "b", "b"], ["a", "b", "b", "a"])
        self.assertAlmostEqual(result["accuracy"], 0.5)
        self.assertAlmostEqual(result["macro_f1"], 0.5)
        self.assertAlmostEqual(result["balanced_accuracy"], 0.5)

    def test_shared_class_tv_distance(self):
        seen = pd.DataFrame({"true_ec3": ["a", "a", "a", "b"]})
        unseen = pd.DataFrame({"true_ec3": ["a", "b", "b", "b"]})
        rows, tv, shared_count = _class_distribution_tv(
            "phylum", "30-40%", seen, unseen
        )
        self.assertAlmostEqual(tv, 0.5)
        self.assertEqual(shared_count, 2)
        self.assertEqual(len(rows), 2)

    def test_reweighting_matches_unseen_common_class_mix(self):
        seen = pd.DataFrame({
            "true_ec3": ["a", "a", "a", "b"],
            "pred_ec3": ["a", "a", "b", "b"],
        })
        unseen = pd.DataFrame({
            "true_ec3": ["a", "b", "b", "b"],
            "pred_ec3": ["a", "b", "b", "a"],
        })
        result = _reweighted_accuracy(seen, unseen)
        self.assertAlmostEqual(result["class_mix_adjusted_seen_accuracy"], 11 / 12)
        self.assertAlmostEqual(result["unseen_common_support_accuracy"], 0.75)
        self.assertAlmostEqual(result["adjusted_gap_unseen_minus_seen"], -1 / 6)
        self.assertAlmostEqual(result["weight_sum"], 1)

    def test_bootstrap_is_reproducible(self):
        frame = pd.DataFrame({
            "true_ec3": ["a", "a", "b", "b"],
            "pred_ec3": ["a", "b", "b", "a"],
        })
        left = bootstrap_metrics(frame, 30, np.random.default_rng(355))
        right = bootstrap_metrics(frame, 30, np.random.default_rng(355))
        for metric in left:
            np.testing.assert_array_equal(left[metric], right[metric])


if __name__ == "__main__":
    unittest.main()
