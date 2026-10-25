"""Regression checks for identity selection and explicit missing-hit handling."""
import tempfile
import unittest
from pathlib import Path

from nearest_train_identity import reduce_hits


class IdentityReductionTests(unittest.TestCase):
    def test_identity_maximum_can_differ_from_first_and_best_score(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "hits.tsv"
            path.write_text(
                "q1\tt1\t0.4\t40\t100\t0.9\t0.9\t1e-40\t200\n"
                "q1\tt2\t0.6\t60\t100\t0.9\t0.9\t1e-20\t100\n",
                encoding="utf-8")
            best, counts, first_below = reduce_hits(path, ["q1", "q2"], {"t1", "t2"})
            self.assertEqual(best.iloc[0].nearest_train_accession, "t2")
            self.assertEqual(best.iloc[0].max_train_identity, 0.6)
            self.assertEqual(counts["q2"], 0)
            self.assertNotIn("q2", set(best.accession))
            self.assertEqual(first_below, 1)

    def test_target_must_belong_to_final_training_set(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "hits.tsv"
            path.write_text("q1\tforeign\t0.4\t40\t100\t0.9\t0.9\t1e-20\t100\n",
                            encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "final training"):
                reduce_hits(path, ["q1"], {"t1"})

    def test_empty_hit_table_does_not_invent_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "hits.tsv"
            path.write_text("", encoding="utf-8")
            best, counts, first_below = reduce_hits(path, ["q1", "q2"], {"t1"})
            self.assertTrue(best.empty)
            self.assertEqual(counts.to_dict(), {"q1": 0, "q2": 0})
            self.assertEqual(first_below, 0)

    def test_reported_identity_is_maximized_without_extra_precision(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "hits.tsv"
            path.write_text(
                "q1\tt1\t0.6\t60\t100\t0.9\t0.9\t1e-40\t200\n"
                "q1\tt2\t0.6\t301\t501\t0.9\t0.9\t1e-20\t100\n",
                encoding="utf-8")
            best, _, _ = reduce_hits(path, ["q1"], {"t1", "t2"})
            self.assertEqual(best.iloc[0].nearest_train_accession, "t1")
            self.assertEqual(best.iloc[0].max_train_identity, 0.6)


if __name__ == "__main__":
    unittest.main()
