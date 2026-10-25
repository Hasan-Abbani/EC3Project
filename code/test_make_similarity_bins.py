"""Check scientific boundary conventions and no-hit separation."""
import unittest
import pandas as pd
from make_similarity_bins import assign


class SimilarityBinTests(unittest.TestCase):
    def test_exact_boundaries_and_perfect_identity(self):
        values = pd.Series(["0", "0.188", "0.299", "0.3", "0.4", "0.5", "0.6", "0.7", "1", ""])
        status = pd.Series(["hit"] * 9 + ["no_qualifying_hit"])
        actual = assign(values, status, [0,30,40,50,60,70,100]).tolist()
        self.assertEqual(actual, ["<30%","<30%","<30%","30-40%","40-50%",
                                  "50-60%","60-70%","70-100%","70-100%","no_hit"])

    def test_missing_hit_is_not_zero(self):
        for value, status in [("0", "no_qualifying_hit"), ("", "hit")]:
            with self.assertRaises(ValueError):
                assign(pd.Series([value]),pd.Series([status]),[0,100])

    def test_invalid_identities_rejected(self):
        for value in ["-0.001", "1.001", "NaN", "Infinity"]:
            with self.assertRaises(ValueError):
                assign(pd.Series([value]),pd.Series(["hit"]),[0,100])


if __name__ == "__main__":
    unittest.main()
