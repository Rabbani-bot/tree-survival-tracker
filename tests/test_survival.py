import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from survival import kaplan_meier, logrank, generate  # noqa: E402


class TestSurvival(unittest.TestCase):
    def test_km_textbook(self):
        # times 1..5, deaths at 1,3,5; censored at 2,4 -> S = 0.8, 0.533, 0
        t = np.array([1, 2, 3, 4, 5]); e = np.array([1, 0, 1, 0, 1])
        km = kaplan_meier(t, e)
        self.assertAlmostEqual(km["survival"].iloc[1], 0.8)
        self.assertAlmostEqual(km["survival"].iloc[2], 0.8 * (2 / 3))
        self.assertAlmostEqual(km["survival"].iloc[3], 0.0)

    def test_km_monotone(self):
        df = generate()
        km = kaplan_meier(df["time_months"].to_numpy(), df["died"].to_numpy())
        self.assertTrue((np.diff(km["survival"]) <= 0).all())

    def test_logrank_identical_groups_not_significant(self):
        df = generate()
        df["random_group"] = np.where(np.arange(len(df)) % 2 == 0, "x", "y")
        _, p = logrank(df, "random_group")
        self.assertGreater(p, 0.01)


if __name__ == "__main__":
    unittest.main()
