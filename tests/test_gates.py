import pathlib
import sys
import unittest
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).parents[1] / "scripts"))

from part6_smoke_gate import judge_smoke


class GateTests(unittest.TestCase):
    def test_smoke_requires_all_metrics(self):
        table = pd.DataFrame({"perturbation": ["KO"], "max_abs_cosine_diff": [0.0]})
        with self.assertRaises(SystemExit):
            judge_smoke(table, perturbations=["KO"])


    def test_smoke_accepts_finite_metrics(self):
        table = pd.DataFrame({"perturbation": ["KO"], "max_abs_cosine_diff": [1e-7], "max_abs_embedding_rerun_diff": [1e-7]})
        self.assertTrue(judge_smoke(table, perturbations=["KO"])["pass"])

