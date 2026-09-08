import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from part5_eligibility import design_matrix, flag_arm
from part6_sign_tests import sign_tests


class ScientificContractTests(unittest.TestCase):
    def test_missing_covariate_is_not_imputed(self):
        with self.assertRaises(SystemExit):
            design_matrix(pd.DataFrame({"x": [1, np.nan, 3]}), [], ["x"])

    def test_accessions_do_not_replace_sources(self):
        status, _ = flag_arm(n_units=30, n_datasets=5, n_blocks=1, rdf=20,
            rank=3, n_cols=3, exposure_sd=1, n_formal=12, n_exploratory=8,
            min_datasets_formal=3, min_rdf_formal=6, min_rdf_fit=4,
            sd_floor=.5, part4_flag="")
        self.assertEqual(status, "exploratory")

    def test_duplicate_donor_is_not_a_second_sign(self):
        frame = pd.DataFrame({"target_symbol": ["G"]*5, "endpoint": ["E"]*5,
            "perturbation": ["KO"]*5, "dataset_donor_id": ["d"]*5,
            "median_delta_axis": [.1]*5})
        with self.assertRaises(SystemExit):
            sign_tests(frame, target="G", family_size=1)


if __name__ == "__main__":
    unittest.main()

