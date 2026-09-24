import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from part5_eligibility import design_matrix, evaluate, flag_arm, main as eligibility_main
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

    def _eligibility_frame(self, n_donors=12, subtypes=1, *, eligible=True, source=True):
        rows = []
        for i in range(n_donors):
            for s in range(subtypes):
                row = {
                    "subtype": f"S{s}",
                    "dataset": f"D{i % 3}",
                    "dataset_donor_id": f"d{i:03d}",
                    "unit_id": f"d{i:03d}_S{s}",
                    "n_cells": 40,
                    "z_log1p_n_cells": float(i % 4) + 0.1,
                    "z_log1p_mean_umi": float((i * 2) % 5) + 0.2 * s,
                    "jeffreys_per_10pct": float(i) * 0.5 + 0.4 * s,
                }
                if eligible is not None:
                    row["eligible"] = "true"
                if source:
                    row["source_block"] = f"B{i % 3}"
                rows.append(row)
        return pd.DataFrame(rows)

    def _evaluate(self, frame):
        return evaluate(
            frame,
            arm_key="subtype",
            categorical=["subtype", "dataset"],
            numeric=["z_log1p_n_cells", "z_log1p_mean_umi", "jeffreys_per_10pct"],
            exposure="jeffreys_per_10pct",
            dataset_key="dataset",
            source_key="source_block",
            n_formal=12,
            n_exploratory=8,
            min_datasets_formal=3,
            min_rdf_formal=6,
            min_rdf_fit=4,
            sd_floor=0.5,
            drop_single=True,
            part4=None,
        )

    def test_missing_eligible_column_is_refused(self):
        frame = self._eligibility_frame()
        frame = frame.drop(columns=["eligible"])
        with self.assertRaises(SystemExit):
            self._evaluate(frame)

    def test_missing_source_block_is_refused_unless_declared(self):
        frame = self._eligibility_frame(source=False)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            meta = root / "meta.csv"
            frame.to_csv(meta, index=False)
            out = root / "eligibility.csv"
            with self.assertRaises(SystemExit):
                eligibility_main([str(meta), "--out", str(out)])
            self.assertFalse(out.exists())
            self.assertEqual(
                eligibility_main([str(meta), "--source-key", "dataset", "--out", str(out)]),
                0,
            )
            audit = json.loads(out.with_suffix(".audit.json").read_text(encoding="utf-8"))
            self.assertEqual(audit["source_key"], "dataset")

    def test_all_arm_counts_units_and_gates_on_donors(self):
        table = self._evaluate(self._eligibility_frame(n_donors=8, subtypes=2))
        pooled = table.loc[table["arm"].eq("ALL")].iloc[0]
        self.assertEqual(int(pooled.n_units), 16)
        self.assertEqual(int(pooled.n_donors), 8)
        self.assertEqual(pooled.status, "exploratory")

    def test_duplicate_donor_is_not_a_second_sign(self):
        frame = pd.DataFrame({"target_symbol": ["G"]*5, "endpoint": ["E"]*5,
            "perturbation": ["KO"]*5, "dataset_donor_id": ["d"]*5,
            "median_delta_axis": [.1]*5})
        with self.assertRaises(SystemExit):
            sign_tests(frame, target="G", family_size=1)


if __name__ == "__main__":
    unittest.main()

