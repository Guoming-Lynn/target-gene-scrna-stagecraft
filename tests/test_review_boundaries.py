import itertools
import math
import os
import sys
import tempfile
import subprocess
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from part5_verdict import guarded_verdict
from part6_sign_tests import median_order_statistic_interval, sign_tests
from protocol_chronology import check_chronology
from part6_controls import select_controls


class ReviewBoundaries(unittest.TestCase):
    def test_frozen_project_control_exclusion(self):
        frame = pd.DataFrame({"gene_symbol": ["TARGET_GENE"] + [f"FEATURE{i}" for i in range(20)],
                              "detection_fraction": [.5] * 21,
                              "mean_raw_counts_per_cell": [10.] + [9. + i * .1 for i in range(20)],
                              "in_model_vocabulary": [True] * 21})
        selected, candidates, _ = select_controls(frame, target="TARGET_GENE", endpoint_union=set(),
                                                   excluded_genes={"FEATURE10"})
        self.assertEqual(len(selected), 10)
        self.assertNotIn("FEATURE10", selected.gene_symbol.tolist())
        self.assertEqual(candidates.set_index("gene_symbol").loc["FEATURE10", "exclusion_reason"],
                         "frozen_project_exclusion")

    def test_figure_cli_parses_arguments(self):
        script = Path(__file__).resolve().parents[1] / "scripts" / "part6_figures.py"
        result = subprocess.run([sys.executable, str(script), "--help"],
                                capture_output=True, text=True, timeout=60, check=False)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("forest", result.stdout)

    def test_freeze_rejects_existing_output_and_detects_changes(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            protocol = root / "PROTOCOL.md"
            protocol.write_text("frozen analysis")
            check_chronology(protocol, root, freeze=True)
            self.assertIn("LOCAL_SEQUENCE_CONSISTENT", check_chronology(protocol, root))
            protocol.write_text("changed after looking")
            with self.assertRaises(ValueError):
                check_chronology(protocol, root)

    def test_preexisting_artifact_cannot_be_backdated_by_freeze(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            protocol = root / "PROTOCOL.md"
            protocol.write_text("protocol")
            (root / "02_tables").mkdir()
            output = root / "02_tables/gene.csv"
            output.write_text("result")
            with self.assertRaises(ValueError):
                check_chronology(protocol, root, freeze=True)
            output.unlink()
            check_chronology(protocol, root, freeze=True)
            output.write_text("result")
            os.utime(output, (1, 1))
            with self.assertRaises(ValueError):
                check_chronology(protocol, root)

    def test_interval_small_n_must_be_unbounded(self):
        low, high, coverage, status = median_order_statistic_interval([1, 2, 3, 4, 5])
        self.assertEqual((low, high), (-math.inf, math.inf))
        self.assertEqual((coverage, status), (1, "UNBOUNDED_SMALL_N"))

    def test_interval_exact_coverage_by_enumeration(self):
        # Enumerate all equally probable sign patterns under a continuous null.
        for n in (6, 8, 10, 12):
            covered = 0
            for signs in itertools.product((-1, 1), repeat=n):
                values = [sign * (i + 1) for i, sign in enumerate(signs)]
                low, high, coverage, _ = median_order_statistic_interval(values)
                covered += low <= 0 <= high
            self.assertAlmostEqual(covered / 2**n, coverage)
            self.assertGreaterEqual(coverage, .95)

    def test_zero_donors_retained_for_median_interval(self):
        df = pd.DataFrame({"target_symbol": ["G"] * 6, "endpoint": ["E"] * 6,
                           "perturbation": ["KO"] * 6, "dataset_donor_id": list("abcdef"),
                           "median_delta_axis": [0, 0, 0, 0, 0, 1]})
        result = sign_tests(df, target="G", family_size=1).iloc[0]
        self.assertEqual(result["median_ci_status"], "FINITE")
        self.assertEqual(result["median_ci_lower"], 0)
        self.assertTrue(math.isnan(result["raw_p_value"]))

    def test_external_tokens_rejected_even_in_unreachable_row(self):
        for token in ("EXTERNALLY_REPLICATED", "EXTERNALLY_CONSISTENT_UNDERPOWERED"):
            with self.assertRaises(SystemExit):
                guarded_verdict([{"token": "INCONCLUSIVE"}, {"token": token}], {})

    def test_formal_is_not_calibration(self):
        verdict = guarded_verdict([{"token": "FROZEN_PASS"}], {
            "source_dependent": False, "formal_gates_pass": True, "lodo_complete": True})
        self.assertEqual(verdict["verdict"], "FROZEN_PASS")
        self.assertFalse(verdict["scientifically_calibrated"])
        self.assertFalse(verdict["whole_pipeline_calibrated"])
        self.assertEqual(verdict["inference_scope"], "WITHIN_ARM_ONLY")
        self.assertEqual(verdict["project_multiplicity_control"], "NOT_ESTABLISHED")

    def test_joint_mode_cannot_get_custom_formal_pass(self):
        verdict = guarded_verdict([{"token": "FROZEN_PASS"}], {
            "estimand_mode": "joint_common_slope", "source_dependent": False,
            "formal_gates_pass": True, "lodo_complete": True, "evidence_ceiling": "formal"})
        self.assertEqual(verdict["verdict"], "INCONCLUSIVE")
        self.assertEqual(verdict["evidence_ceiling"], "exploratory")


if __name__ == "__main__":
    unittest.main()
