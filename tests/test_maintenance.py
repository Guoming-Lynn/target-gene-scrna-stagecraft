import hashlib
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import hash_inputs as hashes
import check_protocol as protocol
import cluster_review_tables as review
import part2_figures as part2
import part3_figures as part3
from plotting_style import load_plotting_config, save_figure
import validate_figure_manifest as figure_manifest


class Maintenance(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def fixture(self):
        obs = pd.DataFrame({"dataset": ["A", "B", "A", "B"], "donor_id": ["1"]*4,
                            "cluster": ["0", "0", "1", "1"]}, index=list("abcd"))
        obj = ad.AnnData(np.zeros((4, 1)), obs=obs, var=pd.DataFrame(index=["FEATURE"]))
        obj.layers["counts"] = np.array([[0], [2], [1], [0]])
        obj.layers["normalized"] = np.array([[0.], [4.], [3.], [0.]])
        obj.obsm["X_umap"] = np.array([[0., 0.], [1., 2.], [3., 1.], [4., 4.]])
        return obj

    def test_hash_verification_mutation_missing_and_invalid_digest(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "data.bin"
            path.write_bytes(b"abc")
            expected = hashlib.sha256(b"abc").hexdigest()
            self.assertEqual(hashes.sha256(path), expected)
            for obj in ({"path": "data.bin", "sha256": expected.upper()},
                        {"input": "data.bin", "input_sha256": expected}):
                self.assertEqual(hashes.verify_entries(obj, root)[0]["status"], "verified")
            for obj in ({"path": "data.bin", "sha256": "replace"},
                        {"path": "absent", "sha256": expected}, {}):
                with self.assertRaises(ValueError):
                    hashes.verify_entries(obj, root)
            path.write_bytes(b"changed")
            with self.assertRaises(ValueError):
                hashes.verify_entries({"path": "data.bin", "sha256": expected}, root)

    def test_hash_scan_bare_filenames_and_empty_manifest_fail(self):
        found = []
        hashes.collect({"inputs": ["a.csv", "b.bin"], "input": "c.h5ad"}, found)
        self.assertEqual(found, ["a.csv", "b.bin", "c.h5ad"])
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            folder = root / "00_protocol_manifest"
            folder.mkdir()
            manifest = folder / "inputs.json"
            manifest.write_text('{}')
            self.assertEqual(hashes.manifest_root(manifest), root.resolve())
            self.assertEqual(hashes.main(["hash_inputs", str(manifest)]), 1)

    def test_protocol_requires_headings_and_bodies(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "PROTOCOL.md"
            sections = ["Question", "Hard boundaries", "Inputs", "Statistical unit", "Model", "Verdict", "Not done"]
            body = "Frozen design details and explicit assumptions for this stage. " * 5
            valid = "\n".join(f"## {heading}\n{body}" for heading in sections)
            valid += "\nEvidence ceiling exploratory. Forbidden sentences. Upstream is read-only. NOT_ESTIMABLE."
            path.write_text(valid)
            self.assertEqual(protocol.main(["check_protocol", str(path)]), 0)
            for invalid in (valid.replace("## ", ""), valid.replace("## Model\n" + body, "## Model\n"), "short"):
                path.write_text(invalid)
                self.assertEqual(protocol.main(["check_protocol", str(path)]), 1)

    def test_cluster_qc_donors_and_strict_marker_filter(self):
        obj = self.fixture()
        markers = pd.DataFrame({"cluster": ["0", "0", "1"], "gene": ["FEATURE", "OTHER", "THIRD"],
                                "rank": [1, 2, 1], "score": [1., -1., 1.], "logfoldchange": [1.]*3,
                                "pval_adj": [.01, .01, .2], "pct_nz_group": [.5]*3, "pct_nz_reference": [.1]*3})
        filtered = review.strict_positive(markers)
        self.assertEqual(filtered.gene.tolist(), ["FEATURE"])
        qc = review.cluster_qc(obj, "cluster", filtered, {"FEATURE"}, "FEATURE")
        self.assertEqual(qc.n_donors.tolist(), [2, 2])
        self.assertEqual(qc.n_cells.tolist(), [2, 2])
        self.assertEqual(qc.proposed_label.tolist(), ["", ""])
        np.testing.assert_allclose(qc.target_gene_detection_fraction, [.5, .5])

    def test_part2_summary_and_feature_artist_alignment(self):
        obj = self.fixture()
        summary = part2.summarize_groups(obj, "FEATURE", "cluster")
        self.assertEqual(summary.n_donors.tolist(), [2, 2])
        np.testing.assert_allclose(summary.positive_fraction, [.5, .5])
        np.testing.assert_allclose(summary.mean_normalized_all, [2, 1.5])
        before = obj.layers["normalized"].copy()
        part2.plot_umap_feature(obj, "FEATURE", load_plotting_config(), clip_percentile=None)
        fig = plt.gcf()
        artist = fig.axes[0].collections[0]
        np.testing.assert_allclose(artist.get_array(), [0, 0, 3, 4])
        np.testing.assert_allclose(artist.get_offsets(), obj.obsm["X_umap"][[0, 3, 2, 1]])
        fig.canvas.draw()
        self.assertGreater(np.asarray(fig.canvas.buffer_rgba()).std(), 0)
        np.testing.assert_array_equal(obj.layers["normalized"], before)

    def test_review_export_retains_cluster_without_positive_markers(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.fixture().write_h5ad(root / "input.h5ad")
            ranked = pd.DataFrame({"cluster": ["0", "1"], "gene": ["FEATURE", "FEATURE"],
                                   "rank": [1, 1], "score": [2., -2.], "logfoldchange": [1., -1.],
                                   "pval_adj": [.01, .9]})
            argv = ["cluster_review_tables", str(root / "input.h5ad"), "--leiden-key", "cluster",
                    "--strict-positive", "--out", str(root / "out")]
            with patch.object(sys, "argv", argv), patch.object(review, "markers", return_value=ranked):
                self.assertEqual(review.main(), 0)
            table = pd.read_csv(root / "out/cluster_top20_by_cluster.csv", keep_default_na=False,
                                dtype={"cluster": str})
            self.assertEqual(table.cluster.tolist(), ["0", "1"])
            self.assertEqual(table.n_cells.tolist(), [2, 2])
            self.assertEqual(table.top20_markers.tolist(), ["FEATURE", ""])
            self.assertEqual(table.proposed_label.tolist(), ["", ""])

    def test_part3_marker_values_missing_layer_and_numeric_membership(self):
        obj = self.fixture()
        evidence, missing = part3.marker_matrix(obj, "cluster", ["FEATURE", "ABSENT"])
        self.assertEqual(missing, ["ABSENT"])
        np.testing.assert_allclose(evidence.detection, [.5, .5])
        np.testing.assert_allclose(evidence.mean_normalized, [2, 1.5])
        del obj.layers["normalized"]
        with self.assertRaises(ValueError):
            part3.marker_matrix(obj, "cluster", ["FEATURE"])
        membership = pd.DataFrame({"parent_cluster": [0, 0, 1], "child_cluster": [0, 1, 1], "n_cells": [2, 3, 4]})
        with patch.object(part3, "save_figure") as save:
            part3.plot_membership_heatmap(membership, load_plotting_config(), Path("unused"))
            fig = save.call_args.args[0]
            np.testing.assert_array_equal(fig.axes[0].images[0].get_array(), [[2, 3], [0, 4]])
            fig.canvas.draw()
            self.assertGreater(np.asarray(fig.canvas.buffer_rgba()).std(), 0)

    def test_figure_bundle_is_immutable_and_statistics_manifest_is_valid(self):
        statistics = {
            "estimand": "donor-unit association",
            "biological_unit": "dataset × donor_id",
            "n_definition": "eligible donor-units",
            "n": 12,
            "model_or_test": "limma-voom QW",
            "effect_scale": "log2 fold change",
            "error_bar": "95% CI",
            "multiple_testing_family": "all target-excluded genes",
            "adjustment": "BH",
            "claim_ceiling": "association only",
        }
        with tempfile.TemporaryDirectory() as tmp:
            stem = Path(tmp) / "F05_08"
            fig, ax = plt.subplots()
            ax.plot([0, 1], [0, 1])
            # Export drivers vary across local desktop installations; this
            # contract test isolates the writer's paths and immutability logic.
            with patch.object(fig, "savefig", side_effect=lambda path, **_kwargs: Path(path).touch()):
                written = save_figure(
                    fig, stem, load_plotting_config(), parameters={"plot": "test"},
                    statistics=statistics, formats=["png"]
                )
            self.assertEqual(figure_manifest.main([str(written["statistics"])]), 0)
            self.assertTrue(written["parameters"].is_file())
            self.assertTrue(written["statistics"].is_file())
            with self.assertRaises(FileExistsError):
                save_figure(plt.figure(), stem, load_plotting_config())

    def test_figure_rejects_incomplete_statistics_before_writing(self):
        with tempfile.TemporaryDirectory() as tmp:
            stem = Path(tmp) / "F05_09"
            with self.assertRaises(ValueError):
                save_figure(plt.figure(), stem, load_plotting_config(), statistics={"estimand": "x"})
            self.assertEqual(list(Path(tmp).iterdir()), [])


if __name__ == "__main__":
    unittest.main()
