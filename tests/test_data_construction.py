"""Small independent numerical oracles for data construction, including file I/O."""
import json
import sys
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import part5_pseudobulk as pb
import part5_source_blocks as sources
import part6_axes as axes
import part6_endpoints as endpoints
import part6_token_audit as tokens
import part6_eligibility as eligibility


class DataConstruction(unittest.TestCase):
    def fixture(self, sparse_counts=False):
        # Per donor: five cells; target totals 3, 0, 10 and detections 2, 0, 5.
        counts = np.array([[1, 0, 2], [2, 1, 2], [3, 0, 2], [4, 2, 2], [5, 0, 2],
                           [2, 0, 1]] + [[2, 0, 1]] * 4 + [[1, 2, 3]] * 5)
        obs = pd.DataFrame({"dataset": ["A"]*10 + ["B"]*5,
                            "donor_id": ["d1"]*5 + ["d2"]*5 + ["d1"]*5,
                            "subtype": ["S"]*15}, index=[f"c{i}" for i in range(15)])
        obj = ad.AnnData(X=np.zeros_like(counts), obs=obs,
                        var=pd.DataFrame(index=["G1", "TARGET_FEATURE", "target_feature"]))
        obj.layers["counts"] = sparse.csr_matrix(counts) if sparse_counts else counts
        return obj

    def build(self, obj, **overrides):
        args = dict(gene="TARGET_FEATURE", group_key="subtype", dataset_key="dataset",
                    donor_key="donor_id", unit_col="dataset_donor_id", min_cells=5,
                    include_group_in_unit=True)
        args.update(overrides)
        return pb.build_pseudobulk(obj, **args)

    def assert_gold(self, counts, meta, genes):
        self.assertEqual(genes.gene.tolist(), ["G1", "target_feature"])
        self.assertEqual(meta.unit_id.tolist(), ["A_d1|S", "A_d2|S", "B_d1|S"])
        np.testing.assert_array_equal(counts.toarray(), [[15, 10, 5], [10, 5, 15]])
        np.testing.assert_array_equal(meta.n_cells, [5, 5, 5])
        np.testing.assert_array_equal(meta.gene_umi, [3, 0, 10])
        np.testing.assert_array_equal(meta.n_detected, [2, 0, 5])
        np.testing.assert_allclose(meta.detection_fraction, [0.4, 0, 1])
        np.testing.assert_allclose(meta.jeffreys_per_10pct, [25/6, 5/6, 55/6])
        np.testing.assert_array_equal(meta.total_umi, [28, 15, 30])
        np.testing.assert_allclose(meta.mean_umi_per_cell, [5, 3, 4])
        np.testing.assert_allclose(meta.log2cpm, np.log2([1 + 3e6/28, 1, 1 + 10e6/30]))
        self.assertTrue(meta.eligible.all())

    def test_pseudobulk_hand_counts_dense_sparse_and_cell_permutation(self):
        for sparse_counts in (False, True):
            obj = self.fixture(sparse_counts)
            self.assert_gold(*self.build(obj)[:3])
            order = [4, 9, 14, 0, 5, 10, 1, 6, 11, 2, 7, 12, 3, 8, 13]
            self.assert_gold(*self.build(obj[order].copy())[:3])

    def test_pseudobulk_real_h5ad_to_matrix_market(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.fixture(True).write_h5ad(root / "input.h5ad")
            argv = [str(root / "input.h5ad"), "--gene", "TARGET_FEATURE", "--min-cells", "5",
                    "--out", str(root / "out")]
            self.assertEqual(pb.main(argv), 0)
            out = root / "out"
            self.assert_gold(mmread(out / "counts.mtx").tocsr(), pd.read_csv(out / "metadata.csv"),
                             pd.read_csv(out / "genes.csv"))
            audit = json.loads((out / "pseudobulk_audit.json").read_text())
            self.assertFalse(audit["target_in_outcome_matrix"])
            self.assertEqual(audit["outcome_gene_manifest"], "genes.csv")
            self.assertFalse((out / "target_excluded_genes.csv").exists())
            with self.assertRaises(SystemExit):
                pb.main(argv)

    def test_pseudobulk_rejects_ambiguous_identities_and_counts(self):
        for kind in ("barcode", "gene", "missing", "collision", "split", "fractional"):
            obj = self.fixture()
            if kind == "barcode":
                obj.obs_names = ["same"]*15
            elif kind == "gene":
                obj.var_names = ["G1", "TARGET_FEATURE", "TARGET_FEATURE"]
            elif kind == "missing":
                obj.obs.loc["c0", "donor_id"] = None
            elif kind == "collision":
                obj.obs["dataset_donor_id"] = "merged"
            elif kind == "split":
                obj.obs["dataset_donor_id"] = obj.obs_names
            else:
                obj.layers["counts"] = obj.layers["counts"].astype(float)
                obj.layers["counts"][0, 0] = 0.5
            with self.subTest(kind=kind), self.assertRaises(SystemExit):
                self.build(obj)

    def test_subtype_split_and_explicit_pooling(self):
        obj = self.fixture()
        obj.obs.loc["c0", "subtype"] = "S2"
        counts, meta, _, _ = self.build(obj)
        self.assertEqual(meta.n_cells.tolist(), [1, 4, 5, 5])
        np.testing.assert_array_equal(counts.toarray(), [[1, 14, 10, 5], [2, 8, 5, 15]])
        counts, meta, _, _ = self.build(obj, include_group_in_unit=False)
        self.assertEqual(meta.subtype.tolist(), ["POOLED", "S", "S"])
        np.testing.assert_array_equal(counts.toarray(), [[15, 10, 5], [10, 5, 15]])

    def test_source_variance_and_loo_identity(self):
        share = sources.between_within_fraction(np.array([0., 2., 4., 6.]), np.array(["A", "A", "B", "B"]))
        self.assertAlmostEqual(share["between_fraction"], .8)
        self.assertAlmostEqual(share["within_fraction"], .2)
        meta = pd.DataFrame({"dataset_donor_id": ["d1", "d1", "d2"],
                             "dataset": ["A", "A", "B"], "source_block": ["X", "X", "Y"]})
        loo = sources.relabel_loo(pd.DataFrame({"omitted_donor": ["d2", "d1"]}), meta, "dataset")
        self.assertEqual(loo.source_block.tolist(), ["Y", "X"])
        with self.assertRaises(SystemExit):
            sources.relabel_loo(pd.DataFrame({"omitted_donor": ["absent"]}), meta, "dataset")
        meta.loc[1, "source_block"] = "Y"
        with self.assertRaises(SystemExit):
            sources.relabel_loo(loo, meta, "dataset")

    def test_source_cli_frozen_map_and_boolean_eligibility(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            meta = pd.DataFrame({"dataset": ["A", "A", "B", "B", "B"],
                                 "donor_id": ["1", "2", "1", "2", "3"],
                                 "jeffreys_per_10pct": [0, 2, 4, 6, 100],
                                 "eligible": ["true", "1", "true", "1", "false"], "n_cells": [5]*5})
            meta.to_csv(root / "meta.csv", index=False)
            (root / "map.yaml").write_text("A: STUDY\nB: STUDY\n")
            argv = [str(root / "meta.csv"), "--map", str(root / "map.yaml"), "--out", str(root / "out")]
            self.assertEqual(sources.main(argv), 0)
            mapped = pd.read_csv(root / "out/metadata_with_source_block.csv")
            self.assertEqual(len(mapped), 4)
            self.assertEqual(mapped.source_block.unique().tolist(), ["STUDY"])
            ranges = pd.read_csv(root / "out/exposure_range_by_dataset.csv")
            np.testing.assert_allclose(ranges["mean"], [1, 5])
            np.testing.assert_allclose(ranges.sd, [np.sqrt(2)]*2)
            self.assertEqual(ranges.n_donors.tolist(), [2, 2])
            (root / "map.yaml").write_text("A: STUDY\n")
            with self.assertRaises(ValueError):
                sources.main(argv)

    def test_axis_geometry_coincident_direction_and_heldout_invariance(self):
        v = 1 / np.sqrt(2)
        donor_axes = {"A": np.array([v, v]), "B": np.array([1., 0]), "C": np.array([0., 1])}
        loo = axes.loo_axis(donor_axes, "A")
        np.testing.assert_allclose(loo, [v, v])
        axes.assert_no_self_leakage(donor_axes, "A", loo)
        donor_axes["A"] = np.array([-1., 0])
        np.testing.assert_allclose(axes.loo_axis(donor_axes, "A"), loo)
        np.testing.assert_allclose(axes.donor_axis(np.array([[2., 0], [0, 3.]]), np.array([0]), np.array([1])), [v, -v])
        np.testing.assert_allclose(axes.delta_axis(np.array([[0., 5.]]), np.array([[2., 0.]]), np.array([1., 0.])), [-1])
        with self.assertRaises(ValueError):
            axes.loo_axis({"A": np.array([0., 1.]), "B": np.array([1., 0.]), "C": np.array([-1., 0.])}, "A")

    def test_build_axes_uses_donor_equal_weights_and_unique_cells(self):
        cells = pd.DataFrame({"cell_id": [f"c{i}" for i in range(8)], "dataset_donor_id": ["A"]*2 + ["B"]*2 + ["C"]*4})
        cls = np.array([[-1, -1], [1, 1], [-1, 0], [1, 0], [0, -1], [0, -1], [0, 1], [0, 1]], float)
        scores = pd.Series([0, 1, 0, 1, 0, 0, 1, 1])
        result = axes.build_axes(cls, cells, scores, min_side=1, min_cells=2, min_training=2)
        np.testing.assert_allclose(result["A"], [1/np.sqrt(2)]*2)
        cells.loc[1, "cell_id"] = "c0"
        with self.assertRaises(SystemExit):
            axes.build_axes(cls, cells, scores, min_side=1, min_cells=2, min_training=2)

    def test_endpoint_coverage_target_exclusion_and_duplicate_rejection(self):
        args = dict(endpoint_id="E", role="primary", members=["TARGET_FEATURE", "A", "B", "C"],
                    model_visible={"TARGET_FEATURE", "A", "B"}, target="TARGET_FEATURE", min_frac=2/3,
                    min_genes=2, path="set.txt", sha256="test")
        row = endpoints.coverage_row(**args)
        self.assertEqual((row["n_frozen"], row["n_model_visible"], row["n_target_dropped"]), (3, 2, 1))
        self.assertAlmostEqual(row["coverage_frac"], 2/3)
        self.assertEqual(row["status"], "PASS")
        args["min_genes"] = 3
        self.assertEqual(endpoints.coverage_row(**args)["status"], "STOPPED")
        args["role"] = "support"
        self.assertEqual(endpoints.coverage_row(**args)["status"], "NOT_ESTIMABLE")
        args["members"] = ["A", "A", "B"]
        with self.assertRaises(SystemExit):
            endpoints.coverage_row(**args)

    def test_axes_cli_aligns_shuffled_score_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cells = pd.DataFrame({"cell_id": [f"c{i}" for i in range(110)],
                                  "dataset_donor_id": np.repeat([f"d{i}" for i in range(11)], 10)})
            cls = np.tile(np.array([[-1., 0.]]*5 + [[1., 0.]]*5), (11, 1))
            scores = pd.DataFrame({"cell_id": cells.cell_id, "E": np.tile(np.arange(10), 11)})
            np.save(root / "cls.npy", cls)
            cells.to_csv(root / "cells.csv", index=False)
            scores.iloc[::-1].to_csv(root / "scores.csv", index=False)
            argv = ["--cls", str(root / "cls.npy"), "--cells", str(root / "cells.csv"),
                    "--scores", str(root / "scores.csv"), "--endpoint", "E", "--out", str(root / "axes.npz")]
            self.assertEqual(axes.main(argv), 0)
            with np.load(root / "axes.npz") as result:
                self.assertEqual(len(result.files), 11)
                for key in result.files:
                    np.testing.assert_allclose(result[key], [1., 0.])
            scores.drop(columns="cell_id").to_csv(root / "scores.csv", index=False)
            with self.assertRaises(SystemExit):
                axes.main(argv)

    def test_token_classes_and_cap_boundary_with_string_booleans(self):
        frame = pd.DataFrame({"cell_id": list("abcdef"), "raw_count": [1, 0, 1, 1, 0, 1],
                              "final_token_present": ["true", "false", "false", "false", "true", "false"],
                              "sequence_length": [10, 10, 4096, 4096, 10, 4095],
                              "pretruncation_rank": [1, np.nan, 4095, 4094, np.nan, 4095]})
        result = tokens.audit_ledger(frame)
        self.assertEqual(result.token_class.tolist(), [tokens.KO_OK, tokens.OE_ONLY, tokens.TRUNCATION,
                                                       tokens.ILLEGAL, tokens.ILLEGAL, tokens.ILLEGAL])
        self.assertEqual(tokens.summary(result), dict(n_cells=6, n_ko_eligible=1, n_oe_primary_eligible=3, n_truncation=1, n_illegal=3))
        for field, value in [("raw_count", -1), ("raw_count", np.nan), ("pretruncation_rank", np.inf),
                             ("final_token_present", "maybe"), ("cell_id", "b")]:
            bad = frame.copy()
            bad.loc[0, field] = value
            with self.assertRaises(SystemExit):
                tokens.audit_ledger(bad)

    def test_eligibility_counts_distinct_successful_cells(self):
        rows = []
        for donor, pert, pop, count in [("A", "KO", "PRIMARY", 5), ("B", "KO", "PRIMARY", 4),
                                        ("C", "OE", "PRIMARY", 9), ("D", "OE", "SYMMETRY", 5)]:
            for i in range(count):
                rows.append(dict(cell_id=f"{donor}{i}", dataset_donor_id=donor, target_symbol="TARGET_FEATURE",
                                 perturbation=pert, analysis_population=pop, endpoint="E", delta_axis=i, run_status="RUN"))
        cells = pd.DataFrame(rows)
        ledger, valid = eligibility.donor_eligibility(cells)
        self.assertEqual(ledger.n_success.tolist(), [5, 4, 9, 5])
        self.assertEqual(valid.dataset_donor_id.tolist(), ["A", "D"])
        self.assertEqual(valid.median_delta_axis.tolist(), [2, 2])
        self.assertEqual(valid.q1_delta_axis.tolist(), [1, 1])
        self.assertEqual(valid.q3_delta_axis.tolist(), [3, 3])
        with self.assertRaises(SystemExit):
            eligibility.donor_eligibility(pd.concat([cells, cells.iloc[:1]]))
        cells["run_status"] = "FAILED"
        ledger, valid = eligibility.donor_eligibility(cells)
        self.assertTrue(valid.empty)
        self.assertEqual(valid.columns.tolist(), ledger.columns.tolist())


if __name__ == "__main__":
    unittest.main()
