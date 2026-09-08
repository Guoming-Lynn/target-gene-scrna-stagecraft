import json
import sys
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from part3_prepare_removal import clean_analysis_state, load_decisions, partition, main as removal_main
from part4_identifiability import identifiability_table
from project_arm_inventory import inventory


class StageHandoffs(unittest.TestCase):
    def units(self):
        return pd.DataFrame({"subtype": ["S"]*12, "dataset": ["A", "B"]*6,
                             "dataset_donor_id": [f"d{i}" for i in range(12)], "n_cells": [30]*12,
                             "positive_fraction": np.linspace(.1,.8,12),
                             "log2_cpm_plus_1": np.linspace(1,4,12), "eligible": [True]*12})

    def test_duplicate_donor_does_not_inflate_forecast(self):
        units = self.units()
        with self.assertRaises(SystemExit):
            identifiability_table(pd.concat([units, units.iloc[:1]]), "subtype")

    def test_false_strings_and_cell_floor_retained_as_low_n(self):
        units = self.units()
        units["eligible"] = "False"
        result = identifiability_table(units, "subtype").iloc[0]
        self.assertEqual(result.flag, "LOW_N")
        self.assertEqual(result.n_units_eligible, 0)
        units["eligible"] = "True"
        units["n_cells"] = 19
        self.assertEqual(identifiability_table(units, "subtype").iloc[0].n_units_eligible, 0)

    def test_bad_exposure_and_missing_identity_rejected(self):
        for field, value in [("positive_fraction", np.nan), ("positive_fraction", 2),
                             ("dataset_donor_id", None), ("eligible", "maybe")]:
            units = self.units()
            units[field] = units[field].astype(object)
            units.loc[0,field] = value
            with self.assertRaises(SystemExit):
                identifiability_table(units, "subtype")

    def test_source_map_changes_scope_not_donor_count(self):
        units = self.units()
        units["source_block"] = "ONE"
        pooled = identifiability_table(units, "subtype", source_key="source_block")
        row = pooled.loc[pooled.source.eq("ALL")].iloc[0]
        self.assertEqual((row.flag,row.n_units_eligible), ("WITHIN_SOURCE_RANGE",12))

    def test_partition_aligns_barcodes_and_clears_old_analysis(self):
        raw = ad.AnnData(np.arange(12).reshape(4,3).astype(float))
        raw.obs_names = ["a","b","c","d"]
        raw.layers["counts"] = raw.X.copy()
        raw.var["highly_variable"] = True
        raw.varm["PCs"] = np.ones((3,2))
        raw.obsm["X_umap"] = np.ones((4,2))
        raw.raw = raw.copy()
        clustered = raw[[3,1,0,2]].copy()
        clustered.obs["leiden"] = ["1","0","0","1"]
        removed, retained, _, _ = partition(raw, clustered, "leiden", {"1"})
        self.assertEqual(removed.obs_names.tolist(), ["c","d"])
        clean_analysis_state(retained)
        np.testing.assert_array_equal(retained.X, raw.layers["counts"][:2])
        self.assertFalse(retained.varm.keys())
        self.assertFalse(retained.obsm.keys())
        self.assertNotIn("highly_variable",retained.var)
        self.assertIsNone(retained.raw)
        with self.assertRaises(SystemExit):
            partition(raw, clustered, "leiden", {"0","1"})
        with self.assertRaises(SystemExit):
            partition(raw[:3], clustered, "leiden", {"1"})

    def test_delete_needs_review_evidence(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp)/"decisions.csv"
            frame = pd.DataFrame({"selected_leiden_column":["leiden"]*2,"cluster_id":["0","1"],
                                  "decision":["KEEP","DELETE"],"reason":["retained","contamination"],
                                  "reviewer":["reviewer",""],"review_date":["2026-01-01"]*2})
            frame.to_csv(p,index=False)
            with self.assertRaises(ValueError):
                load_decisions(p,"leiden",{"0","1"})

    def test_removal_writes_counts_and_complete_partition(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            raw = ad.AnnData(np.arange(12).reshape(4,3).astype(float))
            raw.obs_names = ["a","b","c","d"]
            raw.layers["counts"] = raw.X.copy()
            raw.write_h5ad(root/"raw.h5ad")
            clustered = raw[[3,1,0,2]].copy()
            clustered.obs["leiden"] = ["1","0","0","1"]
            clustered.write_h5ad(root/"clustered.h5ad")
            pd.DataFrame({"selected_leiden_column":["leiden"]*2,"cluster_id":["0","1"],
                          "decision":["KEEP","DELETE"],"reason":["retained","contamination"],
                          "reviewer":["reviewer"]*2,"review_date":["2026-01-01"]*2}).to_csv(root/"decision.csv",index=False)
            args = ["--parent-raw",str(root/"raw.h5ad"),"--parent-clustered",str(root/"clustered.h5ad"),
                    "--leiden-key","leiden","--decision",str(root/"decision.csv"),
                    "--removed-out",str(root/"removed.h5ad"),"--child-raw",str(root/"child.h5ad"),
                    "--tables-out",str(root/"tables")]
            self.assertEqual(removal_main(args),0)
            child = ad.read_h5ad(root/"child.h5ad")
            np.testing.assert_array_equal(child.layers["counts"], raw.X[:2])
            ledger = pd.read_csv(root/"tables/parent_round_partition_manifest.csv")
            self.assertEqual(len(ledger),4)
            self.assertEqual(ledger.status.value_counts().to_dict(),{"retained":2,"removed":2})
            with self.assertRaises(FileExistsError):
                removal_main(args)

    def test_inventory_exposes_missing_and_failed_arms(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root/"pass.json").write_text(json.dumps({"status":"SUCCESS","verdict":"FROZEN_PASS"}))
            (root/"failed.json").write_text(json.dumps({"status":"FAILED"}))
            manifest = root/"arms.yaml"
            manifest.write_text('arms:\n- {id: A, verdict: pass.json}\n- {id: B, verdict: failed.json}\n- {id: C, verdict: absent.json}\n')
            result = inventory(manifest)
            self.assertEqual(result["status"], "INCOMPLETE_DECLARED_INVENTORY")
            self.assertEqual((result["n_declared_arms"],result["n_results_present"],result["n_frozen_pass"]),(3,2,1))
            self.assertEqual(result["project_multiplicity_control"], "NOT_ESTABLISHED")
            manifest.write_text('arms:\n- {id: A, verdict: pass.json}\n- {id: B, verdict: pass.json}\n')
            with self.assertRaises(ValueError):
                inventory(manifest)


if __name__ == "__main__":
    unittest.main()
