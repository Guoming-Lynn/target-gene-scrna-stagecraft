import ast
import json
import tempfile
import unittest
import unittest.mock
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

import stagecraft
import check_stage_layout
import generate_toy_data
import part1_init
import simulation_contract
from demo_toy_run import main as demo_main


class OpenSourceContract(unittest.TestCase):
    def test_versions_align(self):
        skill = (ROOT / "SKILL.md").read_text(encoding="utf-8")
        self.assertIn(f'version: "{stagecraft.__version__}"', skill)
        pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
        self.assertIn(f'version = "{stagecraft.__version__}"', pyproject)
        citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
        self.assertIn(f"version: {stagecraft.__version__}", citation)
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        self.assertIn(f"v{stagecraft.__version__}", readme)
        readme_cn = (ROOT / "README.zh-CN.md").read_text(encoding="utf-8")
        self.assertIn(f"v{stagecraft.__version__}", readme_cn)

    def test_dev_pins_match_and_geneformer_manifest_is_not_updated(self):
        dev = (ROOT / "requirements-dev.txt").read_text(encoding="utf-8")
        pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
        pin = "ruff==0.16.8"
        self.assertIn(pin, dev)
        self.assertIn(f'"{pin}"', pyproject)
        dependabot = (ROOT / ".github/dependabot.yml").read_text(encoding="utf-8")
        self.assertIn("references/part6-requirements.txt", dependabot)
        pins = (ROOT / "references/part6-requirements.txt").read_text(encoding="utf-8")
        self.assertIn("anndata==0.9.2", pins)

    def test_argparse_module_docs_encode_on_windows_console(self):
        paths = sorted((ROOT / "scripts").glob("*.py")) + [ROOT / "quickstart.py"]
        for path in paths:
            text = path.read_text(encoding="utf-8")
            if "description=__doc__" not in text:
                continue
            doc = ast.get_docstring(ast.parse(text)) or ""
            try:
                doc.encode("cp1252")
            except UnicodeEncodeError as exc:
                self.fail(f"{path.name} argparse description is not cp1252-safe: {exc}")

    def test_simulation_contract_is_manifest_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "manifest.json"
            self.assertEqual(
                simulation_contract.main(["--seed", "7", "--declared-replicates", "100", "--out", str(out)]),
                0,
            )
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["status"], "MANIFEST_ONLY")
            self.assertEqual(payload["replicates"], 100)
            self.assertIn("does not sample data", payload["note"])

    def test_check_stage_layout_has_help_and_rejects_missing_dir(self):
        with tempfile.TemporaryDirectory() as tmp:
            missing = Path(tmp) / "absent"
            self.assertEqual(check_stage_layout.main([str(missing)]), 1)

    def test_part1_init_writes_protocol_scaffold(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "stage"
            self.assertEqual(part1_init.main(["--out", str(out), "--gene", "DEMO"]), 0)
            protocol = (out / "00_protocol_manifest" / "PROTOCOL.md").read_text(encoding="utf-8")
            self.assertIn("DEMO", protocol)
            self.assertIn("NOT_ESTIMABLE", protocol)
            self.assertTrue((out / "02_tables").is_dir())
            with self.assertRaises(SystemExit):
                part1_init.main(["--out", str(out), "--gene", "DEMO"])

    def test_part1_init_rejects_blank_gene_before_writing(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "stage"
            with self.assertRaises(SystemExit):
                part1_init.main(["--out", str(out), "--gene", "  "])
            self.assertFalse(out.exists())

    def test_identity_csv_keeps_leading_zeros_and_bools_are_strict(self):
        from stagecraft.io import parse_bool_column, read_identity_csv
        import pandas as pd

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "meta.csv"
            pd.DataFrame({"donor_id": ["001", "002"], "eligible": ["true", "false"]}).to_csv(path, index=False)
            frame = read_identity_csv(path)
            self.assertEqual(frame.donor_id.tolist(), ["001", "002"])
            path.write_text("sample_id\n007\n008\n", encoding="utf-8")
            self.assertEqual(read_identity_csv(path).sample_id.tolist(), ["007", "008"])
            self.assertEqual(parse_bool_column(frame.eligible, "eligible").tolist(), [True, False])
            with self.assertRaises(SystemExit):
                parse_bool_column(pd.Series(["Yes", "No"]), "eligible")

    def test_generate_toy_data_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "toy.h5ad"
            self.assertEqual(generate_toy_data.main(["--out", str(path), "--n-donors", "2", "--cells-per-unit", "2"]), 0)
            with self.assertRaises(SystemExit):
                generate_toy_data.main(["--out", str(path)])

    def test_demo_toy_run_consumes_toy_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "demo"
            self.assertEqual(demo_main(["--out", str(out), "--seed", "7"]), 0)
            report = json.loads((out / "demo_report.json").read_text(encoding="utf-8"))
            self.assertEqual(report["status"], "PASS")
            self.assertFalse(report["formal_analysis"])
            self.assertTrue((out / "toy.h5ad").is_file())
            self.assertTrue((out / "02_tables" / "eligibility.csv").is_file())
            manifest = json.loads((out / "simulation_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["status"], "MANIFEST_ONLY")

    def test_multi_file_publish_removes_partials_after_a_writer_failure(self):
        from stagecraft.io import publish_new_files, require_new

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            dest = root / "a.txt"

            def write(partials):
                partials[0].write_text("ok\n", encoding="utf-8")
                raise RuntimeError("boom")

            with self.assertRaises(RuntimeError):
                publish_new_files([dest], write)
            self.assertFalse(dest.exists())
            self.assertEqual(list(root.iterdir()), [])
            reserved = require_new(root / "b.txt")
            self.assertEqual(reserved.read_bytes(), b"")
            with self.assertRaises(SystemExit):
                require_new(reserved)

    def test_publish_rolls_back_a_renamed_final(self):
        import os
        from stagecraft.io import publish_new_files

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = root / "a.txt"
            second = root / "b.txt"
            real = os.rename
            calls = {"n": 0}

            def rename(src, dst):
                calls["n"] += 1
                if calls["n"] > 1:
                    raise OSError("boom")
                return real(src, dst)

            def write(partials):
                partials[0].write_text("a\n", encoding="utf-8")
                partials[1].write_text("b\n", encoding="utf-8")

            with self.assertRaises(OSError):
                with unittest.mock.patch("stagecraft.io.os.rename", rename):
                    publish_new_files([first, second], write)
            self.assertFalse(first.exists())
            self.assertFalse(second.exists())
            self.assertEqual(list(root.iterdir()), [])


if __name__ == "__main__":
    unittest.main()
