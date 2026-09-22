import ast
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

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


if __name__ == "__main__":
    unittest.main()
