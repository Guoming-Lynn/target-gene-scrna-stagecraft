"""Claim lint, stage report, and stage status helpers."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from protocol_chronology import check_chronology
from stagecraft import EXIT_GATE, EXIT_USAGE
from stagecraft.claims import lint_text, load_rules
import claim_lint
import stage_report
import stage_status


class ClaimLintTests(unittest.TestCase):
    def test_every_packaged_example_hits_its_rule(self):
        rules, settings = load_rules()
        for rule in rules:
            with self.subTest(rule=rule.id):
                hits = {item.rule_id for item in lint_text(rule.example, rules, settings)}
                self.assertIn(rule.id, hits)

    def test_combined_sentence_reports_line_and_column(self):
        text = "The association was independently replicated in 12/12 donors."
        findings = {item.rule_id: item for item in lint_text(text)}
        self.assertEqual(
            set(findings),
            {"independent_replication", "donor_fraction_as_replication"},
        )
        self.assertEqual(findings["independent_replication"].line, 1)
        self.assertEqual(
            findings["independent_replication"].column,
            text.find("independently") + 1,
        )
        self.assertEqual(findings["donor_fraction_as_replication"].column, text.find("12/12") + 1)

    def test_negation_suppresses_english_and_chinese(self):
        self.assertEqual(lint_text("This is not independently replicated."), [])
        self.assertEqual(lint_text("不能写成独立重复"), [])

    def test_code_and_allow_lines_are_ignored(self):
        text = "\n".join(
            [
                "See `independently replicated` in the spec.",
                "```",
                "independently replicated",
                "```",
                "independently replicated claim-lint: allow",
            ]
        )
        self.assertEqual(lint_text(text), [])

    def test_chinese_causal_sentence_is_flagged(self):
        hits = {item.rule_id for item in lint_text("该基因驱动了纤维化")}
        self.assertIn("causal_language", hits)

    def test_harmless_phrases_are_not_flagged(self):
        self.assertEqual(lint_text("upregulated genes"), [])
        self.assertEqual(lint_text("negative controls"), [])

    def test_nonsignificant_equivalence_ignores_its_own_negation(self):
        hits = {item.rule_id for item in lint_text("The result was not significant and therefore equivalent.")}
        self.assertIn("nonsignificant_as_absence", hits)

    def test_invalid_rules_raise(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            duplicate = root / "dup.yaml"
            duplicate.write_text(
                "rules:\n"
                "  - {id: a, category: c, pattern: 'x', allowed: 'y', example: 'x'}\n"
                "  - {id: a, category: c, pattern: 'z', allowed: 'y', example: 'z'}\n",
                encoding="utf-8",
            )
            with self.assertRaises(ValueError) as caught:
                load_rules(duplicate)
            self.assertIn("a", str(caught.exception))
            broken = root / "bad.yaml"
            broken.write_text(
                "rules:\n  - {id: bad, category: c, pattern: '(', allowed: 'y', example: 'x'}\n",
                encoding="utf-8",
            )
            with self.assertRaises(ValueError) as caught:
                load_rules(broken)
            self.assertIn("bad", str(caught.exception))

    def test_cli_exits_and_refuses_a_second_json_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            clean = root / "clean.md"
            clean.write_text("Donor-unit association inside one arm.\n", encoding="utf-8")
            dirty = root / "dirty.md"
            dirty.write_text("The association was independently replicated.\n", encoding="utf-8")
            missing = root / "absent.md"
            self.assertEqual(claim_lint.main([str(clean)]), 0)
            self.assertEqual(claim_lint.main([str(dirty)]), 0)
            self.assertEqual(claim_lint.main([str(dirty), "--strict"]), EXIT_GATE)
            out = root / "findings.json"
            self.assertEqual(claim_lint.main([str(dirty), "--json", str(out)]), 0)
            self.assertEqual(json.loads(out.read_text(encoding="utf-8"))["n_findings"], 1)
            with self.assertRaises(SystemExit):
                claim_lint.main([str(dirty), "--json", str(out)])
            with self.assertRaises(SystemExit) as caught:
                claim_lint.main([str(missing)])
            self.assertEqual(caught.exception.code, EXIT_USAGE)


class StageReportTests(unittest.TestCase):
    def _part5(self, root: Path) -> Path:
        stage = root / "arm_a"
        (stage / "05_logs").mkdir(parents=True)
        (stage / "02_tables").mkdir()
        verdict = {
            "verdict": "SINGLE_SOURCE_DEPENDENT",
            "evidence_ceiling": "exploratory",
            "source_dependent": True,
            "n_donors": 14,
            "n_source_blocks": 3,
            "reason": "Dominant source holdout was not estimable.",
            "interpretation_limit": "Read the effect and the interval.",
            "reproduction_status": "NOT_APPLICABLE",
            "calibration_status": "NOT_ESTABLISHED_FOR_THIS_ANALYSIS",
            "scientifically_calibrated": False,
        }
        (stage / "05_logs" / "verdict.json").write_text(json.dumps(verdict), encoding="utf-8")
        audit = {
            "precision_status": "DESCRIPTIVE_CI_PRECISION_ONLY",
            "precision_summary": {
                "median_ci_half_width": 0.2,
                "p90_ci_half_width": 0.4,
                "effect_floor": 0.1,
                "fraction_ci_half_width_above_effect_floor": 0.5,
                "interpretation": "Marginal coefficient CI precision only.",
            },
        }
        (stage / "05_logs" / "model_audit.json").write_text(json.dumps(audit), encoding="utf-8")
        pd.DataFrame(
            {
                "gene": ["GENE_A", "GENE_B"],
                "logFC": [0.123456, 0.01],
                "q_bh": [0.2, 0.01],
                "robust_primary": ["FALSE", "TRUE"],
            }
        ).to_csv(stage / "02_tables" / "gene_evidence.csv", index=False)
        pd.DataFrame(
            {"subset": ["DROP_A"], "held_out_block": ["A"], "status": ["NOT_ESTIMABLE"], "rdf": [3]}
        ).to_csv(stage / "02_tables" / "fold_audit.csv", index=False)
        return stage

    def test_part5_report_keeps_the_token_and_lints_clean(self):
        with tempfile.TemporaryDirectory() as tmp:
            stage = self._part5(Path(tmp))
            self.assertEqual(stage_report.main([str(stage)]), 0)
            text = (stage / "06_reports" / "PART5_REPORT.md").read_text(encoding="utf-8")
            headings = [line for line in text.splitlines() if line.startswith("## ")]
            self.assertEqual(
                headings,
                [
                    "## 1. Conclusion",
                    "## 2. Reproduction anchor",
                    "## 3. Eligibility and coverage",
                    "## 4. Primary numbers",
                    "## 5. Robustness",
                    "## 6. Wording consequences",
                    "## 7. Explicitly not shown or claimed",
                ],
            )
            after = text.split("## 1. Conclusion", 1)[1].splitlines()
            first = next(line for line in after if line.strip())
            self.assertTrue(first.startswith("Verdict: `SINGLE_SOURCE_DEPENDENT`."))
            self.assertIn("0.1235", text)
            self.assertIn("14 donor-units across 3 source blocks", text)
            self.assertIn("`Cross-dataset or cross-source consistency`", text)
            self.assertEqual(lint_text(text), [])

    def test_second_report_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            stage = self._part5(Path(tmp))
            stage_report.main([str(stage)])
            path = stage / "06_reports" / "PART5_REPORT.md"
            before = path.read_text(encoding="utf-8")
            with self.assertRaises(SystemExit):
                stage_report.main([str(stage)])
            self.assertEqual(path.read_text(encoding="utf-8"), before)

    def test_missing_verdict_exits_without_writing(self):
        with tempfile.TemporaryDirectory() as tmp:
            stage = Path(tmp) / "arm"
            stage.mkdir()
            with self.assertRaises(SystemExit) as caught:
                stage_report.main([str(stage)])
            self.assertEqual(caught.exception.code, EXIT_GATE)
            self.assertFalse((stage / "06_reports").exists())

    def test_missing_gene_table_is_reported_as_absent(self):
        with tempfile.TemporaryDirectory() as tmp:
            stage = self._part5(Path(tmp))
            (stage / "02_tables" / "gene_evidence.csv").unlink()
            self.assertEqual(stage_report.main([str(stage)]), 0)
            text = (stage / "06_reports" / "PART5_REPORT.md").read_text(encoding="utf-8")
            self.assertIn("`02_tables/gene_evidence.csv` not found; not reported.", text)

    def test_part6_report_names_the_unpaired_mirror(self):
        with tempfile.TemporaryDirectory() as tmp:
            stage = Path(tmp) / "arm6"
            (stage / "05_logs").mkdir(parents=True)
            verdict = {
                "verdict": "PASS_WITH_LIMITATIONS",
                "evidence_class": "exploratory_embedding_only",
                "evidence_ceiling": "exploratory",
                "tags": ["KO_OE_UNPAIRED"],
                "reason": "Run finished under the frozen family.",
            }
            (stage / "05_logs" / "verdict.json").write_text(json.dumps(verdict), encoding="utf-8")
            self.assertEqual(stage_report.main([str(stage)]), 0)
            text = (stage / "06_reports" / "PART6_REPORT.md").read_text(encoding="utf-8")
            self.assertIn("# Part 6 report:", text)
            self.assertIn("`A paired KO/OE mechanistic mirror`", text)


class StageStatusTests(unittest.TestCase):
    def test_protocol_only_points_at_config(self):
        with tempfile.TemporaryDirectory() as tmp:
            stage = Path(tmp) / "arm"
            manifest = stage / "00_protocol_manifest"
            manifest.mkdir(parents=True)
            (manifest / "PROTOCOL.md").write_text("frozen analysis\n", encoding="utf-8")
            out = stage / "status.json"
            self.assertEqual(stage_status.main([str(stage), "--json", str(out)]), 0)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["next"]["id"], "config")

    def test_later_artifacts_without_an_earlier_step_warn(self):
        with tempfile.TemporaryDirectory() as tmp:
            stage = Path(tmp) / "arm"
            logs = stage / "05_logs"
            logs.mkdir(parents=True)
            (logs / "verdict.json").write_text('{"verdict": "INCONCLUSIVE"}\n', encoding="utf-8")
            out = stage / "status.json"
            self.assertEqual(stage_status.main([str(stage), "--json", str(out)]), 0)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertTrue(any("Later artifacts exist without step" in item for item in payload["warnings"]))

    def test_frozen_part5_skips_pathways_and_points_at_the_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            stage = Path(tmp) / "arm"
            manifest = stage / "00_protocol_manifest"
            manifest.mkdir(parents=True)
            protocol = manifest / "PROTOCOL.md"
            protocol.write_text("frozen analysis\n", encoding="utf-8")
            (manifest / "analysis_config.yaml").write_text("part: 5\n", encoding="utf-8")
            check_chronology(protocol, stage, freeze=True)
            pseudobulk = stage / "03_pseudobulk"
            pseudobulk.mkdir()
            for name in ("counts.mtx", "metadata.csv", "genes.csv"):
                (pseudobulk / name).write_text("x\n", encoding="utf-8")
            tables = stage / "02_tables"
            tables.mkdir()
            (tables / "metadata_with_source_block.csv").write_text("a\n", encoding="utf-8")
            (tables / "eligibility.csv").write_text("a\n", encoding="utf-8")
            (tables / "gene_evidence.csv").write_text("gene\n", encoding="utf-8")
            logs = stage / "05_logs"
            logs.mkdir()
            (logs / "model_audit.json").write_text("{}\n", encoding="utf-8")
            (logs / "verdict.json").write_text('{"verdict": "INCONCLUSIVE"}\n', encoding="utf-8")
            out = stage / "status.json"
            self.assertEqual(stage_status.main([str(stage), "--json", str(out)]), 0)
            payload = json.loads(out.read_text(encoding="utf-8"))
            status = {item["id"]: item["status"] for item in payload["steps"]}
            self.assertEqual(status["pathways"], "SKIPPED")
            self.assertEqual(status["report"], "MISSING")
            for step_id, step_status in status.items():
                if step_id in {"pathways", "report"}:
                    continue
                self.assertEqual(step_status, "DONE", step_id)
            self.assertEqual(payload["next"]["id"], "report")
            self.assertEqual(
                payload["chronology"],
                "LOCAL_SEQUENCE_CONSISTENT_NOT_TRUSTED_PREREGISTRATION",
            )

    def test_input_audit_directory_selects_part6(self):
        with tempfile.TemporaryDirectory() as tmp:
            stage = Path(tmp) / "arm"
            (stage / "00_input_audit").mkdir(parents=True)
            out = stage / "status.json"
            self.assertEqual(stage_status.main([str(stage), "--json", str(out)]), 0)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["part"], 6)

    def test_missing_directory_returns_1(self):
        with tempfile.TemporaryDirectory() as tmp:
            missing = Path(tmp) / "absent"
            self.assertEqual(stage_status.main([str(missing)]), 1)

    def test_json_has_steps_next_and_warnings(self):
        with tempfile.TemporaryDirectory() as tmp:
            stage = Path(tmp) / "arm"
            stage.mkdir()
            out = stage / "status.json"
            self.assertEqual(stage_status.main([str(stage), "--json", str(out)]), 0)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(set(payload) >= {"steps", "next", "warnings"}, True)


if __name__ == "__main__":
    unittest.main()
