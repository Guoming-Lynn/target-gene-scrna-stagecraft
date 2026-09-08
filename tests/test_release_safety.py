import hashlib
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import package_skill
import part6_eligibility
import part6_endpoints
from part5_verdict import load_table
from part6_verdict import guarded_verdict


class ReleaseSafety(unittest.TestCase):
    def audit(self):
        return dict(blocker=False, smoke_pass=True, sham_pass=True, run_finished=True,
                    ko_primary_bh_pass=True, ko_eligible_cells=50,
                    n_sign_tests_run=2, family_size=2, n_eligible_donors_ko=10)

    def test_part6_requires_complete_typed_evidence(self):
        rows = load_table(ROOT / 'scripts/part6_verdict_table.example.yaml')
        self.assertEqual(guarded_verdict(rows, self.audit())['verdict'], 'EMBEDDING_SHIFT_CONSISTENT')
        for key in self.audit():
            for value in (None, 'true', -1):
                audit = self.audit()
                audit[key] = value
                with self.subTest(key=key, value=value):
                    self.assertEqual(guarded_verdict(rows, audit)['verdict'], 'INCONCLUSIVE')
        for audit in ({'ko_primary_bh_pass': True}, {'run_finished': True}):
            self.assertEqual(guarded_verdict(rows, audit)['verdict'], 'INCONCLUSIVE')

    def test_custom_table_cannot_bypass_gates(self):
        rows = [{'token': 'EMBEDDING_SHIFT_CONSISTENT'}]
        for key, value, expected in (
            ('blocker', True, 'STOPPED'), ('smoke_pass', False, 'STOPPED'),
            ('sham_pass', False, 'SHAM_DRIFT'), ('run_finished', False, 'INCONCLUSIVE'),
            ('ko_eligible_cells', 0, 'TOKEN_UNOBSERVABLE'),
            ('n_eligible_donors_ko', 4, 'NOT_ESTIMABLE'),
            ('n_sign_tests_run', 0, 'NOT_ESTIMABLE'),
            ('family_size', 1, 'INCONCLUSIVE'), ('ko_primary_bh_pass', False, 'INCONCLUSIVE')):
            audit = self.audit()
            audit[key] = value
            self.assertEqual(guarded_verdict(rows, audit)['verdict'], expected)
        audit = self.audit()
        audit['ko_primary_bh_pass'] = False
        rows = load_table(ROOT / 'scripts/part6_verdict_table.example.yaml')
        self.assertEqual(guarded_verdict(rows, audit)['verdict'], 'PASS_WITH_LIMITATIONS')

    def test_eligibility_preserves_both_outputs_and_rejects_alias(self):
        frame = pd.DataFrame({'example': [1]})
        for existing in ('ledger.csv', 'donor_effects_eligible.csv', None):
            with tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                out = root / ('ledger.csv' if existing else 'donor_effects_eligible.csv')
                if existing:
                    (root / existing).write_bytes(b'frozen')
                with patch.object(part6_eligibility.pd, 'read_csv', return_value=frame), patch.object(
                    part6_eligibility, 'donor_eligibility', return_value=(frame, frame)
                ), self.assertRaises(SystemExit):
                    part6_eligibility.main(['input.csv', '--out', str(out)])
                if existing:
                    self.assertEqual((root / existing).read_bytes(), b'frozen')
                self.assertEqual(len(list(root.iterdir())), int(existing is not None))

    def test_endpoint_members_preserved_before_any_write(self):
        for existing in ('coverage.csv', 'coverage_members.csv'):
            with tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                (root / 'sets.yaml').write_text('{}')
                (root / 'genes.txt').write_text('A')
                (root / existing).write_bytes(b'frozen')
                with patch.object(part6_endpoints, 'evaluate_sets', return_value=pd.DataFrame()), self.assertRaises(SystemExit):
                    part6_endpoints.main(['--sets', str(root/'sets.yaml'), '--model-genes', str(root/'genes.txt'),
                                          '--target', 'G', '--out', str(root/'coverage.csv')])
                self.assertEqual((root / existing).read_bytes(), b'frozen')
                self.assertEqual(len(list(root.iterdir())), 3)

    def test_archive_contains_only_reviewed_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / 'skill'
            root.mkdir()
            (root / 'SKILL.md').write_text('---\nmetadata: {version: "1.0"}\n---\n')
            (root / 'release-files.txt').write_text('SKILL.md\nrelease-files.txt\n')
            for name in ('.env', 'patients.csv', 'atlas.h5ad', 'weights.safetensors', 'unreviewed.py'):
                (root / name).write_bytes(b'private')
            out = package_skill.package(root)
            with zipfile.ZipFile(out) as archive:
                self.assertEqual(set(archive.namelist()), {'skill/SKILL.md', 'skill/release-files.txt'})
                self.assertIsNone(archive.testzip())
            self.assertIn(hashlib.sha256(out.read_bytes()).hexdigest(), out.with_suffix('.sha256.txt').read_text())
            with self.assertRaises(SystemExit):
                package_skill.package(root)
            for text in ('../outside.txt\n', 'missing.py\n', 'SKILL.md\nSKILL.md\n', ''):
                (root / 'release-files.txt').write_text(text)
                with self.assertRaises(SystemExit):
                    package_skill.release_files(root)

    def test_repository_release_manifest_resolves(self):
        files = package_skill.release_files(ROOT)
        self.assertIn(ROOT / 'LICENSE', files)
        self.assertIn(ROOT / 'scripts/part6_verdict.py', files)

    def test_packager_accepts_an_explicit_output_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / 'skill'
            root.mkdir()
            (root / 'SKILL.md').write_text('---\nmetadata: {version: "1.0"}\n---\n')
            (root / 'release-files.txt').write_text('SKILL.md\nrelease-files.txt\n')
            output = Path(tmp) / 'release'
            output.mkdir()
            with patch.object(package_skill, '__file__', str(root / 'scripts' / 'package_skill.py')):
                # The library entry point keeps repository selection explicit;
                # this verifies callers can choose a writable destination.
                archive = package_skill.package(root, output)
            self.assertEqual(archive.parent, output)
