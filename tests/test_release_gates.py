import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from check_environment import main as environment_main
from part5_verdict import guarded_verdict
from part6_verdict import main as verdict_main

class ReleaseGates(unittest.TestCase):
    def test_missing_r_is_fatal_for_part5(self):
        with patch('check_environment.check', return_value={'status': 'PASS'}), patch('check_environment.shutil.which', return_value=None):
            self.assertEqual(environment_main(['--stage', 'part5', '--rscript', '']), 2)

    def test_broken_r_is_fatal(self):
        def probe(command, *args):
            return {'status': 'FAILED' if command[0] == 'bad-r' else 'PASS'}
        with patch('check_environment.check', side_effect=probe):
            self.assertEqual(environment_main(['--stage', 'part5', '--rscript', 'bad-r']), 2)

    def test_source_dependence_overrides_custom_pass(self):
        result = guarded_verdict([{'token': 'FROZEN_PASS'}], {'source_dependent': True})
        self.assertEqual(result['verdict'], 'SINGLE_SOURCE_DEPENDENT')

    def test_missing_lodo_cannot_pass(self):
        result = guarded_verdict([{'token': 'FROZEN_PASS'}], {'source_dependent': False, 'formal_gates_pass': True})
        self.assertEqual(result['verdict'], 'INCONCLUSIVE')

    def test_part6_ceiling_cannot_be_upgraded(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            (root / 'audit.json').write_text(json.dumps({'evidence_ceiling': 'causal', 'can_only_downgrade': False}))
            (root / 'table.yaml').write_text('rows:\n  - token: PASS\n    when: {}\n')
            verdict_main([str(root / 'audit.json'), '--table', str(root / 'table.yaml'), '--out', str(root / 'out.json')])
            result = json.loads((root / 'out.json').read_text())
            self.assertEqual(result['evidence_ceiling'], 'exploratory_embedding_only')
            self.assertTrue(result['ko_oe_unpaired'])
            self.assertFalse(result['causal_inference'])
            self.assertEqual(result['perturbation_biological_validity'], 'NOT_ESTABLISHED_BY_EMBEDDING_SHIFT')
