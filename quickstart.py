"""Run helper smoke tests; does not execute a formal R or Geneformer analysis."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent

from stagecraft import EXIT_GATE, EXIT_OK  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--timeout', type=float, default=300)
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=False)
    commands = [
        ('environment', ['scripts/check_environment.py', '--stage', 'smoke']),
        ('toy_demo', ['scripts/demo_toy_run.py', '--out', str(args.out.resolve() / 'toy_demo')]),
        ('part5_helpers', ['scripts/_part5_smoke.py']),
        ('part6_helpers', ['scripts/_part6_smoke.py']),
    ]
    results = []
    for name, command in commands:
        print('RUNNING: ' + name, flush=True)
        try:
            code = subprocess.run([sys.executable, '-u'] + command, cwd=ROOT, timeout=args.timeout).returncode
            status = 'PASS' if code == 0 else 'FAILED'
        except subprocess.TimeoutExpired:
            code, status = 2, 'TIMEOUT'
        except OSError as exc:
            print(str(exc), flush=True)
            code, status = 2, 'LAUNCH_FAILED'
        results.append(dict(component=name, status=status, exit_code=code))
        print(f'{name}: {status} (exit {code})', flush=True)
        if code:
            break
    passed = len(results) == len(commands) and all(r['status'] == 'PASS' for r in results)
    report = dict(status='PASS' if passed else 'FAILED', formal_analysis=False, results=results)
    (args.out / 'quickstart_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    return EXIT_OK if passed else EXIT_GATE

if __name__ == '__main__':
    raise SystemExit(main())
