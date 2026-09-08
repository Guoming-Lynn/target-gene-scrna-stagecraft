"""Run helper smoke tests; does not execute a formal R or Geneformer analysis."""
import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--timeout', type=float, default=300)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    commands = [
        ('environment', ['scripts/check_environment.py', '--stage', 'smoke']),
        ('toy', ['scripts/generate_toy_data.py', '--out', str(args.out.resolve() / 'toy.h5ad')]),
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
    return 0 if passed else 2

if __name__ == '__main__':
    raise SystemExit(main())
