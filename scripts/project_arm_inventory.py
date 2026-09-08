"""Report every declared Part 5 arm, including missing or failed runs.

An inventory is a disclosure check, not a project-level multiplicity correction
or proof that no undeclared analysis was conducted elsewhere.
"""
import argparse
import hashlib
import json
from pathlib import Path

import yaml


def inventory(manifest_path):
    path = Path(manifest_path).resolve()
    manifest = yaml.safe_load(path.read_text(encoding="utf-8"))
    arms = manifest.get("arms") if isinstance(manifest, dict) else None
    if not isinstance(arms, list) or not arms:
        raise ValueError("Declare a nonempty arms list before running Part 5")
    ids, paths, records = set(), set(), []
    for arm in arms:
        arm_id = str(arm.get("id", "")).strip()
        relative = arm.get("verdict")
        if not arm_id or arm_id in ids or not isinstance(relative, str) or not relative.strip():
            raise ValueError("Every arm needs a unique id and a verdict path")
        target = (path.parent / relative).resolve()
        target.relative_to(path.parent)
        if target in paths:
            raise ValueError("Different arms cannot reuse a verdict file")
        ids.add(arm_id)
        paths.add(target)
        row = {"arm_id": arm_id, "verdict_path": relative, "status": "MISSING_RESULT",
               "verdict": None, "sha256": None}
        if target.is_file():
            data = target.read_bytes()
            result = json.loads(data)
            if not isinstance(result, dict) or not isinstance(result.get("status"), str):
                raise ValueError("Invalid verdict object: " + relative)
            token = result.get("verdict")
            if result["status"] == "SUCCESS" and (not isinstance(token, str) or not token.strip()):
                raise ValueError("Successful verdict file needs a verdict token")
            row.update(status=result["status"], verdict=token, sha256=hashlib.sha256(data).hexdigest())
        records.append(row)
    return {"status": "COMPLETE_DECLARED_INVENTORY" if all(r["status"] != "MISSING_RESULT" for r in records)
            else "INCOMPLETE_DECLARED_INVENTORY", "n_declared_arms": len(records),
            "n_results_present": sum(r["status"] != "MISSING_RESULT" for r in records),
            "n_frozen_pass": sum(str(r["verdict"]).startswith("FROZEN_PASS") for r in records),
            "project_multiplicity_control": "NOT_ESTABLISHED", "undeclared_arms": "NOT_DETECTABLE",
            "manifest_sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "arms": records}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = inventory(args.manifest)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2)
    print(report["status"])
    return 0 if report["status"] == "COMPLETE_DECLARED_INVENTORY" else 2


if __name__ == "__main__":
    raise SystemExit(main())
