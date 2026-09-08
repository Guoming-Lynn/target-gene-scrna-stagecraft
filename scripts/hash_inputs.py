#!/usr/bin/env python3
"""SHA-256 every existing path listed in a YAML or JSON manifest.

Usage:
    python scripts/hash_inputs.py path/to/analysis_config.yaml
    python scripts/hash_inputs.py path/to/input_manifest.json

Prints path<TAB>sha256<TAB>missing|ok. Exits 1 if any listed existing-key path is missing.
Looks for keys named path, paths, input, inputs, or values that look like files.
"""
from __future__ import annotations

import hashlib
import argparse
import json
import re
import sys
from pathlib import Path

try:
    import yaml  # type: ignore
except ImportError:
    yaml = None


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def collect(obj, out: list[str]) -> None:
    if isinstance(obj, dict):
        for k, v in obj.items():
            if str(k).lower() in {"path", "file", "h5ad", "counts", "metadata", "input", "inputs", "paths"} and isinstance(v, str):
                out.append(v)
            elif str(k).lower() in {"inputs", "paths"} and isinstance(v, list):
                for entry in v:
                    if isinstance(entry, str):
                        out.append(entry)
                    else:
                        collect(entry, out)
            else:
                collect(v, out)
    elif isinstance(obj, list):
        for x in obj:
            collect(x, out)
    elif isinstance(obj, str):
        p = obj.replace("\\", "/")
        if "/" in p and any(p.endswith(ext) for ext in (".h5ad", ".csv", ".tsv", ".mtx", ".json", ".yaml", ".yml", ".gz", ".txt", ".R", ".py")):
            out.append(obj)


def load(path: Path):
    raw = path.read_text(encoding="utf-8")
    if path.suffix.lower() in {".json"}:
        return json.loads(raw)
    if yaml is None:
        raise SystemExit("PyYAML required for yaml manifests: pip install pyyaml")
    return yaml.safe_load(raw)


def hashed_entries(obj):
    """Read explicit {path, sha256} and sibling <name>_sha256 contracts."""
    if isinstance(obj, dict):
        for key in ("path", "file", "h5ad"):
            if isinstance(obj.get(key), str) and "sha256" in obj:
                yield obj[key], obj["sha256"]
                break
        for key, value in obj.items():
            if isinstance(key, str) and key.endswith("_sha256"):
                path_key = key[:-7]
                if not isinstance(obj.get(path_key), str):
                    raise ValueError(f"{key} needs a sibling {path_key} path")
                yield obj[path_key], value
            yield from hashed_entries(value)
    elif isinstance(obj, list):
        for value in obj:
            yield from hashed_entries(value)


def manifest_root(manifest: Path) -> Path:
    parent = manifest.resolve().parent
    return parent.parent if parent.name == "00_protocol_manifest" else parent


def verify_entries(obj, root: Path) -> list[dict]:
    entries = list(hashed_entries(obj))
    if not entries:
        raise ValueError("No explicit SHA-256 entries; freeze an input manifest first.")
    results = []
    for name, expected in entries:
        if not isinstance(expected, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", expected):
            raise ValueError(f"Missing/invalid frozen SHA-256 for {name}")
        path = Path(name)
        if not path.is_absolute():
            path = root / path
        if not path.is_file():
            raise ValueError(f"Missing input: {path}")
        actual = sha256(path)
        if actual != expected.lower():
            raise ValueError(f"SHA-256 mismatch: {path}")
        results.append({"path": str(path.resolve()), "sha256": actual, "status": "verified"})
    return results


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--verify", action="store_true", help="Require and compare explicit frozen SHA-256 entries")
    args = parser.parse_args(argv[1:])
    manifest = args.manifest
    obj = load(manifest)
    if args.verify:
        try:
            results = verify_entries(obj, manifest_root(manifest))
        except ValueError as exc:
            print(str(exc), file=sys.stderr)
            return 1
        print(json.dumps(results, indent=2))
        return 0
    paths: list[str] = []
    collect(obj, paths)
    if not paths:
        print("No input paths found in manifest", file=sys.stderr)
        return 1
    # unique, stable order
    seen = []
    for p in paths:
        if p not in seen:
            seen.append(p)
    missing = 0
    root = manifest_root(manifest)
    for p in seen:
        path = Path(p)
        if not path.is_absolute():
            path = root / path
        if not path.is_file():
            cand = (root / p).resolve()
            if cand.is_file():
                path = cand
            else:
                print(f"{p}\t-\tmissing")
                missing += 1
                continue
        print(f"{path}\t{sha256(path)}\tok")
    return 1 if missing else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

