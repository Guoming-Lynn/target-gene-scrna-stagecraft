"""Local freeze receipts detect accidental amendments; not trusted preregistration."""
import hashlib
import json
import time
from pathlib import Path


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_chronology(protocol, stage_root, freeze=False):
    protocol, root = Path(protocol).resolve(), Path(stage_root).resolve()
    manifest = root / "00_protocol_manifest"
    receipt = manifest / "protocol_freeze.json"
    protocol.relative_to(root)
    outputs = [p for folder in ("03_pseudobulk", "02_tables")
               for p in (root / folder).rglob("*") if p.is_file()]
    if freeze:
        if receipt.exists():
            raise ValueError("Freeze receipt already exists; record an amendment in a new stage")
        if outputs:
            raise ValueError("Cannot claim a new pre-output freeze: stage artifacts already exist")
        manifest.mkdir(parents=True, exist_ok=True)
        inputs = set(p for p in manifest.rglob("*") if p.is_file()) | {protocol}
        payload = {"frozen_at_ns": time.time_ns(), "evidence_class": "LOCAL_CLOCK_ONLY",
                   "files": {p.relative_to(root).as_posix(): sha256(p) for p in sorted(inputs)}}
        with receipt.open("x", encoding="utf-8") as stream:
            json.dump(payload, stream, indent=2)
        return "LOCAL_FREEZE_RECORDED_NOT_TRUSTED_PREREGISTRATION"
    if not receipt.exists():
        raise ValueError("No freeze receipt; timing is NOT_VERIFIED. Do not backdate existing results")
    payload = json.loads(receipt.read_text(encoding="utf-8"))
    entries = payload["files"]
    if protocol.relative_to(root).as_posix() not in entries:
        raise ValueError("Protocol not included in receipt")
    current = {p.relative_to(root).as_posix() for p in manifest.rglob("*")
               if p.is_file() and p != receipt} | {protocol.relative_to(root).as_posix()}
    if current != set(entries):
        raise ValueError("Manifest file set changed after freeze")
    for name, expected in entries.items():
        path = (root / name).resolve()
        path.relative_to(root)
        if not path.is_file() or sha256(path) != expected:
            raise ValueError("Frozen input changed: " + name)
    frozen_at = payload["frozen_at_ns"]
    if frozen_at > time.time_ns():
        raise ValueError("Freeze timestamp is in the future")
    if any(p.stat().st_mtime_ns < frozen_at for p in outputs):
        raise ValueError("Artifact predates local freeze; chronology needs manual review")
    return "LOCAL_SEQUENCE_CONSISTENT_NOT_TRUSTED_PREREGISTRATION"
