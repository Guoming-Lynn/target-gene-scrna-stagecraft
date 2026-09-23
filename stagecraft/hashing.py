"""SHA-256 helpers used by provenance writers."""

from __future__ import annotations

import hashlib
from pathlib import Path

PathLike = str | Path


def sha256_file(path: PathLike, block_size: int = 2**20) -> str:
    if isinstance(block_size, bool) or not isinstance(block_size, int) or block_size <= 0:
        raise ValueError("block_size must be a positive integer")
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()
