"""Shared text patterns used by more than one language."""

from __future__ import annotations

from importlib.resources import files


def unlikely_arm_pattern() -> str:
    """Labels that must not enter a formal Part 5 arm."""
    return files("stagecraft").joinpath("unlikely_arm_pattern.txt").read_text(encoding="utf-8").strip()
