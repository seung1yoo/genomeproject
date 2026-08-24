from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class AppConfig:
    reference_build: str
    reference_fasta: Path
    expected_data_types: tuple[str, ...]
    min_gq: int
    min_dp: int
    autosomes_only: bool
    normalize_markers: bool
    workers: int
    output_root: Path


def _require(mapping: dict[str, Any], key: str, section: str) -> Any:
    if key not in mapping:
        raise ValueError(f"Missing config key: {section}.{key}")
    return mapping[key]


def load_config(path: str | Path) -> AppConfig:
    config_path = Path(path)
    with config_path.open(encoding="utf-8") as handle:
        raw = yaml.safe_load(handle) or {}
    reference = raw.get("reference", {})
    inventory = raw.get("inventory", {})
    markers = raw.get("markers", {})
    output = raw.get("output", {})
    return AppConfig(
        reference_build=str(_require(reference, "build", "reference")),
        reference_fasta=Path(_require(reference, "fasta", "reference")),
        expected_data_types=tuple(
            inventory.get("expected_data_types", ["cnv", "cram", "dtc", "gvcf", "sv", "vcf"])
        ),
        min_gq=int(markers.get("min_gq", 20)),
        min_dp=int(markers.get("min_dp", 8)),
        autosomes_only=bool(markers.get("autosomes_only", True)),
        normalize_markers=bool(markers.get("normalize", True)),
        workers=max(1, int(markers.get("workers", 4))),
        output_root=Path(output.get("root", "results")),
    )
