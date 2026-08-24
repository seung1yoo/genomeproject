from __future__ import annotations

from dataclasses import asdict, dataclass


ALLOWED_DATA_TYPES = frozenset({"cnv", "cram", "dtc", "gvcf", "sv", "vcf"})


@dataclass(frozen=True)
class ManifestRow:
    sample_id: str
    data_type: str
    path: str


@dataclass(frozen=True)
class InventoryRow:
    sample_id: str
    data_type: str
    path: str
    exists: bool
    is_file: bool
    size_bytes: int | None
    index_path: str | None
    index_exists: bool | None
    duplicate_path: bool
    duplicate_sample_type: bool
    status: str
    error: str | None

    def as_dict(self) -> dict[str, object]:
        return asdict(self)
