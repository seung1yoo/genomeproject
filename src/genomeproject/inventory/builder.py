from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

from .models import ALLOWED_DATA_TYPES, InventoryRow, ManifestRow


REQUIRED_COLUMNS = ("sample_id", "data_type", "path")


def read_manifest(path: str | Path) -> list[ManifestRow]:
    manifest_path = Path(path)
    with manifest_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [name for name in REQUIRED_COLUMNS if name not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"Manifest is missing columns: {', '.join(missing)}")
        rows: list[ManifestRow] = []
        for line_no, raw in enumerate(reader, start=2):
            sample_id = (raw.get("sample_id") or "").strip()
            data_type = (raw.get("data_type") or "").strip().lower()
            file_path = (raw.get("path") or "").strip()
            if not sample_id or not file_path:
                raise ValueError(f"Empty sample_id/path at line {line_no}")
            if data_type not in ALLOWED_DATA_TYPES:
                raise ValueError(f"Unsupported data_type '{data_type}' at line {line_no}")
            rows.append(ManifestRow(sample_id, data_type, file_path))
    if not rows:
        raise ValueError("Manifest contains no data rows")
    return rows


def candidate_index_paths(file_path: Path, data_type: str) -> list[Path]:
    text = str(file_path)
    if data_type == "cram" or text.lower().endswith(".cram"):
        return [Path(text + ".crai"), file_path.with_suffix(".crai")]
    if data_type in {"vcf", "gvcf", "sv"} and text.lower().endswith(
        (".vcf.gz", ".g.vcf.gz", ".gvcf.gz")
    ):
        return [Path(text + ".tbi"), Path(text + ".csi")]
    if text.lower().endswith(".bcf"):
        return [Path(text + ".csi")]
    return []


def build_inventory(rows: list[ManifestRow]) -> list[InventoryRow]:
    normalized_paths = [str(Path(row.path).expanduser().absolute()) for row in rows]
    path_counts = Counter(normalized_paths)
    pair_counts = Counter((row.sample_id, row.data_type) for row in rows)
    result: list[InventoryRow] = []
    for row, normalized_path in zip(rows, normalized_paths):
        file_path = Path(normalized_path)
        exists = file_path.exists()
        is_file = file_path.is_file() if exists else False
        size: int | None = None
        error: str | None = None
        try:
            if is_file:
                size = file_path.stat().st_size
        except OSError as exc:
            error = str(exc)

        candidates = candidate_index_paths(file_path, row.data_type)
        existing_index = next((item for item in candidates if item.is_file()), None)
        index_path = str(existing_index or candidates[0]) if candidates else None
        index_exists = bool(existing_index) if candidates else None
        duplicate_path = path_counts[normalized_path] > 1
        duplicate_pair = pair_counts[(row.sample_id, row.data_type)] > 1

        problems: list[str] = []
        if not exists:
            problems.append("MISSING")
        elif not is_file:
            problems.append("NOT_FILE")
        if candidates and not index_exists:
            problems.append("UNINDEXED")
        if duplicate_path:
            problems.append("DUPLICATE_PATH")
        if duplicate_pair:
            problems.append("DUPLICATE_SAMPLE_TYPE")
        if error:
            problems.append("STAT_ERROR")
        result.append(
            InventoryRow(
                sample_id=row.sample_id,
                data_type=row.data_type,
                path=normalized_path,
                exists=exists,
                is_file=is_file,
                size_bytes=size,
                index_path=index_path,
                index_exists=index_exists,
                duplicate_path=duplicate_path,
                duplicate_sample_type=duplicate_pair,
                status=";".join(problems) if problems else "OK",
                error=error,
            )
        )
    return result
