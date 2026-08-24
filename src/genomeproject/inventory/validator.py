from __future__ import annotations

from collections import defaultdict

from .models import InventoryRow


def summarize_inventory(
    rows: list[InventoryRow], expected_data_types: tuple[str, ...]
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    samples = sorted({row.sample_id for row in rows})
    by_type: dict[str, list[InventoryRow]] = defaultdict(list)
    present: set[tuple[str, str]] = set()
    for row in rows:
        by_type[row.data_type].append(row)
        if row.exists and row.is_file:
            present.add((row.sample_id, row.data_type))

    summary: list[dict[str, object]] = []
    for data_type in expected_data_types:
        typed = by_type.get(data_type, [])
        unique_samples = {row.sample_id for row in typed if row.exists and row.is_file}
        summary.append(
            {
                "data_type": data_type,
                "manifest_rows": len(typed),
                "file_count": sum(row.exists and row.is_file for row in typed),
                "unique_sample_count": len(unique_samples),
                "expected_sample_count": len(samples),
                "missing_sample_count": len(samples) - len(unique_samples),
                "total_bytes": sum(row.size_bytes or 0 for row in typed),
                "unindexed_count": sum(row.index_exists is False for row in typed),
                "problem_row_count": sum(row.status != "OK" for row in typed),
            }
        )

    completeness: list[dict[str, object]] = []
    for sample_id in samples:
        item: dict[str, object] = {"sample_id": sample_id}
        for data_type in expected_data_types:
            item[data_type] = (sample_id, data_type) in present
        item["complete"] = all(bool(item[data_type]) for data_type in expected_data_types)
        completeness.append(item)
    return summary, completeness
