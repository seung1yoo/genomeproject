from pathlib import Path

import pytest

from genomeproject.inventory.builder import build_inventory, read_manifest
from genomeproject.inventory.validator import summarize_inventory


def test_inventory_detects_missing_index_duplicate_and_completeness(tmp_path: Path) -> None:
    cram = tmp_path / "S1.cram"
    cram.write_bytes(b"cram")
    (tmp_path / "S1.cram.crai").write_bytes(b"index")
    gvcf = tmp_path / "S1.g.vcf.gz"
    gvcf.write_bytes(b"gvcf")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "sample_id\tdata_type\tpath\n"
        f"S1\tcram\t{cram}\n"
        f"S1\tgvcf\t{gvcf}\n"
        f"S2\tgvcf\t{tmp_path / 'missing.g.vcf.gz'}\n",
        encoding="utf-8",
    )

    rows = build_inventory(read_manifest(manifest))
    assert rows[0].status == "OK"
    assert rows[1].status == "UNINDEXED"
    assert "MISSING" in rows[2].status
    summary, completeness = summarize_inventory(rows, ("cram", "gvcf"))
    assert summary[0]["file_count"] == 1
    assert summary[1]["unindexed_count"] == 2
    assert completeness == [
        {"sample_id": "S1", "cram": True, "gvcf": True, "complete": True},
        {"sample_id": "S2", "cram": False, "gvcf": False, "complete": False},
    ]


def test_manifest_rejects_duplicate_unsupported_type(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tdata_type\tpath\nS1\tbam\t/a.bam\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Unsupported data_type"):
        read_manifest(manifest)
