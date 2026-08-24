from pathlib import Path

import duckdb

from genomeproject.storage.exporters import csv_to_parquet


def test_csv_to_parquet_preserves_rows(tmp_path: Path) -> None:
    source = tmp_path / "inventory.csv"
    source.write_text("sample_id,data_type\nS001,cram\nS002,gvcf\n", encoding="utf-8")
    output = tmp_path / "nested" / "inventory.parquet"

    result = csv_to_parquet(source, output)

    assert result == output
    assert output.is_file()
    with duckdb.connect() as connection:
        rows = connection.execute(
            "SELECT sample_id, data_type FROM read_parquet(?) ORDER BY sample_id",
            [str(output)],
        ).fetchall()
    assert rows == [("S001", "cram"), ("S002", "gvcf")]
