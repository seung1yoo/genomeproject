from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

import duckdb


def write_csv(rows: list[dict[str, object]], path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"Cannot write an empty table: {output}")
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output


def csv_to_parquet(csv_path: str | Path, parquet_path: str | Path) -> Path:
    source = Path(csv_path)
    output = Path(parquet_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with duckdb.connect() as connection:
        connection.execute(
            "COPY (SELECT * FROM read_csv_auto(? , header=true)) TO ? (FORMAT PARQUET)",
            [str(source), str(output)],
        )
    return output


def build_catalog(tables: Iterable[tuple[str, str | Path]], database_path: str | Path) -> Path:
    output = Path(database_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with duckdb.connect(str(output)) as connection:
        for table_name, source in tables:
            if not table_name.replace("_", "").isalnum():
                raise ValueError(f"Unsafe table name: {table_name}")
            connection.execute(
                f"CREATE OR REPLACE TABLE {table_name} AS SELECT * FROM read_csv_auto(?, header=true)",
                [str(source)],
            )
    return output
