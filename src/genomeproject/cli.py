from __future__ import annotations

import argparse
import json
import shutil
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import yaml

from .config import AppConfig, load_config
from .inventory.builder import build_inventory, read_manifest
from .inventory.validator import summarize_inventory
from .markers.frequency import aggregate_frequencies
from .markers.gvcf_reader import read_sample_gvcf
from .markers.models import GenotypeObservation, Marker
from .markers.normalizer import normalize_markers, read_markers, validate_reference
from .storage.exporters import build_catalog, csv_to_parquet, write_csv


def _run_id(value: str | None) -> str:
    return value or datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def _resolved_config(config: AppConfig, manifest: Path, markers: Path | None = None) -> dict[str, object]:
    payload: dict[str, object] = {
        "reference": {"build": config.reference_build, "fasta": str(config.reference_fasta.absolute())},
        "inventory": {"expected_data_types": list(config.expected_data_types)},
        "markers": {
            "min_gq": config.min_gq,
            "min_dp": config.min_dp,
            "autosomes_only": config.autosomes_only,
            "normalize": config.normalize_markers,
            "workers": config.workers,
        },
        "inputs": {"manifest": str(manifest.absolute())},
    }
    if markers:
        payload["inputs"]["markers"] = str(markers.absolute())  # type: ignore[index]
    return payload


def _write_run_config(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False)


def command_inventory_build(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    manifest_path = Path(args.manifest)
    output_dir = Path(args.output) if args.output else config.output_root / _run_id(args.run_id) / "inventory"
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = build_inventory(read_manifest(manifest_path))
    summary, completeness = summarize_inventory(rows, config.expected_data_types)

    inventory_csv = write_csv([row.as_dict() for row in rows], output_dir / "file_inventory.csv")
    summary_csv = write_csv(summary, output_dir / "inventory_summary.csv")
    completeness_csv = write_csv(completeness, output_dir / "sample_completeness.csv")
    for source in (inventory_csv, summary_csv, completeness_csv):
        csv_to_parquet(source, source.with_suffix(".parquet"))
    catalog = build_catalog(
        [
            ("file_inventory", inventory_csv),
            ("inventory_summary", summary_csv),
            ("sample_completeness", completeness_csv),
        ],
        output_dir / "inventory.duckdb",
    )
    _write_run_config(output_dir.parent / "run_config.yaml", _resolved_config(config, manifest_path))
    print(json.dumps({"rows": len(rows), "samples": len(completeness), "catalog": str(catalog)}, ensure_ascii=False))
    return 0


def command_inventory_summary(args: argparse.Namespace) -> int:
    catalog = Path(args.catalog)
    if not catalog.is_file():
        raise FileNotFoundError(f"Catalog not found: {catalog}")
    with duckdb.connect(str(catalog), read_only=True) as connection:
        columns = [item[0] for item in connection.execute("DESCRIBE inventory_summary").fetchall()]
        rows = connection.execute("SELECT * FROM inventory_summary ORDER BY data_type").fetchall()
    print("\t".join(columns))
    for row in rows:
        print("\t".join("" if value is None else str(value) for value in row))
    return 0


def _prepare_markers(args: argparse.Namespace, config: AppConfig) -> list[Marker]:
    markers = read_markers(args.markers, autosomes_only=config.autosomes_only)
    validate_reference(markers, config.reference_fasta)
    should_normalize = config.normalize_markers and not getattr(args, "skip_normalization", False)
    if should_normalize:
        markers = normalize_markers(markers, config.reference_fasta)
        validate_reference(markers, config.reference_fasta)
    return markers


def command_markers_validate(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    markers = _prepare_markers(args, config)
    print(json.dumps({"markers": len(markers), "build": config.reference_build}, ensure_ascii=False))
    return 0


def _read_worker(payload: tuple[str, str, list[Marker], int, int]) -> list[GenotypeObservation]:
    return read_sample_gvcf(*payload)


def command_markers_frequency(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    manifest_path = Path(args.manifest)
    marker_path = Path(args.markers)
    markers = _prepare_markers(args, config)
    manifest_rows = [row for row in read_manifest(manifest_path) if row.data_type == "gvcf"]
    if not manifest_rows:
        raise ValueError("Manifest contains no gvcf rows")
    pairs = [(row.sample_id, row.data_type) for row in manifest_rows]
    if len(pairs) != len(set(pairs)):
        raise ValueError("Each sample must have exactly one gvcf row")

    run_id = _run_id(args.run_id)
    output_dir = Path(args.output) if args.output else config.output_root / run_id / "marker_frequency"
    output_dir.mkdir(parents=True, exist_ok=True)
    observations: list[GenotypeObservation] = []
    payloads = [
        (row.sample_id, row.path, markers, config.min_gq, config.min_dp) for row in manifest_rows
    ]
    with ProcessPoolExecutor(max_workers=config.workers) as executor:
        futures = [executor.submit(_read_worker, payload) for payload in payloads]
        for future in as_completed(futures):
            observations.extend(future.result())
    observations.sort(key=lambda item: (item.chrom, item.pos, item.ref, item.alt, item.sample_id))

    audit_dir = output_dir / "audit_by_chrom"
    for chrom in sorted({item.chrom for item in observations}):
        chrom_rows = [item.as_dict() for item in observations if item.chrom == chrom]
        source = write_csv(chrom_rows, audit_dir / f"chrom={chrom}" / "observations.csv")
        csv_to_parquet(source, source.with_suffix(".parquet"))
    frequency_csv = write_csv(aggregate_frequencies(markers, observations), output_dir / "marker_frequency.csv")
    csv_to_parquet(frequency_csv, frequency_csv.with_suffix(".parquet"))
    _write_run_config(
        output_dir.parent / "run_config.yaml", _resolved_config(config, manifest_path, marker_path)
    )
    print(
        json.dumps(
            {"run_id": run_id, "samples": len(manifest_rows), "markers": len(markers), "output": str(output_dir)},
            ensure_ascii=False,
        )
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="genomeproject", description="WGS inventory and marker-frequency EDA")
    groups = parser.add_subparsers(dest="group", required=True)

    inventory = groups.add_parser("inventory")
    inventory_commands = inventory.add_subparsers(dest="command", required=True)
    inventory_build = inventory_commands.add_parser("build")
    inventory_build.add_argument("--manifest", required=True)
    inventory_build.add_argument("--config", default="configs/default.yaml")
    inventory_build.add_argument("--output")
    inventory_build.add_argument("--run-id")
    inventory_build.set_defaults(handler=command_inventory_build)
    inventory_summary = inventory_commands.add_parser("summary")
    inventory_summary.add_argument("--catalog", required=True)
    inventory_summary.set_defaults(handler=command_inventory_summary)

    marker_group = groups.add_parser("markers")
    marker_commands = marker_group.add_subparsers(dest="command", required=True)
    marker_validate = marker_commands.add_parser("validate")
    marker_validate.add_argument("--markers", required=True)
    marker_validate.add_argument("--config", default="configs/default.yaml")
    marker_validate.add_argument("--skip-normalization", action="store_true")
    marker_validate.set_defaults(handler=command_markers_validate)
    marker_frequency = marker_commands.add_parser("frequency")
    marker_frequency.add_argument("--markers", required=True)
    marker_frequency.add_argument("--manifest", required=True)
    marker_frequency.add_argument("--config", default="configs/default.yaml")
    marker_frequency.add_argument("--output")
    marker_frequency.add_argument("--run-id")
    marker_frequency.add_argument("--skip-normalization", action="store_true")
    marker_frequency.set_defaults(handler=command_markers_frequency)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.handler(args))
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
