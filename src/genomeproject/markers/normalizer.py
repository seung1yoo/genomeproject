from __future__ import annotations

import csv
import shutil
import subprocess
import tempfile
from pathlib import Path

import pysam

from .models import Marker


REQUIRED_COLUMNS = ("marker_id", "chrom", "pos", "ref", "alt")


def canonical_autosome(chrom: str) -> int | None:
    value = chrom[3:] if chrom.lower().startswith("chr") else chrom
    try:
        number = int(value)
    except ValueError:
        return None
    return number if 1 <= number <= 22 else None


def read_markers(path: str | Path, autosomes_only: bool = True) -> list[Marker]:
    marker_path = Path(path)
    with marker_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [column for column in REQUIRED_COLUMNS if column not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"Marker table is missing columns: {', '.join(missing)}")
        result: list[Marker] = []
        seen: set[tuple[str, int, str, str]] = set()
        seen_ids: set[str] = set()
        for line_no, raw in enumerate(reader, start=2):
            marker_id = (raw.get("marker_id") or "").strip()
            chrom = (raw.get("chrom") or "").strip()
            ref = (raw.get("ref") or "").strip().upper()
            alt = (raw.get("alt") or "").strip().upper()
            try:
                pos = int(raw.get("pos") or "")
            except ValueError as exc:
                raise ValueError(f"Invalid pos at line {line_no}") from exc
            if not marker_id or not chrom or pos < 1 or not ref or not alt:
                raise ValueError(f"Invalid marker at line {line_no}")
            if "," in alt:
                raise ValueError(f"One ALT per row is required at line {line_no}")
            if marker_id in seen_ids:
                raise ValueError(f"Duplicate marker_id at line {line_no}: {marker_id}")
            if autosomes_only and canonical_autosome(chrom) is None:
                raise ValueError(f"Only autosomes 1-22 are supported: {chrom} at line {line_no}")
            key = (chrom, pos, ref, alt)
            if key in seen:
                raise ValueError(f"Duplicate marker allele at line {line_no}: {chrom}:{pos}:{ref}:{alt}")
            seen.add(key)
            seen_ids.add(marker_id)
            result.append(Marker(marker_id, chrom, pos, ref, alt))
    if not result:
        raise ValueError("Marker table contains no data rows")
    return result


def validate_reference(markers: list[Marker], fasta_path: str | Path) -> None:
    fasta = Path(fasta_path)
    if not fasta.is_file():
        raise FileNotFoundError(f"Reference FASTA not found: {fasta}")
    with pysam.FastaFile(str(fasta)) as reference:
        references = set(reference.references)
        for marker in markers:
            contig = resolve_contig(marker.chrom, references)
            observed = reference.fetch(contig, marker.pos - 1, marker.pos - 1 + len(marker.ref)).upper()
            if observed != marker.ref:
                raise ValueError(
                    f"Reference mismatch for {marker.marker_id}: expected {marker.ref}, FASTA has {observed}"
                )


def resolve_contig(chrom: str, references: set[str]) -> str:
    candidates = [chrom]
    if chrom.lower().startswith("chr"):
        candidates.append(chrom[3:])
    else:
        candidates.append(f"chr{chrom}")
    matches = [candidate for candidate in candidates if candidate in references]
    if len(matches) != 1:
        raise ValueError(f"Cannot uniquely resolve contig '{chrom}' against input header/reference")
    return matches[0]


def normalize_markers(markers: list[Marker], fasta_path: str | Path) -> list[Marker]:
    executable = shutil.which("bcftools")
    if not executable:
        raise RuntimeError("bcftools is required for marker normalization but was not found in PATH")
    with pysam.FastaFile(str(fasta_path)) as reference:
        reference_lengths = dict(zip(reference.references, reference.lengths))
    resolved_markers = [
        (marker, resolve_contig(marker.chrom, set(reference_lengths))) for marker in markers
    ]
    with tempfile.TemporaryDirectory(prefix="genomeproject-markers-") as temp_dir:
        source = Path(temp_dir) / "markers.vcf"
        with source.open("w", encoding="utf-8") as handle:
            handle.write("##fileformat=VCFv4.2\n")
            for contig in dict.fromkeys(contig for _, contig in resolved_markers):
                handle.write(f"##contig=<ID={contig},length={reference_lengths[contig]}>\n")
            handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for marker, contig in resolved_markers:
                handle.write(
                    f"{contig}\t{marker.pos}\t{marker.marker_id}\t{marker.ref}\t{marker.alt}\t.\t.\t.\n"
                )
        process = subprocess.run(
            [executable, "norm", "-f", str(fasta_path), "-m", "-any", str(source)],
            check=False,
            text=True,
            capture_output=True,
        )
        if process.returncode:
            raise RuntimeError(f"bcftools norm failed: {process.stderr.strip()}")
        normalized: list[Marker] = []
        for line in process.stdout.splitlines():
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            normalized.append(Marker(fields[2], fields[0], int(fields[1]), fields[3], fields[4]))
        return normalized
