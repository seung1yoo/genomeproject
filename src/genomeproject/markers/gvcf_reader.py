from __future__ import annotations

from pathlib import Path
from typing import Any

import pysam

from .models import GenotypeObservation, Marker
from .normalizer import resolve_contig


REFERENCE_ALTS = {"<NON_REF>", "<*>"}


def _integer(value: Any) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _gt_text(gt: tuple[int | None, ...] | None, phased: bool) -> str | None:
    if gt is None:
        return None
    separator = "|" if phased else "/"
    return separator.join("." if allele is None else str(allele) for allele in gt)


def _observation(
    sample_id: str,
    marker: Marker,
    source: str,
    status: str,
    min_gq: int,
    min_dp: int,
    sample_call: Any | None = None,
    target_index: int | None = None,
) -> GenotypeObservation:
    gt = tuple(sample_call.get("GT") or ()) if sample_call is not None else ()
    gq = _integer(sample_call.get("GQ")) if sample_call is not None else None
    dp = _integer(sample_call.get("DP")) if sample_call is not None else None
    phased = bool(getattr(sample_call, "phased", False)) if sample_call is not None else False
    called_alleles = [allele for allele in gt if allele is not None]
    raw_called = bool(gt) and len(called_alleles) == len(gt)
    dosage = sum(allele == target_index for allele in called_alleles) if target_index is not None else 0
    exclusion: str | None = None
    if not raw_called:
        exclusion = "NO_CALL"
    elif gq is None:
        exclusion = "MISSING_GQ"
    elif gq < min_gq:
        exclusion = "LOW_GQ"
    elif dp is None:
        exclusion = "MISSING_DP"
    elif dp < min_dp:
        exclusion = "LOW_DP"
    qc_called = exclusion is None
    return GenotypeObservation(
        sample_id=sample_id,
        marker_id=marker.marker_id,
        chrom=marker.chrom,
        pos=marker.pos,
        ref=marker.ref,
        alt=marker.alt,
        status=status,
        gt=_gt_text(gt, phased) if gt else None,
        gq=gq,
        dp=dp,
        allele_number=len(called_alleles) if raw_called else 0,
        target_dosage=dosage if raw_called else 0,
        raw_called=raw_called,
        qc_called=qc_called,
        qc_exclusion=exclusion,
        source_gvcf=source,
    )


def read_sample_gvcf(
    sample_id: str,
    gvcf_path: str | Path,
    markers: list[Marker],
    min_gq: int,
    min_dp: int,
) -> list[GenotypeObservation]:
    source = str(gvcf_path)
    output: list[GenotypeObservation] = []
    try:
        variant_file = pysam.VariantFile(source)
    except Exception:
        return [
            _observation(sample_id, marker, source, "OPEN_ERROR", min_gq, min_dp)
            for marker in markers
        ]
    with variant_file:
        samples = list(variant_file.header.samples)
        if len(samples) != 1:
            return [
                _observation(sample_id, marker, source, "SAMPLE_COUNT_ERROR", min_gq, min_dp)
                for marker in markers
            ]
        header_sample = samples[0]
        references = set(variant_file.header.contigs)
        for marker in markers:
            try:
                contig = resolve_contig(marker.chrom, references)
                records = [
                    record
                    for record in variant_file.fetch(contig, marker.pos - 1, marker.pos)
                    if record.start <= marker.pos - 1 < record.stop
                ]
            except (ValueError, OSError):
                output.append(_observation(sample_id, marker, source, "QUERY_ERROR", min_gq, min_dp))
                continue

            exact = [
                record
                for record in records
                if record.pos == marker.pos and record.ref.upper() == marker.ref and marker.alt in (record.alts or ())
            ]
            if len(exact) > 1:
                output.append(_observation(sample_id, marker, source, "AMBIGUOUS_RECORDS", min_gq, min_dp))
                continue
            if exact:
                record = exact[0]
                target_index = list(record.alleles).index(marker.alt)
                output.append(
                    _observation(
                        sample_id,
                        marker,
                        source,
                        "TARGET_VARIANT",
                        min_gq,
                        min_dp,
                        record.samples[header_sample],
                        target_index,
                    )
                )
                continue

            same_position = [record for record in records if record.pos == marker.pos]
            if any(record.ref.upper() != marker.ref for record in same_position):
                output.append(_observation(sample_id, marker, source, "REF_MISMATCH", min_gq, min_dp))
                continue
            callable_records = []
            for record in records:
                alts = set(record.alts or ())
                if record.pos < marker.pos:
                    if alts and alts <= REFERENCE_ALTS:
                        callable_records.append(record)
                elif marker.ref == record.ref.upper() and any(
                    alt in REFERENCE_ALTS for alt in alts
                ):
                    callable_records.append(record)
            if len(callable_records) != 1:
                status = "NO_RECORD" if not records else "UNSUPPORTED_OR_AMBIGUOUS"
                output.append(_observation(sample_id, marker, source, status, min_gq, min_dp))
                continue
            record = callable_records[0]
            call = record.samples[header_sample]
            gt = tuple(call.get("GT") or ())
            non_ref_indices = {
                index for index, allele in enumerate(record.alleles) if allele in REFERENCE_ALTS
            }
            if any(allele in non_ref_indices for allele in gt if allele is not None):
                output.append(_observation(sample_id, marker, source, "NON_REF_GENOTYPE", min_gq, min_dp))
                continue
            status = "REFERENCE_BLOCK" if record.pos < marker.pos or set(record.alts or ()) <= REFERENCE_ALTS else "OTHER_ALT"
            output.append(
                _observation(sample_id, marker, source, status, min_gq, min_dp, call, target_index=-1)
            )
    return output
