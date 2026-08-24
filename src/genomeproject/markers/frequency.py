from __future__ import annotations

from collections import Counter, defaultdict

from .models import GenotypeObservation, Marker


def sample_marker_stats(observations: list[GenotypeObservation]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for record in sorted(
        observations,
        key=lambda item: (item.sample_id, item.chrom, item.pos, item.ref, item.alt),
    ):
        rows.append(
            {
                "sample_id": record.sample_id,
                "marker_id": record.marker_id,
                "chrom": record.chrom,
                "pos": record.pos,
                "ref": record.ref,
                "alt": record.alt,
                "gt": record.gt,
                "ac": record.target_dosage,
                "an": record.allele_number,
                "af": record.target_dosage / record.allele_number
                if record.allele_number
                else None,
                "dp": record.dp,
                "gq": record.gq,
                "raw_called": record.raw_called,
                "qc_called": record.qc_called,
                "status": record.status,
                "qc_exclusion": record.qc_exclusion,
                "source_gvcf": record.source_gvcf,
            }
        )
    return rows


def aggregate_frequencies(
    markers: list[Marker], observations: list[GenotypeObservation]
) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, int, str, str], list[GenotypeObservation]] = defaultdict(list)
    for observation in observations:
        grouped[
            (observation.marker_id, observation.chrom, observation.pos, observation.ref, observation.alt)
        ].append(observation)

    results: list[dict[str, object]] = []
    for marker in markers:
        records = grouped.get(
            (marker.marker_id, marker.chrom, marker.pos, marker.ref, marker.alt), []
        )
        status_counts = Counter(record.status for record in records)
        exclusion_counts = Counter(record.qc_exclusion for record in records if record.qc_exclusion)
        row: dict[str, object] = {
            "marker_id": marker.marker_id,
            "chrom": marker.chrom,
            "pos": marker.pos,
            "ref": marker.ref,
            "alt": marker.alt,
            "cohort_n": len(records),
        }
        for prefix, predicate in (
            ("raw", lambda record: record.raw_called),
            ("qc", lambda record: record.qc_called),
        ):
            called = [record for record in records if predicate(record)]
            allele_number = sum(record.allele_number for record in called)
            allele_count = sum(record.target_dosage for record in called)
            row.update(
                {
                    f"{prefix}_called_n": len(called),
                    f"{prefix}_call_rate": len(called) / len(records) if records else None,
                    f"{prefix}_ac": allele_count,
                    f"{prefix}_an": allele_number,
                    f"{prefix}_af": allele_count / allele_number if allele_number else None,
                    f"{prefix}_carrier_n": sum(record.target_dosage > 0 for record in called),
                    f"{prefix}_hom_ref_n": sum(
                        record.allele_number == 2 and record.target_dosage == 0 for record in called
                    ),
                    f"{prefix}_het_n": sum(
                        record.allele_number == 2 and record.target_dosage == 1 for record in called
                    ),
                    f"{prefix}_hom_alt_n": sum(
                        record.allele_number == 2 and record.target_dosage == 2 for record in called
                    ),
                }
            )
        row["status_counts"] = ";".join(f"{key}={value}" for key, value in sorted(status_counts.items()))
        row["qc_exclusion_counts"] = ";".join(
            f"{key}={value}" for key, value in sorted(exclusion_counts.items())
        )
        results.append(row)
    return results
