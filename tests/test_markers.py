from pathlib import Path

import pysam
import pytest

from genomeproject.markers.frequency import aggregate_frequencies
from genomeproject.markers.gvcf_reader import read_sample_gvcf
from genomeproject.markers.models import Marker
from genomeproject.markers.normalizer import normalize_markers, read_markers, validate_reference


def make_reference(tmp_path: Path) -> Path:
    fasta = tmp_path / "GRCh38.fa"
    fasta.write_text(">chr1\nC" + "A" * 99 + "\n", encoding="utf-8")
    pysam.faidx(str(fasta))
    return fasta


def make_gvcf(tmp_path: Path, sample: str, low_gq: bool = False) -> Path:
    plain = tmp_path / f"{sample}.g.vcf"
    gq = 10 if low_gq else 40
    plain.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=100>\n"
        "##INFO=<ID=END,Number=1,Type=Integer,Description=End>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=GenotypeQuality>\n"
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Depth>\n"
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
        "chr1\t1\t.\tC\t<NON_REF>\t.\t.\tEND=9\tGT:GQ:DP\t0/0:40:20\n"
        f"chr1\t10\t.\tA\tG,<NON_REF>\t.\t.\t.\tGT:GQ:DP\t0/1:{gq}:30\n"
        "chr1\t11\t.\tA\tT,<NON_REF>\t.\t.\t.\tGT:GQ:DP\t0/1:50:30\n"
        "chr1\t12\t.\tA\t<NON_REF>\t.\t.\tEND=20\tGT:GQ:DP\t0/0:10:30\n",
        encoding="utf-8",
    )
    compressed = tmp_path / f"{sample}.g.vcf.gz"
    pysam.tabix_compress(str(plain), str(compressed), force=True)
    pysam.tabix_index(str(compressed), preset="vcf", force=True)
    return compressed


def test_direct_gvcf_read_and_frequency(tmp_path: Path) -> None:
    gvcf1 = make_gvcf(tmp_path, "S1")
    gvcf2 = make_gvcf(tmp_path, "S2", low_gq=True)
    markers = [
        Marker("ref", "chr1", 5, "A", "G"),
        Marker("target", "chr1", 10, "A", "G"),
        Marker("other", "chr1", 11, "A", "G"),
        Marker("low", "chr1", 12, "A", "G"),
        Marker("absent", "chr1", 25, "A", "G"),
    ]
    observations = read_sample_gvcf("S1", gvcf1, markers, 20, 8)
    observations += read_sample_gvcf("S2", gvcf2, markers, 20, 8)
    by_marker_sample = {(item.marker_id, item.sample_id): item for item in observations}

    assert by_marker_sample[("ref", "S1")].status == "REFERENCE_BLOCK"
    assert by_marker_sample[("target", "S1")].target_dosage == 1
    assert by_marker_sample[("other", "S1")].status == "OTHER_ALT"
    assert by_marker_sample[("other", "S1")].target_dosage == 0
    assert by_marker_sample[("low", "S1")].qc_exclusion == "LOW_GQ"
    assert by_marker_sample[("absent", "S1")].status == "NO_RECORD"
    assert not by_marker_sample[("absent", "S1")].raw_called

    frequencies = {row["marker_id"]: row for row in aggregate_frequencies(markers, observations)}
    assert frequencies["target"]["raw_ac"] == 2
    assert frequencies["target"]["raw_an"] == 4
    assert frequencies["target"]["raw_af"] == 0.5
    assert frequencies["target"]["qc_ac"] == 1
    assert frequencies["target"]["qc_an"] == 2
    assert frequencies["absent"]["raw_an"] == 0


def test_marker_table_and_reference_validation(tmp_path: Path) -> None:
    fasta = make_reference(tmp_path)
    marker_table = tmp_path / "markers.tsv"
    marker_table.write_text(
        "marker_id\tchrom\tpos\tref\talt\nM1\t1\t5\tA\tG\n", encoding="utf-8"
    )
    markers = read_markers(marker_table)
    validate_reference(markers, fasta)
    assert markers[0].pos == 5

    marker_table.write_text(
        "marker_id\tchrom\tpos\tref\talt\nM1\tchrX\t5\tA\tG\n", encoding="utf-8"
    )
    with pytest.raises(ValueError, match="Only autosomes"):
        read_markers(marker_table)


def test_normalize_markers_defines_reference_contig(tmp_path: Path) -> None:
    fasta = make_reference(tmp_path)
    markers = [Marker("M1", "1", 5, "A", "G")]

    normalized = normalize_markers(markers, fasta)

    assert normalized == [Marker("M1", "chr1", 5, "A", "G")]
