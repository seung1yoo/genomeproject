from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class Marker:
    marker_id: str
    chrom: str
    pos: int
    ref: str
    alt: str

    @property
    def key(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"


@dataclass(frozen=True)
class GenotypeObservation:
    sample_id: str
    marker_id: str
    chrom: str
    pos: int
    ref: str
    alt: str
    status: str
    gt: str | None
    gq: int | None
    dp: int | None
    allele_number: int
    target_dosage: int
    raw_called: bool
    qc_called: bool
    qc_exclusion: str | None
    source_gvcf: str

    def as_dict(self) -> dict[str, object]:
        return asdict(self)
