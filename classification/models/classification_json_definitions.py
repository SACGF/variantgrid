from typing import Any, Optional, TypedDict


class ClassificationJsonVersionDict(TypedDict):
    version: float
    publish_level: str
    is_published: bool
    can_write: bool


class ClassificationJsonConfigDict(TypedDict):
    evidence_key_overrides: dict[str, dict[str, Any]]
    namespaces: list[str]
    allele_origin_bucket: str


class ClassificationJsonAlleleGenomeBuild(TypedDict, total=True):
    variant_id: int
    c_hgvs: str


class ClassificationJsonAlleleRevolvedDict(TypedDict, total=False):
    allele_id: int
    allele_info_id: int
    allele_info_status: str
    status: str  # remove me?
    include: Optional[bool]
    variant_coordinate: str
    # fields from c.HGVS, applies for preferred only
    transcript: Optional[str]
    gene_symbol: Optional[str]
    c_nomen: Optional[str]
    full: str
    genome_build: Optional[str]
    desired: bool
    normalized: bool
    # if there is a warning icon
    icon: str
    tooltip: str


class ClassificationJsonAlleleDict(TypedDict, total=False):
    resolved: ClassificationJsonAlleleRevolvedDict
    genome_builds: dict[str, ClassificationJsonAlleleGenomeBuild]


