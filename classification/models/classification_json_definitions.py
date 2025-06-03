from typing import Optional, TypedDict, Any

"""
This file was made do give definitions to the JSON returned by the classification API.
A few of the definitions are used, but not all.
"""


class ClassificationJsonLabDict(TypedDict):
    group_name: str
    lab_name: str
    org_name: str


class ClassificationJsonVersionDict(TypedDict):
    version: float
    publish_level: str
    is_published: bool
    can_write: bool


class ClassificationJsonConfigDict(TypedDict):
    evidence_key_overrides: dict[str, dict[str, Any]]
    namespaces: list[str]
    allele_origin_bucket: str


class ClassificationJsonSampleDict(TypedDict):
    id: str
    name: str


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


class ClassificationJsonDictv3(TypedDict, total=False):
    id: int
    lab_record_id: str
    cr_lab_id: str
    data: dict[str, Any]
    messages: Optional[list]
    version: ClassificationJsonVersionDict
    version_published: ClassificationJsonVersionDict
    version_latest: ClassificationJsonVersionDict

    latest_version: ClassificationJsonVersionDict
    config: ClassificationJsonConfigDict
    allele_info: ClassificationJsonAlleleDict
