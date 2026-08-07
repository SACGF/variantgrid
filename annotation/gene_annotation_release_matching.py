"""
Work out which GeneAnnotationRelease matches the VEP a VariantAnnotationVersion was built
against.

VEP reports the gene set it used as version strings (eg RefSeq
'GCF_000001405.40-RS_2025_08 - GCF_000001405.40_GRCh38.p14_genomic.gff'). From those we
derive:

    release token - the identifier cdot uses in its per-GFF release asset names, eg
                    'RS_2025_08'. @see genes.cdot_data_release
    GFF/GTF url   - the file cdot read, which it records on every transcript and which ends
                    up on GeneAnnotationImport.url. This identifies a GeneAnnotationRelease
                    across deployments, unlike GeneAnnotationRelease.version which is a
                    hand-typed label.

Kept separate from the model so the derivation can be tested without a database - creating a
VariantAnnotationVersion builds partitions, which makes table-driven tests over every VEP
version string slow.
"""
import re
from dataclasses import dataclass
from typing import Optional

from genes.gene_annotation_import_urls import canonical_import_url
from genes.models_enums import AnnotationConsortium

# The one RefSeq value that reports a timestamp instead of an annotation release
LEGACY_GRCH37_REFSEQ = {
    "2020-10-26 17:03:42": "105.20201022",
}

# Ensembl stopped updating GRCh37 at release 87
GRCH37_ENSEMBL_RELEASE = "87"

# 'GCF_000001405.40-RS_2025_08' -> 'RS_2025_08'
_REFSEQ_ASSEMBLY_PREFIX_PATTERN = re.compile(r"^GCF_\d+\.\d+-")

_NCBI_ANNOTATION_RELEASES = "https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606"
_ENSEMBL_PUB = "https://ftp.ensembl.org/pub"


@dataclass(frozen=True)
class VepGeneSetVersions:
    """ The VariantAnnotationVersion fields describing which gene set VEP used """
    genome_build_name: str
    annotation_consortium: str
    refseq: str
    vep: Optional[int]
    vep_cache: Optional[int]
    genebuild: str

    @classmethod
    def from_variant_annotation_version(cls, vav) -> "VepGeneSetVersions":
        return cls(genome_build_name=vav.genome_build_id,
                   annotation_consortium=vav.annotation_consortium,
                   refseq=vav.refseq or "",
                   vep=vav.vep,
                   vep_cache=vav.vep_cache,
                   genebuild=vav.genebuild or "")

    @property
    def _refseq_release_and_filename(self) -> tuple[Optional[str], Optional[str]]:
        """ eg ('GCF_000001405.40-RS_2025_08', 'GCF_000001405.40_GRCh38.p14_genomic.gff') """
        release, separator, gff_filename = self.refseq.partition(" - ")
        if not separator:
            return None, None
        release = LEGACY_GRCH37_REFSEQ.get(release, release)
        return release, gff_filename

    @property
    def release_token(self) -> Optional[str]:
        """ The identifier cdot uses in per-GFF release asset names. None if unknown """
        if self.annotation_consortium == AnnotationConsortium.REFSEQ:
            release, _gff_filename = self._refseq_release_and_filename
            if release is None:
                return None
            return _REFSEQ_ASSEMBLY_PREFIX_PATTERN.sub("", release)

        if self.genome_build_name == "GRCh37":
            return GRCH37_ENSEMBL_RELEASE
        if self.genome_build_name == "GRCh38":
            if cache_version := self.vep_cache or self.vep:
                return str(cache_version)
            return None
        # Ensembl rapid-release (T2T) is versioned by geneset date, eg '2022-06' -> '2022_06'
        if self.genebuild:
            return self.genebuild.replace("-", "_")
        return None

    @property
    def gff_url(self) -> Optional[str]:
        """ The GFF/GTF cdot read for this gene set. None where we can't derive it """
        if self.annotation_consortium == AnnotationConsortium.REFSEQ:
            return self._refseq_gff_url()
        return self._ensembl_gtf_url()

    def _refseq_gff_url(self) -> Optional[str]:
        release, gff_filename = self._refseq_release_and_filename
        if release is None:
            return None
        if not gff_filename.endswith(".gz"):
            gff_filename += ".gz"

        # Releases named after the assembly hold the GFF directly; older numbered releases
        # (eg '110', '105.20220307') nest it under an assembly directory
        if release.startswith("GCF_"):
            path = f"{release}/{gff_filename}"
        else:
            assembly_dir = gff_filename.removesuffix("_genomic.gff.gz")
            path = f"{release}/{assembly_dir}/{gff_filename}"
        return f"{_NCBI_ANNOTATION_RELEASES}/{path}"

    def _ensembl_gtf_url(self) -> Optional[str]:
        release = self.release_token
        if release is None:
            return None
        if self.genome_build_name == "GRCh37":
            return (f"{_ENSEMBL_PUB}/grch37/release-{release}/gtf/homo_sapiens/"
                    f"Homo_sapiens.GRCh37.{release}.gtf.gz")
        if self.genome_build_name == "GRCh38":
            return (f"{_ENSEMBL_PUB}/release-{release}/gtf/homo_sapiens/"
                    f"Homo_sapiens.GRCh38.{release}.gtf.gz")
        # Rapid-release URLs embed a GCA accession we don't hold, so leave these to token
        # matching
        return None


def _url_contains_release(url: str, release_token: str) -> bool:
    """ Whether a path segment of url names release_token, eg '.../9606/110/...',
        '.../GCF_000001405.40-RS_2025_08/...', '.../release-116/...' """
    for segment in url.split("/"):
        if segment == release_token or segment.endswith("-" + release_token):
            return True
    return False


@dataclass(frozen=True)
class GeneAnnotationReleaseMatch:
    release_ids: list[int]
    description: str

    @property
    def unique_release_id(self) -> Optional[int]:
        if len(self.release_ids) == 1:
            return self.release_ids[0]
        return None


def match_gene_annotation_release(release_urls_by_id: dict[int, str],
                                  versions: VepGeneSetVersions) -> GeneAnnotationReleaseMatch:
    """ release_urls_by_id: GeneAnnotationRelease pk -> gene_annotation_import.url, already
        restricted to the genome build and annotation consortium being matched """
    if gff_url := versions.gff_url:
        canonical_url = canonical_import_url(gff_url)
        release_ids = [pk for pk, url in release_urls_by_id.items()
                       if canonical_import_url(url) == canonical_url]
        if release_ids:
            return GeneAnnotationReleaseMatch(release_ids, f"GFF url '{canonical_url}'")

    if release_token := versions.release_token:
        release_ids = [pk for pk, url in release_urls_by_id.items()
                       if _url_contains_release(url, release_token)]
        if release_ids:
            return GeneAnnotationReleaseMatch(release_ids, f"release '{release_token}' in url")

    return GeneAnnotationReleaseMatch([], "no match")
