""" Resolves free-text variant entries (one per line) to variants and regions.

    Identifier formats come from the search framework's receivers (dbSNP, HGVS, variant
    coordinates, loci, ClinGen allele IDs, COSMIC), so anything the search box understands
    can be pasted in, as can VEP's "&"-joined dbSNP ids. Regions are the one shape search doesn't cover, so they're parsed here.
"""
import re
from dataclasses import dataclass, field
from typing import Optional

from django.conf import settings
from django.contrib.auth.models import User

from library.genomics import format_chrom
from snpdb.models import GenomeBuild, Variant
from snpdb.models.models_variant import (
    VARIANT_GNOMAD_PATTERN,
    VARIANT_PATTERN,
    VARIANT_SYMBOLIC_PATTERN,
    VARIANT_VCF_PATTERN,
    VariantCoordinate,
)
from snpdb.search import search_data
from snpdb.variant_pk_lookup import VariantPKLookup

REGION_PATTERN = re.compile(r"^(?P<chrom>[^:\s]+)\s*:\s*(?P<start>[0-9,]+)\s*-\s*(?P<end>[0-9,]+)$")
REGION_SEPARATOR = "-"
# VEP joins a variant's dbSNP ids with "&", which is how they appear in grids and exports
DBSNP_IDS_PATTERN = re.compile(r"^rs\d+(&rs\d+)+$", re.IGNORECASE)
# The search box's coordinate formats - chrom and position are groups 1 and 2 of each
COORDINATE_PATTERNS = [
    (VARIANT_PATTERN, VariantCoordinate.from_variant_match),
    (VARIANT_VCF_PATTERN, VariantCoordinate.from_variant_match),
    (VARIANT_GNOMAD_PATTERN, VariantCoordinate.from_variant_match),
    (VARIANT_SYMBOLIC_PATTERN, VariantCoordinate.from_symbolic_match),
]


@dataclass
class VariantTextResolution:
    variant_ids: list[int] = field(default_factory=list)
    regions: list[str] = field(default_factory=list)
    unresolved: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)


def format_region(chrom: str, start: int, end: int) -> str:
    return f"{chrom}:{start}{REGION_SEPARATOR}{end}"


def parse_region(region: str) -> tuple[str, int, int]:
    """ Splits a region stored by resolve_variant_text() back into its parts """
    chrom, _, positions = region.partition(":")
    start, _, end = positions.partition(REGION_SEPARATOR)
    return chrom, int(start), int(end)


def _resolve_region(line: str, genome_build: GenomeBuild, resolution: VariantTextResolution) -> bool:
    if m := REGION_PATTERN.match(line):
        chrom = m.group("chrom")
        start = int(m.group("start").replace(",", ""))
        end = int(m.group("end").replace(",", ""))
        if contig := genome_build.chrom_contig_mappings.get(chrom):
            if start > end:
                resolution.errors.append(f"'{line}': start > end")
            else:
                resolution.regions.append(format_region(contig.name, start, end))
        else:
            resolution.errors.append(f"'{line}': chromosome/contig not in {genome_build}")
        return True
    return False


def _search_variant_ids(user: User, genome_build: GenomeBuild, line: str) -> list[int]:
    """ Search results span genome builds and object types - keep this build's variants """
    variant_ids = set()
    for result in search_data(user, line).results:
        preview = result.preview
        if preview.is_operation or preview.is_error:
            continue  # eg an offer to create a variant we don't have
        if isinstance(preview.obj, Variant):
            variant_ids.add(preview.obj.pk)
    if not variant_ids:
        return []
    qs = Variant.objects.filter(pk__in=variant_ids).filter(Variant.get_contigs_q(genome_build))
    return list(qs.values_list("pk", flat=True))


def _search_line_variant_ids(user: User, genome_build: GenomeBuild, line: str) -> list[int]:
    """ A line of dbSNP ids (rs1&rs2) is searched one id at a time, until one finds variants """
    if DBSNP_IDS_PATTERN.match(line):
        for dbsnp_rs_id in line.split("&"):
            if variant_ids := _search_variant_ids(user, genome_build, dbsnp_rs_id):
                return variant_ids
        return []
    return _search_variant_ids(user, genome_build, line)


def _variant_coordinate(line: str, genome_build: GenomeBuild) -> Optional[VariantCoordinate]:
    for pattern, get_variant_coordinate in COORDINATE_PATTERNS:
        if m := pattern.fullmatch(line):
            chrom = format_chrom(m.group(1), genome_build.reference_fasta_has_chr)
            if Variant.validate(genome_build, chrom, m.group(2)):
                return None
            try:
                return get_variant_coordinate(m, genome_build)
            except ValueError:
                return None
    return None


def _exact_coordinate_variant_ids(user: User, genome_build: GenomeBuild, lines: list[str]) -> dict[str, int]:
    """ Line -> variant for coordinates already in the database, in batched queries - a search per line
        takes ~0.5s (#2531). Lines left out (a REF the search corrects, an rsID, HGVS etc) go to the search """
    if settings.SEARCH_VARIANT_REQUIRE_CLASSIFICATION_FOR_NON_ADMIN and not user.is_superuser:
        return {}  # Search only shows classified variants - leave that to it

    coordinate_by_line = {}
    for line in lines:
        if variant_coordinate := _variant_coordinate(line, genome_build):
            coordinate_by_line[line] = variant_coordinate

    variant_pk_lookup = VariantPKLookup(genome_build)
    variant_hash_by_line = {line: variant_pk_lookup.get_variant_coordinate_hash(variant_coordinate)
                            for line, variant_coordinate in coordinate_by_line.items()}
    variant_ids = variant_pk_lookup.get_variant_ids(variant_hash_by_line.values(), validate_not_null=False)
    variant_id_by_line = {line: variant_id for line, variant_id in zip(variant_hash_by_line, variant_ids) if variant_id}
    return variant_id_by_line


def resolve_variant_text(user: User, genome_build: GenomeBuild, text: str) -> VariantTextResolution:
    """ One entry per line. A line that resolves to nothing is reported rather than dropped """
    resolution = VariantTextResolution()
    lines = [stripped for line in (text or "").splitlines() if (stripped := line.strip())]
    variant_id_by_line = _exact_coordinate_variant_ids(user, genome_build, lines)
    seen_variant_ids = set()
    for line in lines:
        if _resolve_region(line, genome_build, resolution):
            continue
        if variant_id := variant_id_by_line.get(line):
            variant_ids = [variant_id]
        else:
            variant_ids = _search_line_variant_ids(user, genome_build, line)
        if variant_ids:
            for variant_id in variant_ids:
                if variant_id not in seen_variant_ids:
                    seen_variant_ids.add(variant_id)
                    resolution.variant_ids.append(variant_id)
        else:
            resolution.unresolved.append(line)
    return resolution


def get_variant_text_summary(num_entries: int, variant_ids: list[int], regions: list[str],
                             unresolved: list[str]) -> Optional[str]:
    if not num_entries:
        return None
    resolved = num_entries - len(unresolved)
    summary = f"{resolved}/{num_entries} entries resolved to {len(variant_ids)} variants"
    if regions:
        summary += f" and {len(regions)} regions"
    return summary
