""" Resolves free-text variant entries (one per line) to variants and regions.

    Identifier formats come from the search framework's receivers (dbSNP, HGVS, variant
    coordinates, loci, ClinGen allele IDs, COSMIC), so anything the search box understands
    can be pasted in. Regions are the one shape search doesn't cover, so they're parsed here.
"""
import re
from dataclasses import dataclass, field
from typing import Optional

from django.contrib.auth.models import User

from snpdb.models import GenomeBuild, Variant
from snpdb.search import search_data

REGION_PATTERN = re.compile(r"^(?P<chrom>[^:\s]+)\s*:\s*(?P<start>[0-9,]+)\s*-\s*(?P<end>[0-9,]+)$")
REGION_SEPARATOR = "-"


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


def resolve_variant_text(user: User, genome_build: GenomeBuild, text: str) -> VariantTextResolution:
    """ One entry per line. A line that resolves to nothing is reported rather than dropped """
    resolution = VariantTextResolution()
    seen_variant_ids = set()
    for line in (text or "").splitlines():
        if not (line := line.strip()):
            continue
        if _resolve_region(line, genome_build, resolution):
            continue
        if variant_ids := _search_variant_ids(user, genome_build, line):
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
