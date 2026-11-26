""" What kinds of call a sample or a VCF holds (snpdb.enums.VariantContent: Small, CNV, SV, Chrom CN, Fusion,
    Gene CN, Splice), read from the per-sample stats (#2759). DRAGEN writes a small variant, CNV and SV VCF per patient,
    all naming the sample alike, so this is what tells them apart wherever samples and VCFs are listed.

    A sample whose stats predate the counts (CohortGenotypeStats.VARIANT_CONTENTS_STATS_VERSION) has no
    contents until its stats are recomputed, so every caller copes with a sample or VCF missing from the result.

    Entry points: get_sample_variant_contents and get_vcf_variant_contents (the counts), variant_content_exists
    (a grid filter), and variant_contents_chips / variant_contents_text / VARIANT_CONTENTS_HELP (how they are
    shown - chips are drawn by the variant_contents template tag and TableFormat.variantContents, text is for
    an autocomplete result).
"""
from collections.abc import Iterable
from typing import Optional, Union

from django.db.models import Exists, OuterRef, QuerySet, Sum

from library.utils.json_utils import JsonObjType
from snpdb.enums import VariantContent
from snpdb.models.models_cohort_stats import CohortGenotypeStats
from snpdb.models.models_vcf import VCF, Sample

VariantContents = dict[VariantContent, int]

# One line per kind - joined with newlines for a title attribute, <br> for a help popover
VARIANT_CONTENTS_HELP_LINES = [
    "What kinds of call the VCF holds, counted when its sample stats are calculated.",
    *(f"{content.label}: {content.description}" for content in VariantContent),
]
VARIANT_CONTENTS_HELP = "\n".join(VARIANT_CONTENTS_HELP_LINES)

_COUNT_FIELDS = ("variant_count", "structural_count", "whole_chromosome_count", "fusions_count",
                 "gene_copy_number_count", "gene_splice_count")


def latest_sample_variant_stats(samples: Union[Iterable[Sample], QuerySet[Sample]]) -> QuerySet[CohortGenotypeStats]:
    """ One stats row per sample that has the counts. Per-sample rows always sit on the sample's VCF cohort,
        so this is what the original VCF holds whichever cohort the sample is being read through """
    stats_qs = CohortGenotypeStats.variant_contents_qs().filter(sample__in=samples)
    # Latest per sample - a reload in place keeps the collection, so there is normally just one
    return stats_qs.select_related("sample__vcf").order_by("sample_id", "-pk").distinct("sample_id")


def get_sample_variant_contents(samples: Union[Iterable[Sample], QuerySet[Sample]]) -> dict[int, VariantContents]:
    """ sample_id -> count per kind of call the sample holds """
    return {stats.sample_id: stats.get_variant_contents(has_copy_number=bool(stats.sample.vcf.copy_number_field))
            for stats in latest_sample_variant_stats(samples)}


def get_vcf_variant_contents(vcfs: Union[Iterable[VCF], Iterable[int], QuerySet[VCF]]) -> dict[int, VariantContents]:
    """ vcf_id -> count per kind of call, summed over the VCF's samples in SQL as a VCF can hold thousands.
        The sum is calls rather than records, so it says which kinds are there and not how many variants """
    sums = {f"total_{f}": Sum(f) for f in _COUNT_FIELDS}
    totals_qs = (CohortGenotypeStats.variant_contents_qs().filter(sample__vcf__in=vcfs)
                 .values("sample__vcf_id", "sample__vcf__copy_number_field").annotate(**sums))
    vcf_contents = {}
    for totals in totals_qs:
        stats = CohortGenotypeStats(**{f: totals[f"total_{f}"] for f in _COUNT_FIELDS})
        has_copy_number = bool(totals["sample__vcf__copy_number_field"])
        vcf_contents[totals["sample__vcf_id"]] = stats.get_variant_contents(has_copy_number=has_copy_number)
    return vcf_contents


def variant_content_exists(content: VariantContent, **stats_filter) -> Exists:
    """ For a Sample (sample=OuterRef("pk")) or VCF (sample__vcf=OuterRef("pk")) queryset: holds this kind of call """
    stats_qs = CohortGenotypeStats.variant_contents_qs().filter(**stats_filter)
    return Exists(stats_qs.filter(CohortGenotypeStats.get_variant_content_q(content)))


def sample_variant_content_exists(content: VariantContent) -> Exists:
    return variant_content_exists(content, sample=OuterRef("pk"))


def vcf_variant_content_exists(content: VariantContent) -> Exists:
    return variant_content_exists(content, sample__vcf=OuterRef("pk"))


def variant_contents_chips(contents: Optional[VariantContents], counts: bool = True) -> list[JsonObjType]:
    """ A chip per kind held, in VariantContent order. counts=False where the numbers are sums over samples """
    chips = []
    for content in VariantContent:
        if count := (contents or {}).get(content):
            title = content.description
            if counts:
                title += f" - {count:,} {'call' if count == 1 else 'calls'}"
            chips.append({"content": content.value, "label": content.label, "title": title,
                          "count": count if counts else None})
    return chips


def variant_contents_text(contents: Optional[VariantContents]) -> str:
    """ "Small, SV" - for where a chip can't be drawn (an autocomplete result, a node chip's hover) """
    return ", ".join(content.label for content in VariantContent if content in (contents or {}))


def label_with_variant_contents(label: str, contents: Optional[VariantContents]) -> str:
    """ "SMP5008 (dragen.sv.vcf) — Small, SV" - the label alone where nothing has been counted """
    if text := variant_contents_text(contents):
        return f"{label} \u2014 {text}"
    return label
