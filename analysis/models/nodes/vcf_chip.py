""" The "VCF" chip a source node draws for the files its samples come from, with a coloured dot per kind of
    call those files hold (#2759). DRAGEN writes a small variant, CNV and SV VCF per patient, all naming the
    sample alike, so the dots are what tell two otherwise identical nodes apart.

    Entry points: get_vcf_contents (the counts, read from each sample's own VCF stats, so a custom cohort or
    family still shows what the original VCFs hold) and vcf_chip (the chip). Used by SampleNode and every
    cohort-based node (Cohort, Pedigree, Duo, Trio, Quad).
"""
from collections import Counter, defaultdict
from collections.abc import Iterable

from analysis.models.nodes.node_display import NodeChip
from snpdb.enums import VariantContent
from snpdb.models.models_cohort_stats import CohortGenotypeStats
from snpdb.models.models_vcf import Sample

VCF_CHIP_ICON = "fa-solid fa-file-lines"


def get_vcf_contents(samples: Iterable[Sample]) -> dict[int, Counter]:
    """ vcf_id -> kinds of call its samples hold. Read fresh rather than memoized with a card, as the stats
        that say it land after the node is made. A sample whose stats predate the counts leaves its VCF
        without dots until they're recomputed.

        Per-sample stats rows always sit on the sample's VCF cohort, so this is what the original VCF holds
        whichever cohort the node reads the sample through """
    stats_qs = CohortGenotypeStats.objects.filter(
        sample__in=samples, filter_key__isnull=True, passing_filter=False,
        code_version__version__gte=CohortGenotypeStats.VARIANT_CONTENTS_STATS_VERSION)
    # Latest per sample - a reload in place keeps the collection, so there is normally just one
    stats_qs = stats_qs.select_related("sample__vcf").order_by("sample_id", "-pk").distinct("sample_id")
    contents = defaultdict(Counter)
    for stats in stats_qs:
        vcf = stats.sample.vcf
        contents[vcf.pk].update(stats.get_variant_contents(has_copy_number=bool(vcf.copy_number_field)))
    return contents


def vcf_chip(vcfs: list[tuple[int, str]], vcf_contents: dict[int, Counter]) -> NodeChip:
    """ Always one chip for (vcf_id, name) - a card has no room for a pill per caller, so the names are the
        hover, each with what it holds, and the chip's dots are what the files hold between them """
    lines = []
    totals = Counter()
    for vcf_id, name in vcfs:
        contents = vcf_contents.get(vcf_id, {})
        totals.update(contents)
        if len(vcfs) > 1 and contents:
            name += " - " + ", ".join(content.label for content in VariantContent if content in contents)
        lines.append(name)
    dots = tuple(NodeChip(text=content.label, css_class=f"variant-content-{content.value}",
                          title=f"{totals[content]:,} {'call' if totals[content] == 1 else 'calls'}")
                 for content in VariantContent if content in totals)
    return NodeChip(text="VCF", icon=VCF_CHIP_ICON, title="\n".join(lines),
                    count=len(vcfs) if len(vcfs) > 1 else None, dots=dots)
