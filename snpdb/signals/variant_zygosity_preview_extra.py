from collections import defaultdict
from typing import Optional

from django.conf import settings
from django.dispatch import receiver
from django.template.loader import render_to_string
from django.utils.safestring import SafeString
from registration.forms import User

from annotation.annotation_version_querysets import (
    get_variant_queryset_for_latest_annotation_version,
)
from library.preview_request import PreviewKeyValue, preview_extra_signal
from library.utils import first
from snpdb.models import Allele, GenomeBuild, Variant, VariantZygosityCountCollection
from variantopedia.interesting_nearby import interesting_summary


def _variant_zygosity_summary_html(user: User, variant: Variant, genome_build: GenomeBuild) -> Optional[str]:
    qs = get_variant_queryset_for_latest_annotation_version(genome_build)
    qs, _ = VariantZygosityCountCollection.annotate_global_germline_counts(qs)
    qs = qs.filter(pk=variant.pk)

    # FIXME we really want clinvar to be a PreviewKeyValue instead of a plain string
    summary, tag_counts = interesting_summary(qs, user, genome_build, total=False,
                                              clinvar=settings.SEARCH_SUMMARY_VARIANT_SHOW_CLINVAR,
                                              classifications=False,  # handled elsewhere
                                              clinical_significance=False)

    context = {
        "summary": summary,
        "tag_counts_dict": tag_counts
    }
    return render_to_string('variantopedia/tags/search_summary.html', context).strip()


@receiver(preview_extra_signal, sender=Allele)
def allele_preview_zygosity_extra(sender, user: User, obj: Allele, **kwargs):
    """ Counts are stored per variant, and the same mutation is a separate variant in each build,
        so show a row per build. Builds sharing a contig (MT, unplaced scaffolds) share the one
        variant - and so the one set of counts - so they're labelled together """

    builds_with_annotation = set(GenomeBuild.builds_with_annotation())
    builds_by_variant = defaultdict(list)
    for va in obj.variant_alleles():  # Ordered by build name
        if va.genome_build in builds_with_annotation:
            builds_by_variant[va.variant].append(va.genome_build)

    # A single build speaks for itself, but a shared contig covers two builds from the one variant
    show_genome_build = sum(len(gbs) for gbs in builds_by_variant.values()) > 1
    extra = []
    for variant, genome_builds in builds_by_variant.items():
        if html := _variant_zygosity_summary_html(user, variant, genome_builds[0]):
            key = ", ".join(str(gb) for gb in genome_builds) if show_genome_build else None
            extra.append(PreviewKeyValue(key=key, value=SafeString(html), dedicated_row=show_genome_build))
    return extra


@receiver(preview_extra_signal, sender=Variant)
def variant_preview_zygosity_extra(sender, user: User, obj: Variant, **kwargs):
    # choice of GenomeBuild with a single variant is
    genome_build = first(obj.genome_builds)
    if html := _variant_zygosity_summary_html(user, obj, genome_build):
        return [PreviewKeyValue(value=SafeString(html))]
    return None
