from django.conf import settings
from django.dispatch import receiver
from django.template.loader import render_to_string
from django.utils.safestring import SafeString
from registration.forms import User

from annotation.annotation_version_querysets import get_variant_queryset_for_latest_annotation_version
from library.log_utils import get_current_logged_in_user
from library.preview_request import preview_extra_signal, PreviewKeyValue
from library.utils import first
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.models import Allele, Variant, VariantZygosityCountCollection, GenomeBuild
from variantopedia.interesting_nearby import interesting_summary


def _variant_preview_zygosity_extra(variant: Variant, genome_build: GenomeBuild):
    qs = get_variant_queryset_for_latest_annotation_version(first(variant.genome_builds))
    qs, _ = VariantZygosityCountCollection.annotate_global_germline_counts(qs)
    qs = qs.filter(pk=variant.pk)

    # FIXME we really want clinvar to be a PreviewKeyValue instead of a plain string
    summary, tag_counts = interesting_summary(qs, get_current_logged_in_user(), genome_build, total=False,
                               clinvar=settings.SEARCH_SUMMARY_VARIANT_SHOW_CLINVAR,
                               classifications=False,  # handled elsewhere
                               clinical_significance=False)

    context = {
        "summary": summary,
        "tag_counts": tag_counts
    }
    html = render_to_string('variantopedia/tags/search_summary.html', context).strip()
    if html:
        return [PreviewKeyValue(value=SafeString(html))]


@receiver(preview_extra_signal, sender=Allele)
def allele_preview_classifications_extra(sender, user: User, obj: Allele, **kwargs):
    genome_build = GenomeBuildManager.get_current_genome_build()
    if variant := obj.variant_for_build_optional(genome_build):
        return _variant_preview_zygosity_extra(variant, genome_build)


@receiver(preview_extra_signal, sender=Variant)
def allele_preview_classifications_extra(sender, user: User, obj: Variant, **kwargs):
    # choice of GenomeBuild with a single variant is
    genome_build = first(obj.genome_builds)
    if variant := obj.variant_for_build_optional(genome_build):
        return _variant_preview_zygosity_extra(variant, genome_build)