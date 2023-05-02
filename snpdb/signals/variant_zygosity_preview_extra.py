from django.conf import settings
from django.dispatch import receiver
from django.utils.safestring import SafeString
from registration.forms import User

from annotation.annotation_version_querysets import get_variant_queryset_for_latest_annotation_version
from library.log_utils import get_current_logged_in_user
from library.preview_request import preview_extra_signal, PreviewKeyValue
from library.utils import first
from snpdb.models import Allele, Variant, VariantZygosityCountCollection
from variantopedia.interesting_nearby import interesting_summary


def _variant_preview_zygosity_extar(variant: Variant):
    if genome_build := first(variant.genome_builds):
        qs = get_variant_queryset_for_latest_annotation_version(first(variant.genome_builds))
        qs, _ = VariantZygosityCountCollection.annotate_global_germline_counts(qs)
        qs = qs.filter(pk=variant.pk)

        return interesting_summary(qs, get_current_logged_in_user(), genome_build, total=False,
                                   clinvar=settings.SEARCH_SUMMARY_VARIANT_SHOW_CLINVAR,
                                   classifications=settings.SEARCH_SUMMARY_VARIANT_SHOW_CLASSIFICATIONS,
                                   clinical_significance=False)


@receiver(preview_extra_signal, sender=Allele)
def allele_preview_classifications_extra(sender, user: User, obj: Allele, **kwargs):
    # return [PreviewKeyValue(value=SafeString("This is <b>BOLD</b>"))]
    pass