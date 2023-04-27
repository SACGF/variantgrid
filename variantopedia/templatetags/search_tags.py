from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import Count
from django.template import Library

from annotation.annotation_version_querysets import get_variant_queryset_for_latest_annotation_version
from annotation.models import patients_qs_for_ontology_term
from classification.models import Classification
from ontology.models import OntologyTerm
from snpdb.models import Variant, VariantZygosityCountCollection
from variantopedia.interesting_nearby import interesting_summary

register = Library()


@register.inclusion_tag("variantopedia/tags/search_summary.html")
def search_summary(user: User, search_result):
    record = search_result.record
    context = {}
    if isinstance(record, Variant):
        # FIXME do this as a preview
        summary, tag_counts = _variant_interesting_summary(user, record, search_result.genome_build)
        context["summary"] = summary
        context["tag_counts_dict"] = tag_counts
    return context


def _variant_interesting_summary(user: User, variant: Variant, genome_build, clinical_significance=False):
    qs = get_variant_queryset_for_latest_annotation_version(genome_build)
    qs, _ = VariantZygosityCountCollection.annotate_global_germline_counts(qs)
    qs = qs.filter(pk=variant.pk)

    return interesting_summary(qs, user, genome_build, total=False,
                               clinvar=settings.SEARCH_SUMMARY_VARIANT_SHOW_CLINVAR,
                               classifications=settings.SEARCH_SUMMARY_VARIANT_SHOW_CLASSIFICATIONS,
                               clinical_significance=clinical_significance)
