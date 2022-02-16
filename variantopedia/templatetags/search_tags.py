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


@register.simple_tag(takes_context=True)
def search_summary(context, search_result):
    user = context["user"]
    summary = ""
    record = search_result.record
    if isinstance(record, Variant):
        summary = _variant_interesting_summary(user, record, search_result.genome_build)
    elif isinstance(record, OntologyTerm):
        summary = _get_ontology_summary(user, record)

    return summary


def _get_ontology_summary(user: User, ontology_term) -> str:
    ontology_summary_list = []
    terms = [{"term_id": ontology_term.pk}]
    qs = Classification.filter_for_user(user).filter(condition_resolution__resolved_terms__contains=terms)
    if num_classifications := qs.count():
        classification_summary = f"Condition for {num_classifications} classification"
        if num_classifications > 1:
            classification_summary += "s"
        ontology_summary_list.append(classification_summary)

    patients_qs = patients_qs_for_ontology_term(user, ontology_term)
    data = patients_qs.aggregate(num_patients=Count("id", distinct=True), num_samples=Count("sample", distinct=True))
    if num_patients := data.get("num_patients"):
        patient_summary_list = [f"{num_patients} patient"]
        if num_patients > 1:
            patient_summary_list.append("s")
        if num_samples := data.get("num_samples"):
            patient_summary_list.append(f" with {num_samples} samples")
        ontology_summary_list.append("".join(patient_summary_list))

    return ", ".join(ontology_summary_list)


def _variant_interesting_summary(user: User, variant: Variant, genome_build, clinical_significance=False) -> str:
    qs = get_variant_queryset_for_latest_annotation_version(genome_build)
    qs, _ = VariantZygosityCountCollection.annotate_global_germline_counts(qs)
    qs = qs.filter(pk=variant.pk)

    return interesting_summary(qs, user, genome_build, total=False,
                               clinvar=settings.SEARCH_SUMMARY_VARIANT_SHOW_CLINVAR,
                               classifications=settings.SEARCH_SUMMARY_VARIANT_SHOW_CLASSIFICATIONS,
                               clinical_significance=clinical_significance)
