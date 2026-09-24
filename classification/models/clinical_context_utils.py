from collections.abc import Iterable

from django.db.models import F, Q, QuerySet

from classification.models import Classification
from classification.models.classification_grouping import ClassificationGrouping


def classifications_needing_rehoming(classification_qs: QuerySet[Classification] = None) -> QuerySet[Classification]:
    """ Classifications whose derived records sit under a different Allele than the classification itself -
        both are moved lazily on publish, so a bulk allele change (a merge, or the #1361 dedupe) strands them """
    if classification_qs is None:
        classification_qs = Classification.objects.all()

    grouping_allele = "classificationgroupingentry__grouping__allele_origin_grouping__allele"
    wrong_clinical_context = Q(clinical_context__isnull=False) & ~Q(clinical_context__allele=F("allele"))
    wrong_grouping = Q(classificationgroupingentry__isnull=False) & ~Q(**{grouping_allele: F("allele")})
    return classification_qs.filter(allele__isnull=False) \
        .filter(wrong_clinical_context | wrong_grouping).distinct()


def rehome_classifications(classifications: Iterable[Classification], force_recalc_text: str):
    """ Put clinical contexts and groupings back under the Allele their classification now points at """
    classifications = list(classifications)
    if not classifications:
        return 0

    for classification in classifications:
        ClassificationGrouping.assign_grouping_for_classification(classification)
    ClassificationGrouping.update_all_dirty()
    return len(classifications)
