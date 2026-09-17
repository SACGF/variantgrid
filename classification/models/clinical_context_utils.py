from collections.abc import Iterable
from typing import Optional

from django.db.models import F, Q, QuerySet

from classification.models import Classification, ClinicalContext, ClinicalContextRecalcTrigger
from classification.models.classification_grouping import ClassificationGrouping


def _assign_new_cc_reason(classification: Classification) -> Optional[ClinicalContextRecalcTrigger]:
    if existing_cc := classification.clinical_context:
        if existing_cc.allele != classification.allele_object:
            # we have a cc, but the allele doesn't match
            return ClinicalContextRecalcTrigger.VARIANT_SET
        if existing_cc.allele_origin_bucket != classification.allele_origin_bucket:
            # we have a default cc, but for the wrong allele bucket
            return ClinicalContextRecalcTrigger.CLINICAL_GROUPING_SET
    else:
        if classification.allele_object:
            # if we have no cc, but we have an allele (so we should have a cc too)
            return ClinicalContextRecalcTrigger.VARIANT_SET
    return None


def update_clinical_context(classification: Classification):
    """ Files the classification under the ClinicalContext for its allele and allele origin bucket. Agreement
        between labs is calculated by Overlaps from the groupings, so contexts are not recalculated here """
    if not _assign_new_cc_reason(classification):
        return

    existing_clinical_context = classification.clinical_context
    assign_new_clinical_context: Optional[ClinicalContext] = None
    if allele := classification.allele_object:
        assign_new_clinical_context, _ = ClinicalContext.objects.get_or_create(
            allele=allele,
            allele_origin_bucket=classification.allele_origin_bucket,
            name=existing_clinical_context.name if existing_clinical_context else ClinicalContext.default_name
        )

    classification.clinical_context = assign_new_clinical_context
    classification.save(update_fields=['clinical_context'])


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


def rehome_classifications(classifications: Iterable[Classification]):
    """ Put clinical contexts and groupings back under the Allele their classification now points at """
    classifications = list(classifications)
    if not classifications:
        return 0

    for classification in classifications:
        update_clinical_context(classification)
        ClassificationGrouping.assign_grouping_for_classification(classification)
    ClassificationGrouping.update_all_dirty()
    return len(classifications)
