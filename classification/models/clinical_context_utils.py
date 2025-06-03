from collections import defaultdict
from dataclasses import dataclass
from typing import List, Optional, Iterable

from classification.models import Classification, ClinicalContextRecalcTrigger, ClinicalContext
from library.utils import first


@dataclass(frozen=True)
class DirtyClinicalContext:
    clinical_context: ClinicalContext
    cause_code: ClinicalContextRecalcTrigger
    cause: str

    def recalc_and_save(self):
        self.clinical_context.recalc_and_save(cause_code=self.cause_code, cause=self.cause)

    @staticmethod
    def merge(dirty_contexts: ['DirtyClinicalContext']) -> ['DirtyClinicalContext']:
        if len(dirty_contexts) <= 1:
            return dirty_contexts

        by_clinical_contexts = defaultdict(list)
        dc: DirtyClinicalContext
        for dc in dirty_contexts:
            by_clinical_contexts[dc.clinical_context].append(dc)

        flattened_dc: List[DirtyClinicalContext] = []
        for dcc_list in by_clinical_contexts.values():
            if len(dcc_list) == 1:
                flattened_dc.append(dcc_list[0])
            else:
                first_dcc: DirtyClinicalContext = first(dcc_list)
                flattened_dc.append(DirtyClinicalContext(
                    clinical_context=first_dcc.clinical_context,
                    cause_code=first_dcc.cause_code,
                    cause="\n".join([dcc.cause for dcc in dcc_list]),
                ))
        return flattened_dc


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


def _update_clinical_context(
        classification: Classification,
        force_recalc_text: Optional[str] = None) -> List[DirtyClinicalContext]:
    """
    :param classification: The classification to check the clinical context for
    :param force_recalc_text: If this has a value, and the classification is assigned to an allele, it's
    clinical context will recalc_and_save regardless of if the clinical context changed or not, and this
    text will be used as the trigger.
    :return: None
    """

    assign_new_clinical_trigger = _assign_new_cc_reason(classification)
    existing_clinical_context = classification.clinical_context

    if not assign_new_clinical_trigger:
        if force_recalc_text and existing_clinical_context:
            # the current clinical context is just fine, but we've been asked to recalculate it regardless
            # e.g., a classification record has changed classification and was re-published
            existing_clinical_context.recalc_and_save(cause=force_recalc_text, cause_code=ClinicalContextRecalcTrigger.SUBMISSION)
        return []

    contexts_to_recalculate: List[ClinicalContext] = []

    if existing_clinical_context:
        # since we're changing the clinical context, we're going to need to call recalc on the old one after
        # we move over to the new one
        contexts_to_recalculate.append(existing_clinical_context)

    assign_new_clinical_context: Optional[ClinicalContext] = None
    if allele := classification.allele_object:
        assign_new_clinical_context, _ = ClinicalContext.objects.get_or_create(
            allele=allele,
            allele_origin_bucket=classification.allele_origin_bucket,
            name=existing_clinical_context.name if existing_clinical_context else ClinicalContext.default_name
        )
        contexts_to_recalculate.append(assign_new_clinical_context)

    # just make sure the new clinical context doesn't match the old one (should never happen)
    if classification.clinical_context == assign_new_clinical_context:
        raise ValueError(f"Code determined there needed to be a new clinical context but recalculated the old one, trigger was {assign_new_clinical_trigger}, context was {assign_new_clinical_context}")

    # actually assign the new clinical context, and then recalc all affected clinical contexts
    classification.clinical_context = assign_new_clinical_context
    classification.save(update_fields=['clinical_context'])

    dirty_contexts = []

    for needs_recalc in contexts_to_recalculate:
        notes = f"{classification.cr_lab_id}"
        if assign_new_clinical_context == ClinicalContextRecalcTrigger.CLINICAL_GROUPING_SET:
            notes += " changes between germline/somatic"
        elif assign_new_clinical_context == ClinicalContextRecalcTrigger.VARIANT_SET:
            if classification.allele_object:
                notes += " matched to this allele"
            else:
                notes += " unmatched from the allele"
        elif assign_new_clinical_context == ClinicalContextRecalcTrigger.SUBMISSION:
            notes += " submitted"

        dirty_contexts.append(
            DirtyClinicalContext(
                clinical_context=needs_recalc,
                cause=force_recalc_text or notes,
                cause_code=assign_new_clinical_trigger
            )
        )
    return dirty_contexts


def update_clinical_contexts(
        classifications: Iterable[Classification],
        force_recalc_text: Optional[str] = None):

    all_dcs: List[DirtyClinicalContext] = []
    for classification in classifications:
        all_dcs += _update_clinical_context(classification, force_recalc_text)

    all_dcs = DirtyClinicalContext.merge(all_dcs)
    for dc in all_dcs:
        dc.recalc_and_save()


def update_clinical_context(
        classification: Classification,
        force_recalc_text: Optional[str] = None):
    """
    Convenience method for only recalculating a single classification's clinical context
    """
    update_clinical_contexts([classification], force_recalc_text)
