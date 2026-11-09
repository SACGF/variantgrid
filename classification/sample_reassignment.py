"""
Moving a classification to another sample - the one way a record changes sample (#1392).

Re-populating runs the sample autopopulate (get_evidence_fields_for_sample_and_patient, the only writer of those
keys) for both samples: what the new sample supplies is overwritten, and what only the old sample supplied is cleared,
so nothing from the old sample survives and keys neither sample supplies are left alone. Allele origin is never
cleared: a sample that asserts none (mixed tumour) leaves the per-variant origin to the curator. The patch is
SubmissionSource.VARIANT_GRID, so it overwrites values a curator typed - the preview is where they see what goes.

Entry points: sample_reassignment_preview (what would change) and reassign_sample (do it), both called from
classification/views/views_classification_sample.py.
"""
from dataclasses import dataclass
from typing import Any, Optional

from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.db import transaction

from classification.autopopulate_evidence_keys.evidence_from_sample_and_patient import (
    get_evidence_fields_for_sample_and_patient,
)
from classification.enums import SpecialEKeys, SubmissionSource
from classification.models.classification import Classification
from classification.models.classification_utils import ClassificationPatchResponse
from classification.models.evidence_key import EvidenceKeyMap
from snpdb.enums import AlleleOriginBucket
from snpdb.models.models_vcf import Sample


@dataclass(frozen=True)
class SampleEvidenceChange:
    key: str
    label: str
    current: Any  # pretty value on the record now, None when empty
    new: Any  # pretty value the new sample supplies, None when it will be cleared


@dataclass(frozen=True)
class SampleReassignmentPreview:
    sample: Sample
    has_genotype: bool
    changes: list[SampleEvidenceChange]  # only the keys whose value would change
    allele_origin_bucket_before: AlleleOriginBucket
    allele_origin_bucket_after: AlleleOriginBucket

    @property
    def allele_origin_bucket_changes(self) -> bool:
        return self.allele_origin_bucket_before != self.allele_origin_bucket_after


def check_can_reassign_sample(classification: Classification, sample: Sample, user: User):
    if not classification.can_reassign_sample(user):
        raise PermissionDenied(f"You may not change the sample of classification {classification.pk}")
    sample.check_can_view(user)


def sample_has_genotype(classification: Classification, sample: Sample) -> bool:
    return sample.get_genotype(classification.variant) is not None


def _sample_values(classification: Classification, sample: Optional[Sample]) -> dict[str, Any]:
    if sample is None:
        return {}
    return get_evidence_fields_for_sample_and_patient(classification.variant, sample).data


def _sample_patch(classification: Classification, sample: Sample) -> dict[str, Any]:
    """ The new sample's values, and None for each key the record's current sample supplied and the new one doesn't """
    new_values = _sample_values(classification, sample)
    patch = dict(new_values)
    for key in _sample_values(classification, classification.sample):
        if key not in new_values and key != SpecialEKeys.ALLELE_ORIGIN and classification.get(key) is not None:
            patch[key] = None
    return patch


def _apply_sample_patch(classification: Classification, patch: dict[str, Any], user: User,
                        save: bool) -> ClassificationPatchResponse:
    return classification.patch_value(patch, user=user, source=SubmissionSource.VARIANT_GRID,
                                      leave_existing_values=False, save=save)


def sample_reassignment_preview(classification: Classification, sample: Sample,
                                user: User) -> SampleReassignmentPreview:
    check_can_reassign_sample(classification, sample, user)

    # Diffed after normalisation (option labels to keys etc.) on an unsaved copy, so a value that only differs
    # in how the sample spells it isn't reported as a change
    patch = _sample_patch(classification, sample)
    patched = Classification.objects.get(pk=classification.pk)
    _apply_sample_patch(patched, patch, user, save=False)

    e_keys = EvidenceKeyMap.cached()

    def pretty(record: Classification, key: str) -> Optional[str]:
        value = record.get(key)
        return e_keys.get(key).pretty_value(value) if value is not None else None

    changes = []
    for key in patch:
        current = pretty(classification, key)
        new = pretty(patched, key)
        if current != new:
            changes.append(SampleEvidenceChange(key=key, label=e_keys.get(key).pretty_label,
                                                current=current, new=new))

    return SampleReassignmentPreview(
        sample=sample,
        has_genotype=sample_has_genotype(classification, sample),
        changes=changes,
        allele_origin_bucket_before=AlleleOriginBucket(classification.allele_origin_bucket),
        allele_origin_bucket_after=AlleleOriginBucket(patched.allele_origin_bucket),
    )


def reassign_sample(classification: Classification, sample: Sample, user: User,
                    repopulate: bool) -> Optional[ClassificationPatchResponse]:
    check_can_reassign_sample(classification, sample, user)
    with transaction.atomic():
        # computed before re-pointing, as what to clear depends on the record's current sample
        patch = _sample_patch(classification, sample) if repopulate else None
        classification.sample = sample
        classification.save(update_fields=["sample", "modified"])
        if patch is not None:
            return _apply_sample_patch(classification, patch, user, save=True)
    return None
