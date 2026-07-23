""" The two Beacon datasets served as independent resultSets (§5.5):
    - variantgrid_observations   : sample-genotype presence/count (§5.2), permission-scoped
    - variantgrid_classifications: published ShareLevel.PUBLIC classifications (the MME set)

Each builder resolves presence/count/records for one dataset and returns a DatasetResult;
the view stitches them into the response envelope and clamps to the returned granularity.
"""
from dataclasses import dataclass, field
from typing import Optional

from django.conf import settings

from beacon.schema import OBSERVATIONS_DATASET_ID, CLASSIFICATIONS_DATASET_ID
from beacon.variant_mapping import variant_to_g_variant
from classification.enums.classification_enums import SpecialEKeys
from classification.models.evidence_key import EvidenceKeyMap
from mme.contact import lab_contact
from mme.serializers.patient_profile import (
    mme_eligible_classifications,
    classification_ontology_slots,
)
from patients.models import Patient
from snpdb.models import Variant, GenomeBuild, Allele
from snpdb.models.models_zygosity_counts import VariantZygosityCountCollection, VariantZygosityCount
from snpdb.variant_sample_information import VariantZygosityCounts


@dataclass
class DatasetResult:
    """ One dataset's contribution to a g_variants response.
        `reportable_count` is None when a k-anonymity floor suppresses the exact count -
        the resultSet then exposes count 0 with no records (presence without the small
        count), rather than the true value. """
    dataset_id: str
    exists: bool = False
    count: int = 0                        # true count (used for audit + aggregation)
    reportable_count: Optional[int] = 0   # count to expose; None = suppressed
    records: list = field(default_factory=list)  # record-tier detail

    def result_set(self) -> dict:
        # `resultsCount` and `results` are always present: the Beacon resultSet schema
        # requires both, and spec clients (e.g. the EGA beacon-verifier) index `results`
        # unconditionally when exists is true. A floor-suppressed count is exposed as 0.
        return {
            "id": self.dataset_id,
            "setType": "genomicVariant",
            "exists": self.exists,
            "resultsCount": self.reportable_count if self.reportable_count is not None else 0,
            "results": self.records,
        }


def _global_germline_present(variant: Variant) -> Optional[bool]:
    """ Stage-1 negative gate (§5.2): read the precomputed global VariantZygosityCount.
        Returns False when the global het+hom count is 0 (cheap proof of absence, so we
        can skip the CohortGenotype walk), True when present, or None when the global
        collection is unavailable (fall through to the stage-2 permission walk). """
    try:
        vzcc = VariantZygosityCountCollection.get_global_germline_counts()
    except VariantZygosityCountCollection.DoesNotExist:
        return None
    vzc = VariantZygosityCount.objects.filter(collection=vzcc, variant=variant).first()
    if vzc is None:
        return False
    return (vzc.het_count + vzc.hom_count) > 0


def observations_dataset(user, variant: Optional[Variant], genome_build: GenomeBuild) -> DatasetResult:
    """ variantgrid_observations: two-stage lookup (§5.2). Global VZC negative gate, then a
        permission-scoped CohortGenotype count via the extracted VariantZygosityCounts core.
        Exact counts below settings.BEACON_MIN_REPORTABLE_COUNT are suppressed (resultSet
        drops to boolean presence) to limit membership-inference risk. """
    result = DatasetResult(dataset_id=OBSERVATIONS_DATASET_ID)
    if variant is None:
        return result

    # Stage 1 - cheap negative gate. False => provably absent, done.
    if _global_germline_present(variant) is False:
        return result

    # Stage 2 - permission-scoped exact count (anonymous -> public group).
    counts = VariantZygosityCounts(user, variant, genome_build)
    visible = counts.num_visible_alt
    result.count = visible
    result.exists = visible > 0
    if not result.exists:
        result.reportable_count = 0
        return result

    floor = settings.BEACON_MIN_REPORTABLE_COUNT
    if visible < floor:
        # Suppress the exact small count: keep exists=True but drop to boolean presence.
        result.reportable_count = None
        return result

    result.reportable_count = visible
    record = variant_to_g_variant(variant, genome_build)
    record["variantGenotypeCounts"] = {
        "heterozygousCount": counts.num_visible_het,
        "homozygousCount": counts.num_visible_hom_alt,
    }
    result.records = [record]
    return result


def _patient_is_public(patient) -> bool:
    """ §5.5 extra gate: is the linked Patient itself shared publicly (patient-level
        Guardian permission -> public group)? filter_for_user(None) resolves to the
        public group's readable patients. """
    if patient is None:
        return False
    return Patient.filter_for_user(None).filter(pk=patient.pk).exists()


def _classification_record(modification) -> dict:
    """ Record-tier payload for one public classification (§5.5): ACMG significance,
        condition/phenotype, owning-lab contact. Phenotype includes the linked patient's
        HPO terms only when that patient is shared publicly. """
    classification = modification.classification
    cs = modification.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
    cs_label = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).pretty_value(cs) if cs else None

    patient = getattr(classification.sample, "patient", None) if classification.sample else None
    include_patient = _patient_is_public(patient)
    features, disorders = classification_ontology_slots(classification, include_patient_phenotype=include_patient)

    record = {
        "clinicalInterpretations": [{
            "clinicalRelevance": cs_label,
            "conditions": disorders,     # OMIM/Orphanet disorder ids
            "phenotypes": features,      # HPO features (condition + optionally public patient)
        }],
    }
    if contact := lab_contact(classification):
        record["contact"] = contact
    return record


def classifications_dataset(user, allele: Optional[Allele], granularity_is_record: bool) -> DatasetResult:
    """ variantgrid_classifications: published + ShareLevel.PUBLIC classifications for the
        allele (the same consented set MME submits/serves). No k-anonymity floor - each
        record is already a deliberate record-level public share. """
    result = DatasetResult(dataset_id=CLASSIFICATIONS_DATASET_ID)
    if allele is None:
        return result

    qs = mme_eligible_classifications().filter(classification__allele=allele)
    count = qs.count()
    result.count = count
    result.reportable_count = count
    result.exists = count > 0
    if result.exists and granularity_is_record:
        result.records = [_classification_record(m) for m in qs]
    return result
