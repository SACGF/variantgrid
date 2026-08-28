"""
    Derives ClinVar's queryable classification values from the raw VCF text stored on ClinVar rows.

    ClinVar ships three independent classification axes - germline (CLNSIG), somatic clinical impact
    (SCI, the AMP tiers) and oncogenicity (ONC). The import and the backfill command both call this,
    so one implementation defines what a ClinVar row means.
"""
from typing import Optional

from annotation.models.models_enums import ClinVarOncogenicity, ClinVarPathogenicity
from classification.enums import SomaticClinicalSignificance

# Ordered low -> high, matched as substrings so combined forms ("Benign/Likely_benign") and
# count-suffixed conflicting values ("Pathogenic(2)|Uncertain_significance(1)") resolve.
# Keys are case-sensitive, which is what keeps "Benign" clear of "Likely_benign".
CLINSIG_TO_PATHOGENICITY = {
    "Benign": ClinVarPathogenicity.BENIGN,
    "Likely_benign": ClinVarPathogenicity.LIKELY_BENIGN,
    "Uncertain_significance": ClinVarPathogenicity.UNCERTAIN,
    "Likely_pathogenic": ClinVarPathogenicity.LIKELY_PATHOGENIC,
    "Pathogenic": ClinVarPathogenicity.PATHOGENIC,
}

ONC_TO_ONCOGENICITY = {
    "Benign": ClinVarOncogenicity.BENIGN,
    "Likely_benign": ClinVarOncogenicity.LIKELY_BENIGN,
    "Uncertain_significance": ClinVarOncogenicity.UNCERTAIN,
    "Likely_oncogenic": ClinVarOncogenicity.LIKELY_ONCOGENIC,
    "Oncogenic": ClinVarOncogenicity.ONCOGENIC,
}

# ClinVar's marker for "the classification belongs to a haplotype, not this variant" - it appears on
# CLNSIG, ONC and SCI alike, with the real call reported through CLNSIGINCL / ONCINCL / SCIINCL
NO_CLASSIFICATION_FOR_THE_SINGLE_VARIANT = "no_classification_for_the_single_variant"

SCI_TO_SOMATIC_TIER = {
    "Tier_I_-_Strong": SomaticClinicalSignificance.TIER_1,
    "Tier_II_-_Potential": SomaticClinicalSignificance.TIER_2,
    "Tier_III_-_Unknown": SomaticClinicalSignificance.TIER_3,
    "Tier_IV_-_Benign/Likely_benign": SomaticClinicalSignificance.TIER_4,
}


def _highest(value: Optional[str], conflicting_value: Optional[str], ordered_map: dict[str, int]) -> Optional[int]:
    """ Walks ordered_map low -> high keeping the last substring hit. Returns None when the aggregate value
        is conflicting but ClinVar didn't ship the breakdown, so callers can count that case """
    if not value:
        return 0

    if value.startswith("Conflicting_"):
        text = conflicting_value
    else:
        text = value

    if not text:
        return None

    highest = 0
    for aggregate_value, code in ordered_map.items():  # low->high
        if aggregate_value in text:
            highest = code
    return highest


def highest_pathogenicity(clinical_significance: Optional[str],
                          conflicting_clinical_significance: Optional[str]) -> Optional[int]:
    return _highest(clinical_significance, conflicting_clinical_significance, CLINSIG_TO_PATHOGENICITY)


def highest_oncogenicity(oncogenic_classification: Optional[str],
                         oncogenic_conflicting_classification: Optional[str]) -> Optional[int]:
    return _highest(oncogenic_classification, oncogenic_conflicting_classification, ONC_TO_ONCOGENICITY)


def somatic_tier(somatic_clinical_significance: Optional[str]) -> Optional[str]:
    """ None for anything outside ClinVar's 4 tier strings, which covers the
        'no_classification_for_the_single_variant' sentinel """
    if not somatic_clinical_significance:
        return None
    return SCI_TO_SOMATIC_TIER.get(somatic_clinical_significance)
