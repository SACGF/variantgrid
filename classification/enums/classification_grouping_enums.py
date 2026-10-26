"""
ClassificationGrouping vocabulary: ClassificationClassificationBucket (which clinical significance
bucket a grouping falls into), the pathogenic / somatic difference levels between a grouping's
classifications, and ClassificationGroupingSearchTermType (the kinds of term a grouping is indexed by).
"""
from typing import Optional

from django.db.models import IntegerChoices, TextChoices


class ClassificationClassificationBucket(TextChoices):
    BENIGN = "B", "Benign"
    VUS = "V", "VUS"
    PATHOGENIC = "P", "Pathogenic"
    ONCOGENIC = "O", "Oncogenic"  # should oncogenic just be an alias for pathogenic?
    OTHER = "X", "Other"
    CONFLICTING = "C", "Conflicting"
    NO_DATA = "N", "No Data"

    @staticmethod
    def bucket_for_classification(value: Optional[str]):
        if value is None:
            return ClassificationClassificationBucket.NO_DATA
        if value.startswith("VUS"):
            return ClassificationClassificationBucket.VUS
        if value in {"B", "LB"}:
            return ClassificationClassificationBucket.BENIGN
        if value in {"P", "LP"}:
            return ClassificationClassificationBucket.PATHOGENIC
        if value in {"O", "LO"}:
            return ClassificationClassificationBucket.ONCOGENIC
        return ClassificationClassificationBucket.OTHER

        # this goes a little against the buckets that we store directly into


class ClassificationGroupingPathogenicDifference(IntegerChoices):
    NO_DIFF = 0, "No Differences"
    SMALL_DIFF = 1, "Differences"
    CLIN_SIG_DIFFS = 2, "Clinical Significance Difference"


class ClassificationGroupingSomaticDifference(IntegerChoices):
    NO_DIFF = 0, "No Differences"
    AMP_DIFF = 1, "Same Tier, different AMP Level"
    TIER_DIFF = 2, "Different Tier"


class ClassificationGroupingSearchTermType(TextChoices):
    # CLASSIFICATION_ID = "CR_ID", "Classification Record ID"
    # CLASSIFICATION_LAB_ID = "CR_LAB_ID", "Classification Lab Record ID"
    # SOURCE_ID = "SOURCE_ID", "Classification Source ID"
    CONDITION_ID = "CON_ID", "Condition ID"
    CLINVAR_SCV = "SCV", "Clinvar SCV"
    GENE_SYMBOL = "GENE_SYMBOL", "Gene Symbol"
    DISCORDANCE_REPORT = "DR", "Discordance Report"

    @property
    def is_partial_text(self):
        return False
