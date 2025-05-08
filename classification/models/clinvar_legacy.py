from typing import Optional

from django.db import models
from django.db.models import TextChoices
from django.urls import reverse
from model_utils.models import TimeStampedModel

from classification.models import ConditionResolved
from genes.hgvs import CHGVS
from ontology.models import OntologyTerm
from snpdb.models import ClinVarKey, Allele
import re

class ClinVarLegacyAlleleMatchType(TextChoices):
    MATCHED_NOT_ATTEMPTED = "N", "Not Attempted"
    MATCHED_NOT_FOUND = "F", "Not Found"
    MATCHED_ON_CLINVAR_ALLELE = "C", "Matched on ClinVar Allele"
    MATCHED_ON_HGVS = "H", "Matched on HGVS"


CLINICAL_SIGNIFICANCE_MAP = {
    "Benign": "B",
    "Likely benign": "LB",
    "Uncertain significance": "VUS",
    "Likely pathogenic": "LP",
    "Pathogenic": "P"
}


HGVS_RE = re.compile("(?P<c_hgvs>.*?) (?P<p_hgvs>[(]p[.].*[)])")


class ClinVarLegacy(TimeStampedModel):
    class Meta:
        verbose_name = "ClinVar legacy uploads"

    def get_absolute_url(self):
        return reverse('clinvar_legacy_detail', kwargs={'scv': self.scv})

    scv = models.TextField(primary_key=True)
    clinvar_key = models.ForeignKey(ClinVarKey, on_delete=models.CASCADE)
    clinvar_variation_id = models.TextField(null=True, blank=True)
    clinvar_allele_id = models.TextField(null=True, blank=True)
    your_variant_description_hgvs = models.TextField(null=True, blank=True)

    @property
    def your_hgvs_obj(self) -> CHGVS:
        if value := self.your_variant_description_hgvs:
            return CHGVS(value)
        return None

    @property
    def hgvs_obj(self) -> CHGVS:
        if value := self.preferred_variant_name:
            if m := HGVS_RE.match(value):
                return CHGVS(m.group('c_hgvs'))
            return value
        return None

    your_variant_description_chromosome_coordinates = models.TextField(null=True, blank=True)
    preferred_variant_name = models.TextField(null=True, blank=True)
    your_condition_name = models.TextField(null=True, blank=True)
    your_condition_identifier = models.TextField(null=True, blank=True)

    @property
    def condition_obj(self) -> Optional[ConditionResolved]:
        if identifier := self.your_condition_identifier:
            return ConditionResolved(terms=[OntologyTerm.get_or_stub(identifier)])
        return None

    clinvar_condition_name = models.TextField(null=True, blank=True)
    clinical_significance = models.TextField(null=True, blank=True)

    @property
    def classification_code(self):
        return CLINICAL_SIGNIFICANCE_MAP.get(self.clinical_significance)

    date_last_evaluated = models.DateField(null=True, blank=True)
    assertion_criteria = models.TextField(null=True, blank=True)
    submitted_gene = models.TextField(null=True, blank=True)

    allele_match_date = models.DateTimeField(null=True, blank=True)
    allele = models.ForeignKey(Allele, on_delete=models.CASCADE, null=True, blank=True)
    allele_match_type = models.CharField(max_length=1, choices=ClinVarLegacyAlleleMatchType.choices, default=ClinVarLegacyAlleleMatchType.MATCHED_NOT_ATTEMPTED)

    @property
    def has_match_attempted(self) -> bool:
        return not self.allele_match_type == ClinVarLegacyAlleleMatchType.MATCHED_NOT_ATTEMPTED

    # processed_date = models.DateTimeField(null=True, blank=True)


class ClinVarLegacyMatch(TimeStampedModel):
    clinvar_legacy = models.ForeignKey(ClinVarLegacy, on_delete=models.CASCADE)
    clinvar_export = models.ForeignKey('ClinVarExport', on_delete=models.CASCADE)
    matched_allele = models.BooleanField(default=False)
    matched_condition = models.BooleanField(default=False)
    matched_scv = models.BooleanField(default=False)


    """
    VariationID = "VariationID"
    AlleleID = "AlleleID"
    Your_record_id = "Your_record_id"
    SCV = "SCV"
    RCV = "RCV"
    Your_variant_description_HGVS = "Your_variant_description_HGVS"
    Your_variant_description_chromosome_coordinates = "Your_variant_description_chromosome_coordinates"
    Preferred_variant_name = "Preferred_variant_name"
    Your_condition_name = "Your_condition_name"
    Your_condition_identifier = "Your_condition_identifier"
    ClinVar_condition_name = "ClinVar_condition_name"
    Clinical_significance = "Clinical_significance"
    Date_last_evaluated = "Date_last_evaluated"
    Assertion_criteria = "Assertion_criteria"
    Submitted_gene = "Submitted_gene"
    """