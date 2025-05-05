from django.db import models
from django.db.models import TextChoices
from model_utils.models import TimeStampedModel
from snpdb.models import ClinVarKey, Allele


class ClinVarLegacyAlleleMatchType(TextChoices):
    MATCHED_ON_CLINVAR_ALLELE = "C", "Matched on ClinVar Allele"
    MATCHED_ON_HGVS = "H", "Matched on HGVS"


class ClinVarLegacy(TimeStampedModel):
    scv = models.TextField(primary_key=True)
    clinvar_key = models.ForeignKey(ClinVarKey, on_delete=models.CASCADE)
    clinvar_variation_id = models.TextField(null=True, blank=True)
    clinvar_allele_id = models.TextField(null=True, blank=True)
    your_variant_description_hgvs = models.TextField(null=True, blank=True)
    your_variant_description_chromosome_coordinates = models.TextField(null=True, blank=True)
    preferred_variant_name = models.TextField(null=True, blank=True)
    your_condition_name = models.TextField(null=True, blank=True)
    your_condition_identifier = models.TextField(null=True, blank=True)
    clinvar_condition_name = models.TextField(null=True, blank=True)
    clinical_significance = models.TextField(null=True, blank=True)
    date_last_evaluated = models.DateField(null=True, blank=True)
    assertion_criteria = models.TextField(null=True, blank=True)
    submitted_gene = models.TextField(null=True, blank=True)

    allele = models.ForeignKey(Allele, on_delete=models.CASCADE, null=True, blank=True)
    allele_match_type = models.CharField(max_length=1, choices=ClinVarLegacyAlleleMatchType.choices, null=True, blank=True)

    processed_date = models.DateTimeField(null=True, blank=True)



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