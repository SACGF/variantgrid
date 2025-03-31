from django.db.models import TextChoices


class ClinVarExportTypeBucket(TextChoices):
    GERMLINE = "G", "Germline"
    ONCOGENIC = "O", "Oncogenic"
    CLINICAL_IMPACT = "C", "Clinical Impact"