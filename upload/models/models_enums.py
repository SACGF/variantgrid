from django.db import models


class UploadedFileTypes(models.TextChoices):
    BED = 'B', 'BED'
    CLINVAR = 'L', 'Clinvar'
    CUFFDIFF = 'C', 'CuffDiff'
    GENE_LIST = 'G', 'Gene List'
    GENE_COVERAGE = 'O', 'Gene Coverage'
    LIFTOVER = 'I', 'Liftover'
    PED = 'P', 'Pedigree'
    PATIENT_RECORDS = 'R', 'Patient Records'
    VARIANT_CLASSIFICATIONS = 'S', 'Variant Classifications'
    VARIANT_TAGS = "A", "Variant Tags"
    VCF = 'V', 'VCF'
    VCF_INSERT_VARIANTS_ONLY = 'Y', 'VCF - Insert variants only (no samples etc)'


class UploadStepTaskType(models.TextChoices):
    CELERY = 'C', 'Celery'
    SQL = 'Q', 'SQL'
    TOOL = 'T', 'Tool'


class VCFPipelineStage(models.TextChoices):
    """ Some jobs can only run when we've moved to a certain stage of the pipeline """

    INSERT_UNKNOWN_VARIANTS = 'U', 'Insert Unknown Variants'
    DATA_INSERTION = 'D', 'Data Insertion'
    ANNOTATION_COMPLETE = 'A', 'Annotation Complete'
    FINISH = 'F', 'Finish'


class ExpressionType(models.TextChoices):
    CUFFDIFF = 'C', 'CuffDiff'
    # EDGE_R = 'E', 'EdgeR'


class TimeFilterMethod(models.TextChoices):
    DAYS = 'D', "days"
    RECORDS = 'R', "records"


class VCFImportInfoSeverity(models.TextChoices):
    WARNING = 'W', 'WARNING'
    ERROR = 'E', 'ERROR'


class UploadStepOrigin(models.TextChoices):
    USER_ADDITION = 'A', "User Addition"
    IMPORT_TASK_FACTORY = 'I', "Import Task Factory"
