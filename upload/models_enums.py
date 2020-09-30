class UploadedFileTypes:
    BED = 'B'
    CLINVAR = 'L'
    CLINVAR_CITATIONS = 'T'
    CUFFDIFF = 'C'
    GENE_LIST = 'G'
    GENE_COVERAGE = 'O'
    LIFTOVER = 'I'
    PATIENT_RECORDS = 'R'
    PED = 'P'
    VCF = 'V'
    VCF_INSERT_VARIANTS_ONLY = 'Y'
    VARIANT_CLASSIFICATIONS = 'S'

    CHOICES = (
        (BED, 'BED'),
        (CLINVAR, 'Clinvar'),
        (CLINVAR_CITATIONS, 'Clinvar Citations'),
        (CUFFDIFF, 'CuffDiff'),
        (GENE_LIST, 'Gene List'),
        (GENE_COVERAGE, 'Gene Coverage'),
        (LIFTOVER, 'Liftover'),
        (PED, 'Pedigree'),
        (PATIENT_RECORDS, 'Patient Records'),
        (VCF, 'VCF'),
        (VCF_INSERT_VARIANTS_ONLY, 'VCF - Insert variants only (no samples etc)'),
        (VARIANT_CLASSIFICATIONS, 'Variant Classifications'),
    )


class UploadStepTaskType:
    CELERY = 'C'
    SQL = 'Q'
    TOOL = 'T'
    CHOICES = (
        (CELERY, 'Celery'),
        (SQL, 'SQL'),
        (TOOL, 'Tool'),
    )


class VCFPipelineStage:
    """ Some jobs can only run when we've moved to a certain stage of the pipeline """

    INSERT_UNKNOWN_VARIANTS = 'U'
    DATA_INSERTION = 'D'
    ANNOTATION_COMPLETE = 'A'
    FINISH = 'F'

    CHOICES = (
        (INSERT_UNKNOWN_VARIANTS, 'Insert Unknown Variants'),
        (DATA_INSERTION, 'Data Insertion'),
        (ANNOTATION_COMPLETE, 'Annotation Complete'),
        (FINISH, 'Finish'),
    )


class ExpressionType:
    CUFFDIFF = 'C'
    EDGE_R = 'E'
    CHOICES = (
        (CUFFDIFF, 'CuffDiff'),
        #(EDGE_R, 'EdgeR'),
    )


class TimeFilterMethod:
    DAYS = 'D'
    RECORDS = 'R'
    CHOICES = (
        (DAYS, "days"),
        (RECORDS, "records"),
    )


class VCFImportInfoSeverity:
    WARNING = 'W'
    ERROR = 'E'

    CHOICES = [(WARNING, 'WARNING'),
               (ERROR, 'ERROR')]


class UploadStepOrigin:
    USER_ADDITION = 'A'
    IMPORT_TASK_FACTORY = 'I'

    CHOICES = [
        (USER_ADDITION, "User Addition"),
        (IMPORT_TASK_FACTORY, "Import Task Factory"),
    ]
