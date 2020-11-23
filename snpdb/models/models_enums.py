from enum import Enum

from django.core.exceptions import ObjectDoesNotExist


class ImportSource:
    """ Keeps track of where uploaded files came from """

    API = 'A'
    COMMAND_LINE = 'C'
    SEQAUTO = 'S'
    WEB = 'W'
    WEB_UPLOAD = 'U'  # files put in settings.PRIVATE_DATA_ROOT
    CHOICES = (
        (API, 'API'),
        (COMMAND_LINE, 'Command Line'),
        (SEQAUTO, 'SeqAuto'),
        (WEB, 'Web'),
        (WEB_UPLOAD, 'Web Upload'),
    )


class ImportStatus:
    CREATED = 'C'
    IMPORTING = 'I'
    REQUIRES_USER_INPUT = 'R'
    ERROR = 'E'
    SUCCESS = 'S'
    MARKED_FOR_DELETION = 'M'
    DELETING = 'D'

    CHOICES = (
        (CREATED, 'created'),
        (IMPORTING, 'importing'),
        (REQUIRES_USER_INPUT, 'Requires user input'),
        (ERROR, 'error'),
        (SUCCESS, 'success'),
        (MARKED_FOR_DELETION, "Marked For Deletion"),
        (DELETING, 'Deleting'),
    )
    DELETION_STATES = [MARKED_FOR_DELETION, DELETING]


class ProcessingStatus:
    CREATED = 'C'
    PROCESSING = 'P'
    ERROR = 'E'
    SUCCESS = 'S'
    SKIPPED = 'K'
    TERMINATED_EARLY = 'T'
    TIMED_OUT = 'Z'
    CHOICES = (
        (CREATED, 'Created'),
        (PROCESSING, 'Processing'),
        (ERROR, 'Error'),
        (SUCCESS, 'Success'),
        (SKIPPED, 'Skipped'),
        (TERMINATED_EARLY, 'Terminated Early'),
        (TIMED_OUT, 'Timed Out')
    )

    FINISHED_STATES = [
        ERROR, SUCCESS, SKIPPED, TERMINATED_EARLY, TIMED_OUT
    ]


class AnnotationLevel:
    TRANSCRIPT = 'T'
    GENE = 'G'
    CHOICES = (
        (TRANSCRIPT, 'Transcript'),
        (GENE, 'Gene Symbol')
    )


class ColumnAnnotationLevel:
    CLINVAR_LEVEL = 'C'
    DATABASE_LEVEL = 'D'
    GENE_LEVEL = 'G'
    SAMPLE_LEVEL = 'S'
    TRANSCRIPT_LEVEL = 'T'
    VARIANT_LEVEL = 'V'

    CHOICES = [
        (CLINVAR_LEVEL, 'ClinVar'),
        (DATABASE_LEVEL, 'Database'),
        (GENE_LEVEL, 'Gene'),
        (SAMPLE_LEVEL, 'Sample'),
        (TRANSCRIPT_LEVEL, 'Transcript'),
        (VARIANT_LEVEL, 'Variant'),
    ]


class BuiltInFilters:
    TOTAL = "T"
    CLINVAR = "C"
    OMIM = "O"
    IMPACT_HIGH_OR_MODERATE = 'I'
    CLASSIFIED = 'G'  # G = for Genomic Variant Classification
    CLASSIFIED_PATHOGENIC = 'P'
    COSMIC = 'M'  # cosMic

    FILTER_CHOICES = [
        # Don't include total (as that's no filter at all!
        (CLINVAR, 'ClinVar'),
        (OMIM, 'OMIM Phenotype'),
        (IMPACT_HIGH_OR_MODERATE, 'High or Mod impact'),
        (CLASSIFIED, 'Classified'),
        (CLASSIFIED_PATHOGENIC, 'Classified Pathogenic'),
        (COSMIC, 'COSMIC')]
    CHOICES = [(TOTAL, 'Total')] + FILTER_CHOICES
    COLORS = [(TOTAL, "#000000"),
              (CLINVAR, '#ff0000'),
              (OMIM, "#99CD83"),
              (IMPACT_HIGH_OR_MODERATE, "#aaaaff"),
              (CLASSIFIED, "#7c26cb"),
              (CLASSIFIED_PATHOGENIC, "#Ff008b"),
              (COSMIC, "#14559f")]
    DEFAULT_NODE_COUNT_FILTERS = [TOTAL, IMPACT_HIGH_OR_MODERATE, CLINVAR]


class VariantsType:
    UNKNOWN = 'U'
    GERMLINE = 'G'
    MIXED = 'M'
    SOMATIC_ONLY = 'S'  # Eg Tumor/Normal subtraction

    CHOICES = (
        (UNKNOWN, 'Unknown'),
        (GERMLINE, 'Germline'),
        (MIXED, "Mixed (Single Sample)"),
        (SOMATIC_ONLY, "Somatic only (Tumor minus normal)"),
    )

    SOMATIC_TYPES = [MIXED, SOMATIC_ONLY]


class SequenceRole:
    ASSEMBLED_MOLECULE = 'AM'
    UNLOCALIZED_SCAFFOLD = 'ULS'
    UNPLACED_SCAFFOLD = 'UPS'
    ALT_SCAFFOLD = 'ALT'
    FIX_PATCH = 'FP'
    NOVEL_PATCH = 'NP'

    CHOICES = (
        (ASSEMBLED_MOLECULE, "assembled-molecule"),
        (UNLOCALIZED_SCAFFOLD, "unlocalized-scaffold"),
        (UNPLACED_SCAFFOLD, "unplaced-scaffold"),
        (ALT_SCAFFOLD, "alt-scaffold"),
        (FIX_PATCH, "fix-patch"),
        (NOVEL_PATCH, "novel-patch"),
    )


class AssemblyMoleculeType:
    CHROMOSOME = 'C'
    MITOCHONDRION = 'M'

    CHOICES = (
        (CHROMOSOME, "Chromosome"),
        (MITOCHONDRION, "Mitochondrion"),
    )


class ClinGenAlleleRegistryErrorType:
    # @see http://reg.clinicalgenome.org/doc/AlleleRegistry_1.01.xx_api_v1.pdf
    NOT_FOUND = "F"
    AUTHORIZATION_ERROR = "U"
    INCORRECT_REQUEST = "Q"
    HGVS_PARSING_ERROR = "H"
    INCORRECT_HGVS_POSITION = "P"
    INCORRECT_REFERENCE_ALLELE = "I"
    NO_CONSISTENT_ALIGNMENT = "N"
    UNKNOWN_CDS = "C"
    UNKNOWN_GENE = "G"
    UNKNOWN_REFERENCE_SEQUENCE = "R"
    VCF_PARSING_ERROR = "V"
    REQUEST_TOO_LARGE = "L"
    INTERNAL_SERVER_ERROR = "E"

    CHOICES = (
        (NOT_FOUND, "NotFound"),
        (AUTHORIZATION_ERROR, "AuthorizationError"),
        (INCORRECT_REQUEST, "IncorrectRequest"),
        (HGVS_PARSING_ERROR, "HgvsParsingError"),
        (INCORRECT_HGVS_POSITION, "IncorrectHgvsPosition"),
        (INCORRECT_REFERENCE_ALLELE, "IncorrectReferenceAllele"),
        (NO_CONSISTENT_ALIGNMENT, "NoConsistentAlignment"),
        (UNKNOWN_CDS, "UnknownCDS"),
        (UNKNOWN_GENE, "UnknownGene"),
        (UNKNOWN_REFERENCE_SEQUENCE, "UnknownReferenceSequence"),
        (VCF_PARSING_ERROR, "VcfParsingError"),
        (REQUEST_TOO_LARGE, "RequestTooLarge"),
        (INTERNAL_SERVER_ERROR, "InternalServerError")
    )


class AlleleOrigin:
    IMPORTED_TO_DATABASE = 'D'
    IMPORTED_NORMALIZED = 'N'
    LIFTOVER = 'L'
    LIFTOVER_NORMALIZED = 'M'  # This probably shouldn't happen!

    CHOICES = (
        (IMPORTED_TO_DATABASE, 'Imported as this build'),
        (IMPORTED_NORMALIZED, 'Imported (normalized)'),
        (LIFTOVER, 'Liftover'),
        (LIFTOVER_NORMALIZED, 'Liftover (normalized)'),
    )

    @staticmethod
    def variant_origin(variant):
        try:
            origin = variant.variantallele.origin
        except ObjectDoesNotExist:
            # Variant w/o allele must have been imported directly
            if variant.modifiedimportedvariant_set.exists():
                origin = AlleleOrigin.IMPORTED_NORMALIZED
            else:
                origin = AlleleOrigin.IMPORTED_TO_DATABASE
        return origin


class AlleleConversionTool:
    CLINGEN_ALLELE_REGISTRY = 'CA'
    DBSNP = 'DB'
    NCBI_REMAP = 'NR'

    CHOICES = (
        (CLINGEN_ALLELE_REGISTRY, "ClinGen Allele Registry"),
        (DBSNP, "dbSNP API"),
        (NCBI_REMAP, "NCBI Remap")
    )

    @classmethod
    def vcf_tuples_in_destination_build(cls, conversion_tool):
        IN_DEST_BUILD = {cls.CLINGEN_ALLELE_REGISTRY: True,
                         cls.DBSNP: True,
                         cls.NCBI_REMAP: False}
        return IN_DEST_BUILD[conversion_tool]


class ClinGenAlleleExternalRecordType(Enum):
    """ @see http://reg.clinicalgenome.org/doc/AlleleRegistry_1.01.xx_api_v1.pdf
        "Query canonical alleles by identifiers from external records" section """
    DBSNP_ID = "dbSNP.rs"


# Integer, Float, Flag, Character, and String.
class VCFInfoTypes:
    INTEGER = 'I'
    FLOAT = 'F'
    FLAG = 'B'  # Boolean
    CHARACTER = 'C'
    STRING = 'S'
    CHOICES = [
        (INTEGER, 'Integer'),
        (FLOAT, 'Float'),
        (FLAG, 'Flag'),
        (CHARACTER, 'Character'),
        (STRING, 'String'),
    ]