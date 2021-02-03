from enum import Enum

from django.core.exceptions import ObjectDoesNotExist
from django.db import models

from library.utils import Constant


class ImportSource(models.TextChoices):
    """ Keeps track of where uploaded files came from """

    API = 'A', 'API'
    COMMAND_LINE = 'C', 'Command Line'
    SEQAUTO = 'S', 'SeqAuto'
    WEB = 'W', 'Web'
    WEB_UPLOAD = 'U', 'Web Upload'  # files put in settings.PRIVATE_DATA_ROOT


class ImportStatus(models.TextChoices):
    CREATED = 'C', 'created'
    IMPORTING = 'I', 'importing'
    REQUIRES_USER_INPUT = 'R', 'Requires user input'
    ERROR = 'E', 'error'
    SUCCESS = 'S', 'success'
    MARKED_FOR_DELETION = 'M', "Marked For Deletion"
    DELETING = 'D', 'Deleting'

    DELETION_STATES = Constant([e[0] for e in (MARKED_FOR_DELETION, DELETING)])


class ProcessingStatus(models.TextChoices):
    CREATED = 'C', 'Created'
    PROCESSING = 'P', 'Processing'
    ERROR = 'E', 'Error'
    SUCCESS = 'S', 'Success'
    SKIPPED = 'K', 'Skipped'
    TERMINATED_EARLY = 'T', 'Terminated Early'
    TIMED_OUT = 'Z', 'Timed Out'

    FINISHED_STATES = Constant([e[0] for e in (ERROR, SUCCESS, SKIPPED, TERMINATED_EARLY, TIMED_OUT)])


class AnnotationLevel(models.TextChoices):
    TRANSCRIPT = 'T', 'Transcript'
    GENE = 'G', 'Gene Symbol'


class ColumnAnnotationLevel(models.TextChoices):
    CLINVAR_LEVEL = 'C', 'ClinVar'
    DATABASE_LEVEL = 'D', 'Database'
    GENE_LEVEL = 'G', 'Gene'
    HGNC_LEVEL = 'H', 'HGNC'
    SAMPLE_LEVEL = 'S', 'Sample'
    TRANSCRIPT_LEVEL = 'T', 'Transcript'
    UNIPROT_LEVEL = 'U', 'UniProt'
    VARIANT_LEVEL = 'V', 'Variant'


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


class VariantsType(models.TextChoices):
    UNKNOWN = 'U', 'Unknown'
    GERMLINE = 'G', 'Germline'
    MIXED = 'M', "Mixed (Single Sample)"
    SOMATIC_ONLY = 'S', "Somatic only (Tumor minus normal)"

    SOMATIC_TYPES = Constant([e[0] for e in (MIXED, SOMATIC_ONLY)])


class SequenceRole(models.TextChoices):
    ASSEMBLED_MOLECULE = 'AM', "assembled-molecule"
    UNLOCALIZED_SCAFFOLD = 'ULS', "unlocalized-scaffold"
    UNPLACED_SCAFFOLD = 'UPS', "unplaced-scaffold"
    ALT_SCAFFOLD = 'ALT', "alt-scaffold"
    FIX_PATCH = 'FP', "fix-patch"
    NOVEL_PATCH = 'NP', "novel-patch"


class AssemblyMoleculeType(models.TextChoices):
    CHROMOSOME = 'C', "Chromosome"
    MITOCHONDRION = 'M', "Mitochondrion"


class AlleleOrigin(models.TextChoices):
    IMPORTED_TO_DATABASE = 'D', 'Imported as this build'
    IMPORTED_NORMALIZED = 'N', 'Imported (normalized)'
    LIFTOVER = 'L', 'Liftover'
    LIFTOVER_NORMALIZED = 'M', 'Liftover (normalized)'  # This probably shouldn't happen!

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


class AlleleConversionTool(models.TextChoices):
    CLINGEN_ALLELE_REGISTRY = 'CA', "ClinGen Allele Registry"
    DBSNP = 'DB', "dbSNP API"
    NCBI_REMAP = 'NR', "NCBI Remap"

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
class VCFInfoTypes(models.TextChoices):
    INTEGER = 'I', 'Integer'
    FLOAT = 'F', 'Float'
    FLAG = 'B', 'Flag'  # B for Boolean
    CHARACTER = 'C', 'Character'
    STRING = 'S', 'String'


class SuperPopulationCode(models.TextChoices):
    """ https://www.internationalgenome.org/category/population/ """
    AFR = "A", "African"
    AMR = "M", "Ad Mixed American"
    EAS = "E", "East Asian"
    EUR = "U", "European"
    SAS = "S", "South Asian"
