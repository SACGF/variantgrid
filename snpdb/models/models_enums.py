from enum import Enum
from typing import Optional

from django.core.exceptions import ObjectDoesNotExist
from django.db import models

from classification.enums import AlleleOriginBucket
from library.utils import Constant


class UserAwardLevel(models.TextChoices):
    GOLD = "G", "Gold"
    SILVER = "S", "Silver"
    BRONZE = "B", "Bronze"

    @property
    def int_value(self) -> int:
        if self == UserAwardLevel.GOLD:
            return 3
        elif self == UserAwardLevel.SILVER:
            return 2
        elif self == UserAwardLevel.BRONZE:
            return 1

    def __lt__(self, other):
        return self.int_value < other.int_value


class UserAwardKind(models.TextChoices):
    TITLE = "T", "Title"  # Held until someone takes it - crown/medal/trophy
    BADGE = "B", "Badge"  # Permanent once earned, tiered bronze/silver/gold
    KUDOS = "K", "Kudos"  # Hand-given by an admin


class AwardPeriod(models.TextChoices):
    ALL_TIME = "A", "all time"
    MONTH = "M", "this month"
    DAY = "D", "today"

    @property
    def icon(self) -> str:
        return {
            AwardPeriod.ALL_TIME: "fa-crown",
            AwardPeriod.MONTH: "fa-medal",
            AwardPeriod.DAY: "fa-trophy",
        }[self]

    @property
    def rank(self) -> int:
        """ ALL_TIME outranks MONTH outranks DAY """
        return {AwardPeriod.ALL_TIME: 3, AwardPeriod.MONTH: 2, AwardPeriod.DAY: 1}[self]


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

    RUNNING_STATES = Constant([e[0] for e in (CREATED, PROCESSING)])
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
    CLASSIFIED_PATHOGENIC = 'P'  # Germline LP/P
    CLASSIFIED_TIER_1_2 = 'S'  # Somatic Tier I/II
    COSMIC = 'M'  # cosMic

    FILTER_CHOICES = [
        # Don't include total (as that's no filter at all!
        (CLINVAR, 'ClinVar significant'),
        (OMIM, 'OMIM Phenotype'),
        (IMPACT_HIGH_OR_MODERATE, 'High or Mod impact'),
        (CLASSIFIED, 'Classified'),
        (CLASSIFIED_PATHOGENIC, 'Classified Pathogenic'),
        (CLASSIFIED_TIER_1_2, 'Classified Tier I/II'),
        (COSMIC, 'COSMIC')]
    CHOICES = [(TOTAL, 'Total')] + FILTER_CHOICES
    COLORS = [(TOTAL, "#000000"),
              (CLINVAR, '#ff0000'),
              (OMIM, "#99CD83"),
              (IMPACT_HIGH_OR_MODERATE, "#aaaaff"),
              (CLASSIFIED, "#7c26cb"),
              (CLASSIFIED_PATHOGENIC, "#Ff008b"),
              (CLASSIFIED_TIER_1_2, "#cc99cc"),
              (COSMIC, "#14559f")]
    DEFAULT_NODE_COUNT_FILTERS = [TOTAL, IMPACT_HIGH_OR_MODERATE, CLINVAR]


class TagFilter:
    """ Per-tag node counts / grid filters, eg 'tag_Artifact'. These sit alongside BuiltInFilters
        wherever a node count label or 'extra_filters' is used - tag names are alphanumeric so the
        prefix and separator can never collide with a tag, and the label stays a valid URL slug and
        CSS class. A node count is always a single tag, an 'extra_filters' selection can be several
        (comma joined, meaning a variant carrying any of them) """
    PREFIX = "tag_"
    SEPARATOR = ","

    @staticmethod
    def label(tag_id: str) -> str:
        return TagFilter.PREFIX + tag_id

    @staticmethod
    def label_for_tags(tag_ids: list[str]) -> str:
        return TagFilter.SEPARATOR.join(TagFilter.label(tag_id) for tag_id in tag_ids)

    @staticmethod
    def get_tag_id(label: str) -> Optional[str]:
        """ Tag id for a tag node count label, or None if it's not one """
        if label and label.startswith(TagFilter.PREFIX):
            return label[len(TagFilter.PREFIX):]
        return None

    @staticmethod
    def get_tag_ids(label: str) -> list[str]:
        """ Tag ids for a (possibly multi-tag) label, or [] if it's not one """
        tag_ids = []
        for part in (label or "").split(TagFilter.SEPARATOR):
            tag_id = TagFilter.get_tag_id(part)
            if not tag_id:
                return []
            tag_ids.append(tag_id)
        return tag_ids


class VariantsType(models.TextChoices):
    UNKNOWN = 'U', 'Unknown'
    GERMLINE = 'G', 'Germline'
    MIXED = 'M', "Mixed (Single Sample)"
    SOMATIC_ONLY = 'S', "Somatic only (Tumour minus normal)"

    SOMATIC_TYPES = Constant([e[0] for e in (MIXED, SOMATIC_ONLY)])


class SampleFileType(models.TextChoices):
    """ A file linked to a sample """
    BAM = 'B', 'BAM'
    CRAM = 'C', 'CRAM'
    BED = 'E', 'BED'
    VCF = 'V', 'VCF'


class SequenceRole(models.TextChoices):
    ASSEMBLED_MOLECULE = 'AM', "assembled-molecule"
    UNLOCALIZED_SCAFFOLD = 'ULS', "unlocalized-scaffold"
    UNPLACED_SCAFFOLD = 'UPS', "unplaced-scaffold"
    ALT_SCAFFOLD = 'ALT', "alt-scaffold"
    FIX_PATCH = 'FP', "fix-patch"
    NOVEL_PATCH = 'NP', "novel-patch"
    # Not a real sequence role - the shared fake contig gene-level events (gene fusions) are
    # anchored on, so they have somewhere to sit with no coordinate. @see snpdb.gene_level_variants
    VG_GENE_LEVEL_FAKE_CONTIG = 'GL', "VariantGrid gene-level (not a sequence)"


class AssemblyMoleculeType(models.TextChoices):
    CHROMOSOME = 'C', "Chromosome"
    MITOCHONDRION = 'M', "Mitochondrion"


class AlleleOrigin(models.TextChoices):
    IMPORTED_TO_DATABASE = 'D', 'Imported as this build'
    IMPORTED_NORMALIZED = 'N', 'Imported (normalised)'
    LIFTOVER = 'L', 'Liftover'
    LIFTOVER_NORMALIZED = 'M', 'Liftover (normalised)'  # This probably shouldn't happen!

    @staticmethod
    def variant_origin(variant, allele, genome_build):
        try:
            origin = variant.variantallele_set.get(genome_build=genome_build, allele=allele).origin
        except ObjectDoesNotExist:
            # Variant w/o allele must have been imported directly
            if variant.modifiedimportedvariant_set.exists():
                origin = AlleleOrigin.IMPORTED_NORMALIZED
            else:
                origin = AlleleOrigin.IMPORTED_TO_DATABASE
        return origin


class AlleleConversionTool(models.TextChoices):
    SAME_CONTIG = "SC", "Identical Contig/Version"
    CLINGEN_ALLELE_REGISTRY = 'CA', "ClinGen Allele Registry"
    DBSNP = 'DB', "dbSNP API"
    NCBI_REMAP = 'NR', "NCBI Remap"  # This is obsolete as of November 2023
    PICARD = "PC", "Picard LiftoverVCF"
    CROSSMAP = "CM", "CrossMap"
    BCFTOOLS_LIFTOVER = "BL", "BCFtools/liftover"

    @classmethod
    def vcf_tuples_in_destination_build(cls, conversion_tool):
        IN_DEST_BUILD = {
            cls.SAME_CONTIG: True,
            cls.CLINGEN_ALLELE_REGISTRY: True,
            cls.DBSNP: True,
            cls.NCBI_REMAP: False,
            cls.BCFTOOLS_LIFTOVER: False,
        }
        return IN_DEST_BUILD[conversion_tool]


class AlleleOriginFilterDefault(models.TextChoices):
    SHOW_ALL = "A", "All"
    GERMLINE = "G", "Germline"
    SOMATIC = "S", "Somatic"

    @property
    def buckets(self) -> set[AlleleOriginBucket]:
        if self == AlleleOriginFilterDefault.SHOW_ALL:
            return set(AlleleOriginBucket)
        elif self == AlleleOriginFilterDefault.GERMLINE:
            return {AlleleOriginBucket.GERMLINE, AlleleOriginBucket.UNKNOWN}
        elif self == AlleleOriginFilterDefault.SOMATIC:
            return {AlleleOriginBucket.SOMATIC, AlleleOriginBucket.UNKNOWN}


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


class DataState(models.TextChoices):
    """ For files being loaded off disks """
    NON_EXISTENT = 'N', 'Non Existent'
    DELETED = 'D', 'Deleted'
    RUNNING = 'R', 'Running'
    SKIPPED = 'S', 'Skipped'
    ERROR = 'E', 'Error'
    COMPLETE = 'C', 'Complete'

    @staticmethod
    def should_create_new_record(data_state):
        return data_state not in [DataState.DELETED, DataState.SKIPPED]


class CohortGenotypeCollectionType(models.TextChoices):
    COMMON = "C", "Common"
    UNCOMMON = "U", "Uncommon"
