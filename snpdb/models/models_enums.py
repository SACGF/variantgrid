"""
Enums shared across snpdb and the apps built on it: import and processing lifecycles (ImportStatus,
ProcessingStatus, DataState), where uploads came from (ImportSource), allele linking
(AlleleOrigin, AlleleConversionTool, AlleleOriginFilterDefault), contig roles (SequenceRole - including the
gene-level fake contig), sample and file kinds, the grid's built-in and per-tag filter keys
(BuiltInFilters, TagFilter) and the user award scale. Codes are stored in the database: add values,
never renumber. Zygosity lives in patients/models_enums.py.

The classification enums that apps below classification also use live here too - AlleleOriginBucket, ShareLevel
(lab -> organisation -> logged-in users -> public, with which levels count toward discordance) and the
ClinicalSignificance / SomaticClinicalSignificance value helpers. classification.enums re-exports them.
"""
import typing
from enum import Enum
from functools import total_ordering
from typing import Optional, Union

from django.contrib.auth.models import User
from django.core.exceptions import ObjectDoesNotExist
from django.db import models
from django.db.models import TextChoices

from library.guardian_utils import all_users_group, public_group
from library.utils import ChoicesEnum, Constant


class AlleleOriginBucket(TextChoices):
    GERMLINE = "G", "Germline"
    SOMATIC = "S", "Somatic"
    UNKNOWN = "U", "Unknown"


class _ShareLevelData(typing.NamedTuple):
    index: int
    label: str


@total_ordering
class ShareLevel(ChoicesEnum):
    # note the following aren't enums (thus the _ignore_) type is defined here, contents is defined below ShareLevel definition
    _ignore_ = ['ALL_LEVELS', 'DISCORDANT_LEVEL_KEYS', '_DATA']
    ALL_LEVELS: list['ShareLevel']
    DISCORDANT_LEVEL_KEYS: list[str]
    _DATA: dict['ShareLevel', _ShareLevelData]

    # These strings have to be <= 16 characters for choice field
    CURRENT_USER = 'user'
    LAB = 'lab'
    ORGANISATION = 'organisation'
    ALL_USERS = 'logged_in_users'
    PUBLIC = 'public'

    @property
    def _data(self):
        return ShareLevel._DATA[self]

    @property
    def index(self):
        return self._data.index

    @property
    def label(self):
        return self._data.label

    def __lt__(self, other):
        return self.index < other.index

    @property
    def is_discordant_level(self) -> bool:
        return self.value in ShareLevel.DISCORDANT_LEVEL_KEYS

    def has_access(self, lab: 'Lab', user: User) -> bool:
        if user.is_superuser:
            return True
        match self:
            case _ if self.is_discordant_level:
                return True
            case ShareLevel.ORGANISATION:
                return user.groups.filter(name=lab.group_institution).exists()
            case ShareLevel.LAB:
                return user.groups.filter(name=lab.group).exists()
            case _:
                raise ValueError(f"Unhandled share level {self}")

    def group(self, lab: 'Lab', user: User = None):
        groups = {
            ShareLevel.CURRENT_USER: user,
            ShareLevel.LAB: lab.group,
            ShareLevel.ORGANISATION: lab.group_institution,
            ShareLevel.ALL_USERS: all_users_group(),
            ShareLevel.PUBLIC: public_group()
        }
        return groups.get(self)

    def context_label(self, vc: 'Classification') -> str:
        context_labels = {
            ShareLevel.CURRENT_USER: vc.user.username,
            ShareLevel.LAB: vc.lab.name,
            ShareLevel.ORGANISATION: vc.lab.organization.name,
        }
        return context_labels.get(self, self.label)

    # deprecated - just provided as a halfway measure when switching ShareLevel from class to Enum
    @property
    def key(self) -> str:
        return self.value

    @property
    def icon(self) -> str:
        return f'icons/share_level/{self.value}.png'

    @classmethod
    def _missing_(cls, value):
        # Backwards compatibility: old DB rows used 'institution' before the rename to 'organisation'
        if value == 'institution':
            return cls.ORGANISATION
        return None

    def __str__(self):
        return self.value

    @staticmethod
    def same_and_higher(level: 'ShareLevel') -> list['ShareLevel']:
        return [iter_level for iter_level in ShareLevel.ALL_LEVELS if iter_level.index >= level.index]

    @staticmethod
    def from_key(source: Optional[Union['ShareLevel', str, int]]) -> Optional['ShareLevel']:
        if source is None:
            return None

        if isinstance(source, ShareLevel):
            return source

        index = -1
        if isinstance(source, int):
            index = source
        elif str(source).isnumeric():
            index = int(source)

        # Backwards compatibility: old DB rows and old API clients used 'institution'
        if str(source) == 'institution':
            return ShareLevel.ORGANISATION

        for sl in ShareLevel.ALL_LEVELS:
            if sl.value == str(source) or sl.index == index:
                return sl
        raise ValueError('Unknown share level ' + str(source))


ShareLevel._DATA = {
        ShareLevel.CURRENT_USER: _ShareLevelData(index=0, label='Current User'),
        ShareLevel.LAB: _ShareLevelData(index=1, label='Lab'),
        ShareLevel.ORGANISATION: _ShareLevelData(index=2, label='Organisation'),
        ShareLevel.ALL_USERS: _ShareLevelData(index=3, label='App Users'),
        ShareLevel.PUBLIC: _ShareLevelData(index=4, label='3rd Party Databases')
    }
ShareLevel.ALL_LEVELS = [ShareLevel.CURRENT_USER, ShareLevel.LAB, ShareLevel.ORGANISATION, ShareLevel.ALL_USERS, ShareLevel.PUBLIC]
ShareLevel.DISCORDANT_LEVEL_KEYS = {ShareLevel.ALL_USERS.value, ShareLevel.PUBLIC.value}


class ClinicalSignificance:
    OTHER = '0'
    BENIGN = '1'
    LIKELY_BENIGN = '2'
    VUS = '3'
    LIKELY_PATHOGENIC = '4'
    PATHOGENIC = '5'

    CHOICES = [
        (OTHER, 'Other'),
        (BENIGN, 'Benign'),
        (LIKELY_BENIGN, 'Likely Benign'),
        (VUS, 'VUS'),
        (LIKELY_PATHOGENIC, 'Likely Pathogenic'),
        (PATHOGENIC, 'Pathogenic'),
    ]

    LABELS = dict(CHOICES + [(None, "No Data")])

    SHORT_CHOICES = [
        (OTHER, 'O'),
        (BENIGN, 'B'),
        (LIKELY_BENIGN, 'LB'),
        (VUS, 'VUS'),
        (LIKELY_PATHOGENIC, 'LP'),
        (PATHOGENIC, 'P'),
    ]
    SHORT_LABELS = dict(SHORT_CHOICES + [(None, "U")])

    CHART_COLOURS = {
        'B': "#44d", 'LB': "#88d", 'VUS': "#666", 'LP': "#d88", 'P': "#d44", 'O': "#aaa",
    }
    """ Saturated equivalents of the pale .c-pill.cs- fills in global.scss, keyed by SHORT_LABELS """

    DEFAULT_CHART_COLOUR = "#aaa"

    @staticmethod
    def css_class(clinical_significance: Optional[str]) -> str:
        """ Class for the .c-pill.cs-* rules in global.scss, from a CHOICES value """
        if clinical_significance and (short_label := ClinicalSignificance.SHORT_LABELS.get(clinical_significance)):
            return f"cs-{short_label.lower()}"
        return "cs-none"

    @staticmethod
    def chart_colour(clinical_significance: Optional[str]) -> str:
        """ Saturated equivalent of the .cs- fill, for chart marks and swatches """
        short_label = ClinicalSignificance.SHORT_LABELS.get(clinical_significance)
        return ClinicalSignificance.CHART_COLOURS.get(short_label, ClinicalSignificance.DEFAULT_CHART_COLOUR)

    @staticmethod
    def is_significant_change(old_classification: str, new_classification: str) -> bool:
        was_vus_change = old_classification and new_classification and old_classification == 'VUS' and new_classification.startswith('VUS')
        return old_classification != new_classification and not was_vus_change

    @staticmethod
    def distance(old_classification: int|str, new_classification: int|str) -> int | None:
        """ Returns distance (0-4) between classifications (returns None if any are OTHER) """
        d = None
        try:
            oc_val = int(old_classification)
            nc_val = int(new_classification)
            if oc_val and nc_val:
                return oc_val - nc_val
        except ValueError:
            pass
        return d


class SomaticClinicalSignificance:
    """ Values of the 'somatic:clinical_significance' EvidenceKey (AMP tiers) that code needs to reference
        directly - the authoritative option list lives on the EvidenceKey record """
    TIER_1 = "tier_1"
    TIER_1_OR_2 = "tier_1_or_2"
    TIER_2 = "tier_2"
    TIER_3 = "tier_3"
    TIER_4 = "tier_4"

    CHOICES = [
        (TIER_1, "Tier I"),
        (TIER_1_OR_2, "Tier I/II"),
        (TIER_2, "Tier II"),
        (TIER_3, "Tier III"),
        (TIER_4, "Tier IV"),
    ]
    LABELS = dict(CHOICES)

    SHORT_CHOICES = [
        (TIER_1, "I"),
        (TIER_1_OR_2, "I/II"),
        (TIER_2, "II"),
        (TIER_3, "III"),
        (TIER_4, "IV"),
    ]
    SHORT_LABELS = dict(SHORT_CHOICES)

    TIER_1_AND_2_VALUES = [TIER_1, TIER_1_OR_2, TIER_2]
    """ Tiers of strong/potential clinical significance - the somatic analogue of germline LP/P """

    @staticmethod
    def css_class(somatic_clinical_significance: Optional[str]) -> str:
        """ Class for the .cs-* rules in global.scss """
        if somatic_clinical_significance in SomaticClinicalSignificance.LABELS:
            return f"cs-{somatic_clinical_significance}"
        return "cs-none"


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
    BADGE = "B", "Badge"  # Permanent once earned, tiered bronze/silver/gold
    KUDOS = "K", "Kudos"  # Hand-given by an admin


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

    @property
    def icon_css(self) -> str:
        """ Font Awesome classes for the status - the one place the icons are chosen, so the data grids,
            the VCF / sample pages and the sequencing run page all agree (#2175) """
        return {
            ImportStatus.CREATED: "fas fa-spinner fa-spin text-warning",
            ImportStatus.IMPORTING: "fas fa-spinner fa-spin text-warning",
            ImportStatus.REQUIRES_USER_INPUT: "fas fa-exclamation-triangle text-warning",
            ImportStatus.ERROR: "fas fa-times-circle text-danger",
            ImportStatus.SUCCESS: "fas fa-check-circle text-success",
            ImportStatus.MARKED_FOR_DELETION: "fas fa-trash-alt text-muted",
            ImportStatus.DELETING: "fas fa-trash-alt text-muted",
        }[self]


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


class VariantOrigin(models.TextChoices):
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

    ALIGNMENT_TYPES = Constant([e[0] for e in (BAM, CRAM)])  # Reads, ie what IGV is opened with


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


class DuoRelationship(models.TextChoices):
    """ How a Duo's second member is related to the proband - a parent, or a sibling (#1861) """
    MOTHER = 'M', 'Mother'
    FATHER = 'F', 'Father'
    SIBLING = 'S', 'Sibling'
