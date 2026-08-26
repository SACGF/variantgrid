import typing
from enum import Enum
from functools import total_ordering
from typing import Optional, Union

from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import Q, TextChoices

from library.guardian_utils import all_users_group, public_group
from library.utils import ChoicesEnum

CRITERIA_NOT_MET = 'NM'
CRITERIA_NOT_APPLICABLE = 'NA'
CRITERIA_NEUTRAL = 'N'


class AlleleOriginBucket(TextChoices):
    GERMLINE = "G", "Germline"
    SOMATIC = "S", "Somatic"
    UNKNOWN = "U", "Unknown"

    @staticmethod
    def bucket_for_allele_origin(allele_origin: Optional[str]) -> 'AlleleOriginBucket':
        # logic is duplicated in JavaSCript in vc_form.js updateTitle()
        if not allele_origin:
            return AlleleOriginBucket(settings.ALLELE_ORIGIN_NOT_PROVIDED_BUCKET)

        # TODO, it would be good to put the bucket directly into allele origin
        from classification.models import EvidenceKeyMap
        if bucket := EvidenceKeyMap.cached_key(SpecialEKeys.ALLELE_ORIGIN).matched_options(allele_origin)[0].get("bucket"):
            return AlleleOriginBucket(bucket)

        allele_origin = allele_origin.lower()
        if "somatic" in allele_origin:
            return AlleleOriginBucket.SOMATIC
        elif "germline" in allele_origin:
            return AlleleOriginBucket.GERMLINE
        else:
            return AlleleOriginBucket.UNKNOWN


class ReclassificationEventType(TextChoices):
    """ What one row of a classification's significance timeline records @see ReclassificationEvent """
    INITIAL = "I", "Initial Classification"
    RECLASSIFICATION = "R", "Reclassification"
    REEVALUATION = "E", "Re-evaluation"

    @staticmethod
    def review_types() -> list[str]:
        """ A curator looked at the record - a reclassification is a review that moved the call """
        return [ReclassificationEventType.RECLASSIFICATION, ReclassificationEventType.REEVALUATION]


class LabExternalFilter(TextChoices):
    """ Classification grid filter for labs whose records came from elsewhere, e.g. synced down from Shariant """
    ALL = "A", "All"
    INTERNAL = "I", "Internal"
    EXTERNAL = "E", "External"

    def filter_q(self, lab_path: str = "lab") -> Optional[Q]:
        """ Q to apply against the given path to Lab, None if all labs should be included """
        if self in _LAB_EXTERNAL_VALUES:
            return Q(**{f"{lab_path}__external": _LAB_EXTERNAL_VALUES[self]})


# ALL is absent so that it doesn't filter on Lab.external at all
_LAB_EXTERNAL_VALUES = {
    LabExternalFilter.INTERNAL: False,
    LabExternalFilter.EXTERNAL: True
}


class SpecialEKeys:
    AUTOPOPULATE = 'autopopulate'
    VARIANT_COORDINATE = 'variant_coordinate'
    G_HGVS = 'g_hgvs'
    C_HGVS = 'c_hgvs'
    P_HGVS = 'p_hgvs'
    CONDITION = 'condition'
    CLINICAL_SIGNIFICANCE = 'clinical_significance'
    SOMATIC_CLINICAL_SIGNIFICANCE = 'somatic:clinical_significance'
    CURATED_BY = 'curated_by'
    CURATION_DATE = 'curation_date'
    CURATION_VERIFIED_DATE = 'curation_verified_date'
    CURATION_VERIFIED_BY = 'curation_verified_by'
    SAMPLE_DATE = 'sample_date'
    VARIANT_TAGS = 'variant_tags'

    # POPULATED
    # Note: Some fields not here are populated - those with variantgrid_column
    # and the pops - ie "pop_AFR" "pop_NFE" etc.
    AGE = "age"
    AGE_UNITS = "age_units"  # deleted now, but declared fo migrations
    ALLELE_DEPTH = 'allele_depth'
    ALLELE_FREQUENCY = 'allele_frequency'
    ALLELE_ORIGIN = "allele_origin"
    ALOFT = "aloft"
    ANCESTRY = "ancestry"
    CAPTURE_METHOD = "capture_method"
    CLINGEN_ALLELE_ID = "clingen_allele_id"
    CURATION_SYSTEM = "curation_system"
    DB_SNP = "db_rs_id"
    ENSEMBL_GENE_ID = "ensembl_gene_id"
    ENSEMBL_TRANSCRIPT_ID = "ensembl_transcript_id"
    ENTREZ_GENE_ID = "entrez_gene_id"
    FAMILY_ID = "family_id"
    GENE_SYMBOL = 'gene_symbol'
    GENOME_BUILD = "genome_build"
    GENOTYPE_QUALITY = 'genotype_quality'
    GNOMAD_AF = "gnomad_af"
    GNOMAD_OE_LOF = "gnomad_oe_lof"
    INTERNAL_SAMPLES_20X_COVERAGE = "internal_samples_100_percent_20x_gene_coverage"
    LITERATURE = "literature"
    LRG_ID = "lrg_id"
    NUCLEIC_ACID_SOURCE = "nucleic_acid_source"
    PATIENT_ID = "patient_id"
    PHRED_LIKELIHOOD = 'phred_likelihood'
    PHASTCONS = "phastcons"
    PHYLOP = "phylop"
    PUBMED_GENE_SEARCH_COUNT = "pubmed_gene_search_count"
    PUBMED_SEARCH_TERMS = "pubmed_search_terms"
    READ_DEPTH = 'read_depth'
    REFERENCE_DEPTH = 'reference_depth'
    REFSEQ_TRANSCRIPT_ID = "refseq_transcript_id"
    SAMPLE_ID = "sample_id"
    SAMPLE_TYPE = "sample_type"
    SEARCH_TERMS = "search_terms"
    SEQUENCING_PLATFORM = "sequencing_platform"
    SEX = "sex"
    SPECIMEN_ID = "specimen_id"
    SPLICEAI = "spliceai"
    VCF_FILTER = "vcf_filter"
    UNIPROT_ID = "uniprot_id"
    ZYGOSITY = "zygosity"

    # used for Clinvar export
    AFFECTED_STATUS = 'affected_status'
    MODE_OF_INHERITANCE = 'mode_of_inheritance'
    CURATION_CONTEXT = 'curation_context'
    ASSERTION_METHOD = 'assertion_method'
    PATIENT_PHENOTYPE = 'patient_phenotype'
    ORIGIN = 'origin'
    OWNER = 'owner'
    SOURCE_ID = 'source_id'
    SOURCE_URL = 'source_url'

    INTERPRETATION_SUMMARY = 'interpretation_summary'

    ANNOTATION_CONSORTIUM_KEYS = {
        'E': ENSEMBL_TRANSCRIPT_ID,
        'R': REFSEQ_TRANSCRIPT_ID,
    }
    VARIANT_LINKING_HGVS_KEYS = [C_HGVS, G_HGVS]
    VARIANT_LINKING_KEYS = VARIANT_LINKING_HGVS_KEYS + [VARIANT_COORDINATE]

    AMP_LEVELS_TO_LEVEL = {
        'amp:level_a': "A",
        'amp:level_b': "B",
        'amp:level_c': "C",
        'amp:level_d': "D"
    }


class ValidationCode:
    # Basic Errors
    INVALID_DATE = "invalid_date"
    INVALID_FLOAT = "invalid_float"
    INVALID_INTEGER = "invalid_integer"
    INVALID_PERCENT = "invalid_percent"
    INVALID_UNIT = "invalid_unit"
    INVALID_VALUE = "invalid_value"
    MANDATORY = "mandatory"
    USER_NOT_FOUND = "user_not_found"
    PARSE_ERROR = "cant_parse"
    REQUIRES_NOTE = "requires_note"
    TOO_MANY_VALUES = "too_many_values"
    UNKNOWN_KEY = "unknown_key"

    # Advanced Errors (might move them out into the files that raise them)
    NOT_LAB = "not_lab"
    NORMALIZING = 'normalizing'
    MATCHING_ERROR = "matching_error"
    INCONSISTENT_VARIANT = "inconsistent_variant"  # >1 VARIANT_LINKING_KEYS resolved to different coordinate
    UNKNOWN_TRANSCRIPT = "unknown_transcript"
    INTERNAL_ERROR = "internal_error"

    # Missing clinical significance and variant classification
    MISSING_CLINICAL_SIGNIFICANCE = 'missing_significance'
    INVALID_FIELD_FOR_SOMATIC = 'invalid_for_somatic'
    INVALID_FIELD_FOR_GERMLINE = 'invalid_for_germline'

    SOMATIC_MISMATCHED_LEVEL = 'somatic_mismatched_level'


class EvidenceCategory:
    VARIANT = 'V'

    GENE = 'H'

    HEADER_PATIENT = 'HP'
    HEADER_TEST = 'HT'

    POPULATION_DATA = 'P'
    COMPUTATIONAL_AND_PREDICTIVE_DATA = 'CP'
    FUNCTIONAL_DATA = 'F'
    SEGREGATION_DATA = 'S'
    DE_NOVO_DATA = 'D'
    ALLELIC_DATA = 'A'
    OTHER_DATABASE = 'DB'
    OTHER_DATA = 'O'
    SOMATIC_CLINICAL_SIGNIFICANCE_EVIDENCE = "SC"
    INTERPRETATION = 'HI'
    SIGN_OFF = 'SO'
    LITERATURE = 'L'
    # Summary data covers things like literature
    # Things that are typically cross evidence concerns that have been bundled up
    # in one spot

    CHOICES = (
        (VARIANT, 'Variant'),
        (GENE, 'Gene'),
        (HEADER_PATIENT, 'Patient'),
        (HEADER_TEST, 'Test'),

        (POPULATION_DATA, 'Population Data'),
        (COMPUTATIONAL_AND_PREDICTIVE_DATA, 'Computational and Predictive Data'),
        (FUNCTIONAL_DATA, 'Functional Data'),
        (SEGREGATION_DATA, 'Segregation Data'),
        (DE_NOVO_DATA, 'De novo Data'),
        (ALLELIC_DATA, 'Allelic Data'),
        (OTHER_DATABASE, 'Other Database'),
        (OTHER_DATA, 'Other Data'),
        (LITERATURE, 'Literature'),
        (SOMATIC_CLINICAL_SIGNIFICANCE_EVIDENCE, "Somatic Clinical Significance Evidence"),
        (INTERPRETATION, 'Interpretation'),
        (SIGN_OFF, 'Sign Off')
    )


class EvidenceKeyValueType:
    AGE = 'A'
    TEXT_AREA = 'T'
    FREE_ENTRY = 'F'
    SELECT = 'S'
    MULTISELECT = 'M'
    BOOLEAN = 'B'
    DATE = 'D'
    CRITERIA = 'C'
    USER = 'U'
    UNIT = 'N'
    INTEGER = 'I'
    FLOAT = 'L'
    CHOICES = (
        (AGE, 'age'),
        (TEXT_AREA, 'free text (multi-line)'),
        (FREE_ENTRY, 'free text'),
        (SELECT, 'select'),
        (MULTISELECT, 'multi-select'),
        (BOOLEAN, 'boolean'),
        (DATE, 'date'),
        (CRITERIA, 'criteria'),
        (USER, 'user'),
        (UNIT, 'unit (0 to 1)'),
        (INTEGER, 'integer'),
        (FLOAT, 'float')
    )


class _ShareLevelData(typing.NamedTuple):
    index: int
    label: str


@total_ordering
class ShareLevel(ChoicesEnum):
    _ignore_ = ['ALL_LEVELS', 'DISCORDANT_LEVEL_KEYS']
    ALL_LEVELS: list['ShareLevel'] = []
    DISCORDANT_LEVEL_KEYS: list[str] = []

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


class ForceUpdate(str, Enum):
    # NONE = None - when we've got default behaviour
    SOURCE = 'source'  # supports changing source_id and curation_date
    # ALL (in future might support this for updating of genome_build or other immutable fields)


class SubmissionSource(str, Enum):
    FORM = 'form'
    CONSENSUS = 'consensus'  # now is copy from latest
    API = 'api'
    VARIANT_GRID = 'variantgrid'

    @property
    def is_valid_user_source(self) -> bool:
        return self in (SubmissionSource.FORM, SubmissionSource.CONSENSUS, SubmissionSource.API)

    def can_edit(self, immutability_level: 'SubmissionSource'):
        return self.level >= immutability_level.level

    @property
    def level(self) -> int:
        # for legacy reasons, immutable can = True
        if self in (SubmissionSource.FORM, SubmissionSource.CONSENSUS):
            return 1
        if self == SubmissionSource.API:
            return 2
        if self == SubmissionSource.VARIANT_GRID:
            return 3
        return 0


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
        """ Class for the .c-pill.scs-* rules in global.scss """
        if somatic_clinical_significance in SomaticClinicalSignificance.LABELS:
            return f"scs-{somatic_clinical_significance}"
        return "scs-none"


class CriteriaEvaluation:
    NOT_MET = CRITERIA_NOT_MET
    NOT_APPLICABLE = CRITERIA_NOT_APPLICABLE

    # UNSPECIFIED STRENGTH HANDLING

    BENIGN_UNSPECIFIED = 'BX'
    BENIGN_STANDALONE = 'BA'
    BENIGN_STRONG = 'BS'
    BENIGN_MODERATE = 'BM'  # Not a standard ACMG Strength
    BENIGN_SUPPORTING = 'BP'
    NEUTRAL = 'N'
    PATHOGENIC_SUPPORTING = 'PP'
    PATHOGENIC_MODERATE = 'PM'
    PATHOGENIC_STRONG = 'PS'
    PATHOGENIC_VERY_STRONG = 'PVS'
    PATHOGENIC_UNSPECIFIED = 'PX'

    POINTS = {
        BENIGN_STANDALONE: -8,
        BENIGN_STRONG: -4,
        BENIGN_MODERATE: -2,
        BENIGN_SUPPORTING: -1,
        NOT_MET: 0,
        NOT_APPLICABLE: 0,
        NEUTRAL: 0,
        PATHOGENIC_SUPPORTING: 1,
        PATHOGENIC_MODERATE: 2,
        PATHOGENIC_STRONG: 4,
        PATHOGENIC_VERY_STRONG: 8
    }

    CHOICES = (
        (BENIGN_STANDALONE, 'Benign Standalone'),
        (BENIGN_STRONG, 'Benign Strong'),
        (BENIGN_MODERATE, 'Benign Moderate'),
        (BENIGN_SUPPORTING, 'Benign Supporting'),
        (BENIGN_UNSPECIFIED, 'Benign Unspecified Strength'),
        (NEUTRAL, 'Neutral'),
        (PATHOGENIC_SUPPORTING, 'Pathogenic Supporting'),
        (PATHOGENIC_MODERATE, 'Pathogenic Moderate'),
        (PATHOGENIC_STRONG, 'Pathogenic Strong'),
        (PATHOGENIC_VERY_STRONG, "Pathogenic Very Strong"),
        (PATHOGENIC_UNSPECIFIED, 'Pathogenic Unspecified Strength'),
    )

    # Neutral, Not Met, Not Applicable don't count
    ALL_STRENGTHS = [BENIGN_STANDALONE, BENIGN_STRONG, BENIGN_SUPPORTING, BENIGN_UNSPECIFIED,
                     PATHOGENIC_SUPPORTING, PATHOGENIC_MODERATE, PATHOGENIC_STRONG, PATHOGENIC_VERY_STRONG, PATHOGENIC_UNSPECIFIED]

    ####
    # used for drop-downs

    BENIGN_OPTIONS = [
        {'key': CRITERIA_NOT_MET, 'label': 'Not Met', 'index': 1},
        {'key': CRITERIA_NOT_APPLICABLE, 'label': 'Not Applicable', 'index': 9},
        {'key': NEUTRAL, 'label': 'Neutral', 'index': 10},
        {'key': BENIGN_STANDALONE, 'label': 'Benign Standalone', 'index': 2},
        {'key': BENIGN_STRONG, 'label': 'Benign Strong', 'index': 3},
        {'key': BENIGN_MODERATE, 'label': 'Benign Moderate', 'index': 9},  # NOT a STANDARD ACMG STRENGTH
        {'key': BENIGN_SUPPORTING, 'label': 'Benign Supporting', 'index': 4},
        {'key': BENIGN_UNSPECIFIED, 'label': 'Benign Unspecified Strength', 'index': 11},
    ]

    NEUTRAL_OPTIONS = [
        {'key': CRITERIA_NOT_MET, 'label': 'Not Met', 'index': 1},
        {'key': CRITERIA_NOT_APPLICABLE, 'label': 'Not Applicable', 'index': 9},
        {'key': NEUTRAL, 'label': 'Neutral', 'index': 10}
    ]

    PATHOGENIC_OPTIONS = [
        {'key': CRITERIA_NOT_MET, 'label': 'Not Met', 'index': 1},
        {'key': CRITERIA_NOT_APPLICABLE, 'label': 'Not Applicable', 'index': 9},
        {'key': NEUTRAL, 'label': 'Neutral', 'index': 10},
        {'key': PATHOGENIC_SUPPORTING, 'label': 'Pathogenic Supporting', 'index': 5},
        {'key': PATHOGENIC_MODERATE, 'label': 'Pathogenic Moderate', 'index': 6},
        {'key': PATHOGENIC_STRONG, 'label': 'Pathogenic Strong', 'index': 7},
        {'key': PATHOGENIC_VERY_STRONG, 'label': 'Pathogenic Very Strong', 'index': 8},
        {'key': PATHOGENIC_UNSPECIFIED, 'label': 'Pathogenic Unspecified Strength', 'index': 11},
    ]

    @staticmethod
    def is_met(criteria) -> bool:
        return bool(criteria and criteria not in {CriteriaEvaluation.NOT_MET, CriteriaEvaluation.NOT_APPLICABLE, CriteriaEvaluation.NEUTRAL})


class WithdrawReason(TextChoices):
    SHARED_BY_MISTAKE = 'SHARED_BY_MISTAKE', 'This record was shared by mistake'
    DUPLICATE = 'DUPLICATE', 'This record has an exact duplicate'
    SENSITIVE_DATA = 'SENSITIVE_DATA', 'This record contains sensitive data'
    OTHER = 'OTHER', 'Other'
