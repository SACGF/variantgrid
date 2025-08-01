from django.db import models

from library.utils import Constant


class HumanProteinAtlasAbundance(models.TextChoices):
    NOT_DETECTED = 'N', 'Not detected'
    LOW = 'L', 'Low'
    MEDIUM = 'M', 'Medium'
    HIGH = 'H', 'High'

    @staticmethod
    def get_abundance_or_greater_levels(abundance):
        higher_levels = []
        found = False
        for k, _ in HumanProteinAtlasAbundance.choices:
            if k == abundance:
                found = True
            if found:
                higher_levels.append(k)
        return higher_levels


class DetectedHumanProteinAtlasAbundance(models.TextChoices):
    LOW = 'L', 'Low'
    MEDIUM = 'M', 'Medium'
    HIGH = 'H', 'High'


class AnnotationStatus(models.TextChoices):
    CREATED = 'C', "Created"
    DELETING = 'x', "Deleting"
    DUMP_STARTED = 'd', "Dump Started"
    DUMP_COMPLETED = 'D', "Dump Completed"
    ANNOTATION_STARTED = 'a', "Annotation Started"
    ANNOTATION_COMPLETED = 'A', "Annotation Completed"
    UPLOAD_STARTED = 'U', "Upload Started"
    FINISHED = 'F', "Finished"
    ERROR = 'E', "Error"

    @classmethod
    def get_summary_state(cls, annotation_status):
        SUMMARY_STATES = {cls.CREATED: "Queued",
                          cls.FINISHED: "Finished",
                          cls.ERROR: "Error"}
        return SUMMARY_STATES.get(annotation_status, "Running")


class HPOSynonymScope(models.TextChoices):
    BROAD = 'B', 'Broad'
    EXACT = 'E', 'Exact'
    NARROW = 'N', 'Narrow'
    RELATED = 'R', 'Related'


class ClinGenClassification(models.TextChoices):
    DEFINITIVE = 'D', 'Definitive'
    STRONG = 'S', 'Strong'
    MODERATE = 'M', 'Moderate'
    LIMITED = 'L', 'Limited'
    NO_KNOWN_DISEASE_RELATIONSHIP = 'N', 'No Known Disease Relationship'
    REFUTED = 'R', 'Refuted'
    DISPUTED = 'P', 'Disputed'


class VariantAnnotationPipelineType(models.TextChoices):
    """ We have standard long and short  """
    STANDARD = "S", "Standard Short Variant"
    STRUCTURAL_VARIANT = "C", "Structural Variant"


class ColumnAnnotationCategory(models.TextChoices):
    """ Based on categories from:
        https://asia.ensembl.org/info/docs/tools/vep/script/vep_plugins.html """

    CONSERVATION = 'C', "Conservation"
    EXTERNAL_ID = 'E', "External ID"
    FREQUENCY_DATA = 'F', "Frequency Data"
    FUNCTIONAL_EFFECT = 'f', "Functional Effect"
    GENE_ANNOTATIONS = 'G', 'Gene Annotations'
    HGVS = 'H', "HGVS"
    LITERATURE = 'L', 'Literature'
    NEARBY_FEATURES = 'N', "Nearby Features"
    PATHOGENICITY_PREDICTIONS = 'P', "Pathogenicity Predictions"
    PHENOTYPE = 'Y', "Phenotype"
    PROTEIN_DOMAINS = 'D', "Protein Domains"
    SEQUENCE = 'Q', "Sequence"
    SPLICING_PREDICTIONS = 'S', "Splicing Predictions"
    VARIANT_DATA = 'V', "Variant Data"


class VEPPlugin(models.TextChoices):
    DBNSFP = 'd', 'dbNSFP'
    DBSCSNV = 'v', 'dbscSNV'
    GRANTHAM = 'g', 'Grantham'
    LOFTOOL = 'l', 'LoFtool'
    MASTERMIND = 'n', 'Mastermind'
    MAVEDB = 'V', "MaveDb"
    MAXENTSCAN = 'm', 'MaxEntScan'
    NMD = "N", 'NMD'
    SPLICEAI = 'a', 'SpliceAI'
    SPLICEREGION = 's', 'SpliceRegion'


class VEPCustom(models.TextChoices):
    GNOMAD_2 = 'g', 'gnomAD2'
    GNOMAD_3 = 'n', 'gnomAD3'
    GNOMAD_4 = 'o', 'gnomAD4'
    # We split GNOMAD_SV into 2, as we need to call custom twice
    GNOMAD_SV = 'S', 'gnomAD_SV'
    GNOMAD_SV_NAME = 'N', 'gnomAD_SV_name'
    PHASTCONS_100_WAY = '1', 'phastCons100way'
    PHASTCONS_30_WAY = '2', 'phastCons30way'
    PHASTCONS_46_WAY = '3', 'phastCons46way'
    PHYLOP_100_WAY = '4', 'phyloP100way'
    PHYLOP_30_WAY = '5', 'phyloP30way'
    PHYLOP_46_WAY = '6', 'phyloP46way'
    REPEAT_MASKER = 'r', 'RepeatMasker'
    TOPMED = 't', 'TopMed'
    UK10K = 'u', 'UK10k'
    COSMIC = 'c', 'COSMIC'


class VEPSkippedReason(models.TextChoices):
    UNKNOWN_CONTIG = 'c', "Unknown Contig"
    INCOMPLETE = 'i', "Incomplete"
    UNKNOWN = 'u', "Unknown"
    TOO_LONG = 'l', "Too Long"


class ClinVarReviewStatus(models.TextChoices):
    NO_ASSERTION_PROVIDED = "N", "No assertion provided"
    NO_ASSERTION_CRITERIA_PROVIDED = "C", "No assertion criteria provided"
    NO_INTERPRETATION_FOR_THE_SINGLE_VARIANT = "I", "No interpretation for the single variant"
    NO_CLASSIFICATION_PROVIDED = 'n', "No Classification Provided"
    NO_CLASSIFICATIONS_FROM_UNFLAGGED_RECORDS = 'f', "No classifications from unflaggedrecords"
    CRITERIA_PROVIDED_CONFLICTING_INTERPRETATIONS = "F", "Criteria provided - conflicting interpretations"
    CRITERIA_PROVIDED_SINGLE_SUBMITTER = "S", "Criteria provided - single submitter"
    CRITERIA_PROVIDED_MULTIPLE_SUBMITTERS_NO_CONFLICTS = "M", "Criteria provided - multiple submitters w/no conflicts"
    CRITERIA_PROVIDED_MULTIPLE_SUBMITTERS = "m", "Criteria provided - multiple submitters"  # somatic only
    REVIEWED_BY_EXPERT_PANEL = "E", "Reviewed by expert panel"
    PRACTICE_GUIDELINE = "P", "Practice guideline"

    STARS = Constant({
        NO_CLASSIFICATION_PROVIDED[0]: 0,
        NO_ASSERTION_PROVIDED[0]: 0,
        NO_ASSERTION_CRITERIA_PROVIDED[0]: 0,
        NO_INTERPRETATION_FOR_THE_SINGLE_VARIANT[0]: 0,
        NO_CLASSIFICATIONS_FROM_UNFLAGGED_RECORDS[0]: 0,
        CRITERIA_PROVIDED_CONFLICTING_INTERPRETATIONS[0]: 1,
        CRITERIA_PROVIDED_SINGLE_SUBMITTER[0]: 1,
        CRITERIA_PROVIDED_MULTIPLE_SUBMITTERS_NO_CONFLICTS[0]: 2,
        CRITERIA_PROVIDED_MULTIPLE_SUBMITTERS[0]: 2,
        REVIEWED_BY_EXPERT_PANEL[0]: 3,
        PRACTICE_GUIDELINE[0]: 4,
    })

    VCF_MAPPINGS = Constant({
        'no_assertion_provided': NO_ASSERTION_PROVIDED[0],
        'no_assertion_criteria_provided': NO_ASSERTION_CRITERIA_PROVIDED[0],
        'no_interpretation_for_the_single_variant': NO_INTERPRETATION_FOR_THE_SINGLE_VARIANT[0],  # Old value
        'no_classification_for_the_single_variant': NO_INTERPRETATION_FOR_THE_SINGLE_VARIANT[0],  # New one
        'no_classification_provided': NO_CLASSIFICATION_PROVIDED[0],  # New one
        'no_classifications_from_unflagged_records': NO_CLASSIFICATIONS_FROM_UNFLAGGED_RECORDS[0],  # new for oncogenic
        'criteria_provided,_conflicting_interpretations': CRITERIA_PROVIDED_CONFLICTING_INTERPRETATIONS[0],
        'criteria_provided,_conflicting_classifications': CRITERIA_PROVIDED_CONFLICTING_INTERPRETATIONS[0],
        'criteria_provided,_single_submitter': CRITERIA_PROVIDED_SINGLE_SUBMITTER[0],
        'criteria_provided,_multiple_submitters': CRITERIA_PROVIDED_MULTIPLE_SUBMITTERS[0],  # new for somatic
        'criteria_provided,_multiple_submitters,_no_conflicts': CRITERIA_PROVIDED_MULTIPLE_SUBMITTERS_NO_CONFLICTS[0],
        'reviewed_by_expert_panel': REVIEWED_BY_EXPERT_PANEL[0],
        'practice_guideline': PRACTICE_GUIDELINE[0],
    })

    def stars(self):
        return ClinVarReviewStatus.STARS[self.value]

    @staticmethod
    def statuses_gte_stars(min_stars: int) -> list['ClinVarReviewStatus']:
        statuses = []
        for rs, stars in ClinVarReviewStatus.STARS.items():
            if stars >= min_stars:
                statuses.append(rs)
        return statuses


class ManualVariantEntryType(models.TextChoices):
    DBSNP = "d", "dbSNP"
    HGVS = "h", "HGVS"
    VARIANT = "v", "Variant"
    UNKNOWN = "u", "Unknown"


class EssentialGeneCRISPR(models.TextChoices):
    ESSENTIAL = "E", "Essential"
    NON_ESSENTIAL_PHENOTYPE_CHANGING = "N", "Non-essential phenotype-changing"


class EssentialGeneCRISPR2(models.TextChoices):
    ESSENTIAL = "E", "Essential"
    NON_ESSENTIAL_PHENOTYPE_CHANGING = "N", "Non-essential phenotype-changing"
    CONTEXT_SPECIFIC_ESSENTIAL = "S", "Context-Specific essential"


class EssentialGeneGeneTrap(models.TextChoices):
    ESSENTIAL = "E", "Essential"
    NON_ESSENTIAL_PHENOTYPE_CHANGING = "N", "Non-essential phenotype-changing"
    HAP1_SPECIFIC_ESSENTIAL = "H", "HAP1-Specific essential"
    KBM7_SPECIFIC_ESSENTIAL = "K", "KBM7-Specific essential"
