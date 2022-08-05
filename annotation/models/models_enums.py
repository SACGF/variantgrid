from typing import List

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


class CitationSource(models.TextChoices):
    PUBMED = 'P', 'PubMed'
    NCBI_BOOKSHELF = 'N', 'NCBIBookShelf'
    PUBMED_CENTRAL = 'C', 'PubMedCentral'

    CODES = Constant({'PubMed': PUBMED[0],
                      'PMID': PUBMED[0],
                      'NCBIBookShelf': NCBI_BOOKSHELF[0],
                      'PubMedCentral': PUBMED_CENTRAL[0]})


class ClinGenClassification(models.TextChoices):
    DEFINITIVE = 'D', 'Definitive'
    STRONG = 'S', 'Strong'
    MODERATE = 'M', 'Moderate'
    LIMITED = 'L', 'Limited'
    NO_KNOWN_DISEASE_RELATIONSHIP = 'N', 'No Known Disease Relationship'
    REFUTED = 'R', 'Refuted'
    DISPUTED = 'P', 'Disputed'


class VariantClass(models.TextChoices):
    """ https://asia.ensembl.org/info/genome/variation/prediction/classification.html#classes """

    SNV = 'SN', "SNV"
    GENETIC_MARKER = 'GM', "genetic_marker"
    SUBSTITUTION = 'SU', "substitution"
    TANDEM_REPEAT = 'TR', "tandem_repeat"
    ALU_INSERTION = 'AI', "Alu_insertion"
    COMPLEX_STRUCTURAL_ALTERATION = 'CA', "complex_structural_alteration"
    COMPLEX_SUBSTITUTION = 'CS', "complex_substitution"
    COPY_NUMBER_GAIN = 'CG', "copy_number_gain"
    COPY_NUMBER_LOSS = 'CL', "copy_number_loss"
    COPY_NUMBER_VARIATION = 'CN', "copy_number_variation"
    DUPLICATION = 'DU', "duplication"
    INTERCHROMOSOMAL_BREAKPOINT = 'IB', "interchromosomal_breakpoint"
    INTERCHROMOSOMAL_TRANSLOCATION = 'IT', "interchromosomal_translocation"
    INTRACHROMOSOMAL_BREAKPOINT = 'CB', "intrachromosomal_breakpoint"
    INTRACHROMOSOMAL_TRANSLOCATION = 'CT', "intrachromosomal_translocation"
    INVERSION = 'IN', "inversion"
    LOSS_OF_HETEROZYGOSITY = 'LO', "loss_of_heterozygosity"
    MOBILE_ELEMENT_DELETION = 'MD', "mobile_element_deletion"
    MOBILE_ELEMENT_INSERTION = 'MI', "mobile_element_insertion"
    NOVEL_SEQUENCE_INSERTION = 'NI', "novel_sequence_insertion"
    SHORT_TANDEM_REPEAT_VARIATION = 'ST', "short_tandem_repeat_variation"
    TANDEM_DUPLICATION = 'TD', "tandem_duplication"
    TRANSLOCATION = 'TL', "translocation"
    DELETION = 'DE', "deletion"
    INDEL = 'ND', "indel"
    INSERTION = 'IS', "insertion"
    SEQUENCE_ALTERATION = 'SA', "sequence_alteration"
    PROBE = 'PR', "probe"


class ColumnAnnotationCategory(models.TextChoices):
    """ Based on categories from:
        https://asia.ensembl.org/info/docs/tools/vep/script/vep_plugins.html """

    CONSERVATION = 'C', "Conservation"
    EXTERNAL_ID = 'E', "External ID"
    FREQUENCY_DATA = 'F', "Frequency Data"
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
    MAXENTSCAN = 'm', 'MaxEntScan'
    NMD = "N", 'NMD'
    SPLICEAI = 'a', 'SpliceAI'
    SPLICEREGION = 's', 'SpliceRegion'


class VEPCustom(models.TextChoices):
    GNOMAD_2 = 'g', 'gnomAD2'
    GNOMAD_3 = 'n', 'gnomAD3'
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


class ClinVarReviewStatus(models.TextChoices):
    NO_ASSERTION_PROVIDED = "N", "No assertion provided"
    NO_ASSERTION_CRITERIA_PROVIDED = "C", "No assertion criteria provided"
    NO_INTERPRETATION_FOR_THE_SINGLE_VARIANT = "I", "No interpretation for the single variant"
    CRITERIA_PROVIDED_CONFLICTING_INTERPRETATIONS = "F", "Criteria provided - conflicting interpretations"
    CRITERIA_PROVIDED_SINGLE_SUBMITTER = "S", "Criteria provided - single submitter"
    CRITERIA_PROVIDED_MULTIPLE_SUBMITTERS_NO_CONFLICTS = "M", "Criteria provided - multiple submitters w/no conflicts"
    REVIEWED_BY_EXPERT_PANEL = "E", "Reviewed by expert panel"
    PRACTICE_GUIDELINE = "P", "Practice guideline"

    STARS = Constant({
        NO_ASSERTION_PROVIDED[0]: 0,
        NO_ASSERTION_CRITERIA_PROVIDED[0]: 0,
        NO_INTERPRETATION_FOR_THE_SINGLE_VARIANT[0]: 0,
        CRITERIA_PROVIDED_CONFLICTING_INTERPRETATIONS[0]: 1,
        CRITERIA_PROVIDED_SINGLE_SUBMITTER[0]: 1,
        CRITERIA_PROVIDED_MULTIPLE_SUBMITTERS_NO_CONFLICTS[0]: 2,
        REVIEWED_BY_EXPERT_PANEL[0]: 3,
        PRACTICE_GUIDELINE[0]: 4,
    })

    VCF_MAPPINGS = Constant({
        'no_assertion_provided': NO_ASSERTION_PROVIDED[0],
        'no_assertion_criteria_provided': NO_ASSERTION_CRITERIA_PROVIDED[0],
        'no_interpretation_for_the_single_variant': NO_INTERPRETATION_FOR_THE_SINGLE_VARIANT[0],
        'criteria_provided,_conflicting_interpretations': CRITERIA_PROVIDED_CONFLICTING_INTERPRETATIONS[0],
        'criteria_provided,_single_submitter': CRITERIA_PROVIDED_SINGLE_SUBMITTER[0],
        'criteria_provided,_multiple_submitters,_no_conflicts': CRITERIA_PROVIDED_MULTIPLE_SUBMITTERS_NO_CONFLICTS[0],
        'reviewed_by_expert_panel': REVIEWED_BY_EXPERT_PANEL[0],
        'practice_guideline': PRACTICE_GUIDELINE[0],
    })

    def stars(self):
        return ClinVarReviewStatus.STARS[self.value]

    @staticmethod
    def statuses_gte_stars(min_stars: int) -> List['ClinVarReviewStatus']:
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
