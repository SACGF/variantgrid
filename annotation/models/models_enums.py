class HumanProteinAtlasAbundance:
    NOT_DETECTED = 'N'
    LOW = 'L'
    MEDIUM = 'M'
    HIGH = 'H'
    CHOICES = [
        (NOT_DETECTED, 'Not detected'),
        (LOW, 'Low'),
        (MEDIUM, 'Medium'),
        (HIGH, 'High'),
    ]
    DETECTED_CHOICES = [
        (LOW, 'Low'),
        (MEDIUM, 'Medium'),
        (HIGH, 'High'),
    ]

    @staticmethod
    def get_abundance_or_greater_levels(abundance):
        higher_levels = []
        found = False
        for (k, _) in HumanProteinAtlasAbundance.CHOICES:
            if k == abundance:
                found = True
            if found:
                higher_levels.append(k)
        return higher_levels


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


class AnnotationStatus:
    CREATED = 'C'
    DELETING = 'x'
    DUMP_STARTED = 'd'
    DUMP_COMPLETED = 'D'
    ANNOTATION_STARTED = 'a'
    ANNOTATION_COMPLETED = 'A'
    UPLOAD_STARTED = 'U'
    FINISHED = 'F'
    ERROR = 'E'
    CHOICES = [
        (CREATED, "Created"),
        (DELETING, "Deleting"),
        (DUMP_STARTED, "Dump Started"),
        (DUMP_COMPLETED, "Dump Completed"),
        (ANNOTATION_STARTED, "Annotation Started"),
        (ANNOTATION_COMPLETED, "Annotation Completed"),
        (UPLOAD_STARTED, "Upload Started"),
        (FINISHED, "Finished"),
        (ERROR, "Error"),
    ]

    @classmethod
    def get_summary_state(cls, annotation_status):
        SUMMARY_STATES = {cls.CREATED: "Queued",
                          cls.FINISHED: "Finished",
                          cls.ERROR: "Error"}
        return SUMMARY_STATES.get(annotation_status, "Running")


class HPOSynonymScope:
    BROAD = 'B'
    EXACT = 'E'
    NARROW = 'N'
    RELATED = 'R'

    CHOICES = [
        (BROAD, 'Broad'),
        (EXACT, 'Exact'),
        (NARROW, 'Narrow'),
        (RELATED, 'Related'),
    ]


class CitationSource:
    PUBMED = 'P'
    NCBI_BOOKSHELF = 'N'
    PUBMED_CENTRAL = 'C'
    CHOICES = (
        (PUBMED, 'PubMed'),
        (NCBI_BOOKSHELF, 'NCBIBookShelf'),
        (PUBMED_CENTRAL, 'PubMedCentral'),
    )

    CODES = {'PubMed': PUBMED,
             'PMID': PUBMED,
             'NCBIBookShelf': NCBI_BOOKSHELF,
             'PubMedCentral': PUBMED_CENTRAL}


class TranscriptStatus:
    KNOWN = 'K'
    NOVEL = 'N'
    PUTATIVE = 'P'
    CHOICES = (
        (KNOWN, "KNOWN"),
        (NOVEL, "NOVEL"),
        (PUTATIVE, "PUTATIVE")
    )


class GenomicStrand:
    SENSE = '+'
    ANTISENSE = '-'
    CHOICES = (
        (SENSE, "+"),
        (ANTISENSE, "-")
    )


class ClinGenClassification:
    DEFINITIVE = 'D'
    STRONG = 'S'
    MODERATE = 'M'
    LIMITED = 'L'
    NO_REPORTED_EVIDENCE = 'N'
    REFUTED = 'R'
    DISPUTED = 'P'

    CHOICES = (
        (DEFINITIVE, 'Definitive'),
        (STRONG, 'Strong'),
        (MODERATE, 'Moderate'),
        (LIMITED, 'Limited'),
        (NO_REPORTED_EVIDENCE, 'No Reported Evidence'),
        (REFUTED, 'Refuted'),
        (DISPUTED, 'Disputed'),
    )


class VariantClass:
    """ https://asia.ensembl.org/info/genome/variation/prediction/classification.html#classes """

    SNV = 'SN'
    GENETIC_MARKER = 'GM'
    SUBSTITUTION = 'SU'
    TANDEM_REPEAT = 'TR'
    ALU_INSERTION = 'AI'
    COMPLEX_STRUCTURAL_ALTERATION = 'CA'
    COMPLEX_SUBSTITUTION = 'CS'
    COPY_NUMBER_GAIN = 'CG'
    COPY_NUMBER_LOSS = 'CL'
    COPY_NUMBER_VARIATION = 'CN'
    DUPLICATION = 'DU'
    INTERCHROMOSOMAL_BREAKPOINT = 'IB'
    INTERCHROMOSOMAL_TRANSLOCATION = 'IT'
    INTRACHROMOSOMAL_BREAKPOINT = 'CB'
    INTRACHROMOSOMAL_TRANSLOCATION = 'CT'
    INVERSION = 'IN'
    LOSS_OF_HETEROZYGOSITY = 'LO'
    MOBILE_ELEMENT_DELETION = 'MD'
    MOBILE_ELEMENT_INSERTION = 'MI'
    NOVEL_SEQUENCE_INSERTION = 'NI'
    SHORT_TANDEM_REPEAT_VARIATION = 'ST'
    TANDEM_DUPLICATION = 'TD'
    TRANSLOCATION = 'TL'
    DELETION = 'DE'
    INDEL = 'ND'
    INSERTION = 'IS'
    SEQUENCE_ALTERATION = 'SA'
    PROBE = 'PR'

    CHOICES = [
        (SNV, "SNV"),
        (GENETIC_MARKER, "genetic_marker"),
        (SUBSTITUTION, "substitution"),
        (TANDEM_REPEAT, "tandem_repeat"),
        (ALU_INSERTION, "Alu_insertion"),
        (COMPLEX_STRUCTURAL_ALTERATION, "complex_structural_alteration"),
        (COMPLEX_SUBSTITUTION, "complex_substitution"),
        (COPY_NUMBER_GAIN, "copy_number_gain"),
        (COPY_NUMBER_LOSS, "copy_number_loss"),
        (COPY_NUMBER_VARIATION, "copy_number_variation"),
        (DUPLICATION, "duplication"),
        (INTERCHROMOSOMAL_BREAKPOINT, "interchromosomal_breakpoint"),
        (INTERCHROMOSOMAL_TRANSLOCATION, "interchromosomal_translocation"),
        (INTRACHROMOSOMAL_BREAKPOINT, "intrachromosomal_breakpoint"),
        (INTRACHROMOSOMAL_TRANSLOCATION, "intrachromosomal_translocation"),
        (INVERSION, "inversion"),
        (LOSS_OF_HETEROZYGOSITY, "loss_of_heterozygosity"),
        (MOBILE_ELEMENT_DELETION, "mobile_element_deletion"),
        (MOBILE_ELEMENT_INSERTION, "mobile_element_insertion"),
        (NOVEL_SEQUENCE_INSERTION, "novel_sequence_insertion"),
        (SHORT_TANDEM_REPEAT_VARIATION, "short_tandem_repeat_variation"),
        (TANDEM_DUPLICATION, "tandem_duplication"),
        (TRANSLOCATION, "translocation"),
        (DELETION, "deletion"),
        (INDEL, "indel"),
        (INSERTION, "insertion"),
        (SEQUENCE_ALTERATION, "sequence_alteration"),
        (PROBE, "probe"),
    ]


class ColumnAnnotationCategory:
    """ Based on categories from:
        https://asia.ensembl.org/info/docs/tools/vep/script/vep_plugins.html """

    CONSERVATION = 'C'
    EXTERNAL_ID = 'E'
    FREQUENCY_DATA = 'F'
    GENE_ANNOTATIONS = 'G'
    HGVS = 'H'
    LITERATURE = 'L'
    NEARBY_FEATURES = 'N'
    PATHOGENICITY_PREDICTIONS = 'P'
    PHENOTYPE = 'Y'
    PROTEIN_DOMAINS = 'D'
    SEQUENCE = 'Q'
    SPLICING_PREDICTIONS = 'S'
    VARIANT_DATA = 'V'

    CHOICES = [
        (CONSERVATION, "Conservation"),
        (EXTERNAL_ID, "External ID"),
        (FREQUENCY_DATA, "Frequency Data"),
        (GENE_ANNOTATIONS, 'Gene Annotations'),
        (HGVS, "HGVS"),
        (LITERATURE, 'Literature'),
        (NEARBY_FEATURES, "Nearby Features"),
        (PATHOGENICITY_PREDICTIONS, "Pathogenicity Predictions"),
        (PHENOTYPE, "Phenotype"),
        (PROTEIN_DOMAINS, "Protein Domains"),
        (SEQUENCE, "Sequence"),
        (SPLICING_PREDICTIONS, "Splicing Predictions"),
        (VARIANT_DATA, "Variant Data"),
    ]


class VEPPlugin:
    DBNSFP = 'd'
    DBSCSNV = 'v'
    GRANTHAM = 'g'
    LOFTOOL = 'l'
    MASTERMIND = 'n'
    MAXENTSCAN = 'm'
    SPLICEAI = 'a'
    SPLICEREGION = 's'

    CHOICES = [
        (DBNSFP, 'dbNSFP'),
        (DBSCSNV, 'dbscSNV'),
        (GRANTHAM, 'Grantham'),
        (LOFTOOL, 'LoFtool'),
        (MASTERMIND, 'Mastermind'),
        (MAXENTSCAN, 'MaxEntScan'),
        (SPLICEAI, 'SpliceAI'),
        (SPLICEREGION, 'SpliceRegion'),
    ]


class VEPCustom:
    GNOMAD = 'g'
    PHASTCONS_100_WAY = '1'
    PHASTCONS_30_WAY = '2'
    PHASTCONS_46_WAY = '3'
    PHYLOP_100_WAY = '4'
    PHYLOP_30_WAY = '5'
    PHYLOP_46_WAY = '6'
    REPEAT_MASKER = 'r'
    TOPMED = 't'
    UK10K = 'u'
    COSMIC = 'c'

    CHOICES = [
        (GNOMAD, 'gnomAD'),
        (PHASTCONS_100_WAY, 'phastCons100way'),
        (PHASTCONS_30_WAY, 'phastCons30way'),
        (PHASTCONS_46_WAY, 'phastCons46way'),
        (PHYLOP_100_WAY, 'phyloP100way'),
        (PHYLOP_30_WAY, 'phyloP30way'),
        (PHYLOP_46_WAY, 'phyloP46way'),
        (REPEAT_MASKER, 'RepeatMasker'),
        (TOPMED, 'TopMed'),
        (UK10K, 'UK10k'),
        (COSMIC, 'COSMIC'),
    ]


class VEPSkippedReason:
    UNKNOWN_CONTIG = 'c'
    INCOMPLETE = 'i'
    UNKNOWN = 'u'

    CHOICES = [
        (UNKNOWN_CONTIG, "Unknown Contig"),
        (INCOMPLETE, "Incomplete"),
        (UNKNOWN, "Unknown"),
    ]


class ClinVarReviewStatus:
    NO_ASSERTION_PROVIDED = "N"
    NO_ASSERTION_CRITERIA_PROVIDED = "C"
    NO_INTERPRETATION_FOR_THE_SINGLE_VARIANT = "I"
    CRITERIA_PROVIDED_CONFLICTING_INTERPRETATIONS = "C"
    CRITERIA_PROVIDED_SINGLE_SUBMITTER = "S"
    CRITERIA_PROVIDED_MULTIPLE_SUBMITTERS_NO_CONFLICTS = "M"
    REVIEWED_BY_EXPERT_PANEL = "E"
    PRACTICE_GUIDELINE = "P"

    CHOICES = [
        (NO_ASSERTION_PROVIDED, "No assertion provided"),
        (NO_ASSERTION_CRITERIA_PROVIDED, "No assertion criteria provided"),
        (NO_INTERPRETATION_FOR_THE_SINGLE_VARIANT, "No interpretation for the single variant"),
        (CRITERIA_PROVIDED_CONFLICTING_INTERPRETATIONS, "Criteria provided - conflicting interpretations"),
        (CRITERIA_PROVIDED_SINGLE_SUBMITTER, "Criteria provided - single submitter"),
        (CRITERIA_PROVIDED_MULTIPLE_SUBMITTERS_NO_CONFLICTS, "Criteria provided - multiple submitters w/no conflicts"),
        (REVIEWED_BY_EXPERT_PANEL, "Reviewed by expert panel"),
        (PRACTICE_GUIDELINE, "Practice guideline"),
    ]

    STARS = {
        NO_ASSERTION_PROVIDED: 0,
        NO_ASSERTION_CRITERIA_PROVIDED: 0,
        NO_INTERPRETATION_FOR_THE_SINGLE_VARIANT: 0,
        CRITERIA_PROVIDED_CONFLICTING_INTERPRETATIONS: 1,
        CRITERIA_PROVIDED_SINGLE_SUBMITTER: 1,
        CRITERIA_PROVIDED_MULTIPLE_SUBMITTERS_NO_CONFLICTS: 2,
        REVIEWED_BY_EXPERT_PANEL: 3,
        PRACTICE_GUIDELINE: 4,
    }

    VCF_MAPPINGS = {
        'no_assertion_provided': NO_ASSERTION_PROVIDED,
        'no_assertion_criteria_provided': NO_ASSERTION_CRITERIA_PROVIDED,
        'no_interpretation_for_the_single_variant': NO_INTERPRETATION_FOR_THE_SINGLE_VARIANT,
        'criteria_provided,_conflicting_interpretations': CRITERIA_PROVIDED_CONFLICTING_INTERPRETATIONS,
        'criteria_provided,_single_submitter': CRITERIA_PROVIDED_SINGLE_SUBMITTER,
        'criteria_provided,_multiple_submitters,_no_conflicts': CRITERIA_PROVIDED_MULTIPLE_SUBMITTERS_NO_CONFLICTS,
        'reviewed_by_expert_panel': REVIEWED_BY_EXPERT_PANEL,
        'practice_guideline': PRACTICE_GUIDELINE,
    }
