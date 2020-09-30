# Using own file to stop circular dependencies in illuminate_report


class DataState:
    NON_EXISTENT = 'N'
    DELETED = 'D'
    RUNNING = 'R'
    SKIPPED = 'S'
    ERROR = 'E'
    COMPLETE = 'C'

    CHOICES = (
        (NON_EXISTENT, 'Non Existent'),
        (DELETED, 'Deleted'),
        (RUNNING, 'Running'),
        (SKIPPED, 'Skipped'),
        (ERROR, 'Error'),
        (COMPLETE, 'Complete'),
    )

    @staticmethod
    def should_create_new_record(data_state):
        return data_state not in [DataState.DELETED, DataState.SKIPPED]


class QCGraphEnrichmentKitSeparationChoices:
    ALL_ENRICHMENT_KITS = 'A'
    SEPARATED_ENRICHMENT_KITS = 'S'
    SELECTED_ENRICHMENT_KIT = 'O'

    CHOICES = (
        (ALL_ENRICHMENT_KITS, 'All Enrichment Kits'),
        (SEPARATED_ENRICHMENT_KITS, 'Separated Enrichment Kit'),
        (SELECTED_ENRICHMENT_KIT, 'User selected Enrichment Kit'),
    )

    SHOW_ENRICHMENT_KIT = {SELECTED_ENRICHMENT_KIT: True}


class QCGraphType:
    LINE_GRAPH = 'L'
    BOX_PLOT = 'B'

    CHOICES = (
        (LINE_GRAPH, 'Line Graph'),
        (BOX_PLOT, 'Box Plot'),
    )


class QCGraphTypes2:
    BOX_PLOT = 'B'
    BAR_CHART = 'A'

    CHOICES = (
        (BOX_PLOT, 'Box Plot'),
        (BAR_CHART, 'Bar Chart'),
    )


class QCCompareType:
    GOLD_ENRICHMENT_KIT_RUNS = 'G'
    ENRICHMENT_KIT = 'P'
    SEQUENCER = 'S'
    ALL_RUNS = 'A'
    SEQUENCING_RUN = 'R'

    ENRICHMENT_KIT_TYPES = (GOLD_ENRICHMENT_KIT_RUNS, ENRICHMENT_KIT)
    CHOICES = (
        (GOLD_ENRICHMENT_KIT_RUNS, 'Gold Enrichment Kit Runs'),
        (ENRICHMENT_KIT, 'Enrichment Kit'),
        (SEQUENCER, 'Sequencer'),
        (ALL_RUNS, 'All Runs'),
        (SEQUENCING_RUN, 'SequencingRun'),
    )


class SequencingFileType:
    SAMPLE_SHEET = 'S'
    ILLUMINA_FLOWCELL_QC = 'I'
    FASTQ = 'U'  # Unmapped
    FASTQC = 'F'
    BAM = 'B'
    FLAGSTATS = 'T'
    VCF = 'V'
    COMBINED_VCF = 'C'
    QC = 'Q'
    DATA_MIGRATION = 'M'

    CHOICES = (
        (SAMPLE_SHEET, 'SampleSheet'),
        (ILLUMINA_FLOWCELL_QC, 'Illumina_Flowcell_QC'),
        (FASTQ, 'FastQ'),
        (FASTQC, 'FastQC'),
        (BAM, 'Bam'),
        (FLAGSTATS, 'Flagstats'),
        (VCF, 'VCF'),
        (COMBINED_VCF, 'CombinedVCF'),
        (QC, 'QC'),
        (DATA_MIGRATION, "Data Migration")
    )


class JobScriptStatus:
    CREATED = 'C'
    SUBMITTED = 'S'
    FINISHED = 'F'

    CHOICES = (
        (CREATED, 'Created'),
        (SUBMITTED, 'Submitted'),
        (FINISHED, 'Finished'),
    )


class DataGeneration:
    HISEQ = 'H'
    MISEQ = 'M'

    CHOICES = (
        (HISEQ, 'HiSeq'),
        (MISEQ, 'MiSeq'),
    )


class PairedEnd:
    R1 = "R1"
    R2 = "R2"
    CHOICES = ((R1, "R1"),
               (R2, "R2"))


class SequencerRead:
    I1 = "I1"
    I2 = "I2"
    CHOICES = PairedEnd.CHOICES + ((I1, I1), (I2, I2))


class SeqAutoRunStatus:
    CREATED = 'c'
    SCANNING_FILES = 's'
    CREATE_MODELS = 'm'
    SCRIPTS_AND_JOBS = 'j'
    FINISHED = 'F'
    ERROR = 'E'
    CHOICES = [
        (CREATED, "Created"),
        (SCANNING_FILES, "Scanning Files"),
        (CREATE_MODELS, "Create Models"),
        (FINISHED, "Finished"),
        (ERROR, "Error"),
    ]

    COMPLETED_STATES = [FINISHED, ERROR]


class EnrichmentKitType:
    AMPLICON = 'A'
    CAPTURE = 'C'
    VIRTUAL = 'V'

    CHOICES = ((AMPLICON, 'Amplicon'),
               (CAPTURE, 'Capture'),
               (VIRTUAL, 'Virtual'))
