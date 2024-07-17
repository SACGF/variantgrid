# Using own file to stop circular dependencies in illuminate_report
from django.db import models

from library.utils import Constant


class QCGraphTypes2(models.TextChoices):
    BOX_PLOT = 'B', 'Box Plot'
    BAR_CHART = 'A', 'Bar Chart'


class QCCompareType(models.TextChoices):
    GOLD_ENRICHMENT_KIT_RUNS = 'G', 'Gold Enrichment Kit Runs'
    ENRICHMENT_KIT = 'P', 'Enrichment Kit'
    SEQUENCER = 'S', 'Sequencer'
    ALL_RUNS = 'A', 'All Runs'
    SEQUENCING_RUN = 'R', 'SequencingRun'

    ENRICHMENT_KIT_TYPES = Constant((GOLD_ENRICHMENT_KIT_RUNS[0], ENRICHMENT_KIT[0]))


class SequencingFileType(models.TextChoices):
    SAMPLE_SHEET = 'S', 'SampleSheet'
    ILLUMINA_FLOWCELL_QC = 'I', 'Illumina_Flowcell_QC'
    FASTQ = 'U', 'FastQ'  # Unmapped
    FASTQC = 'F', 'FastQC'
    BAM = 'B', 'Bam'
    FLAGSTATS = 'T', 'Flagstats'
    VCF = 'V', 'VCF'
    COMBINED_VCF = 'C', 'CombinedVCF'
    QC = 'Q', 'QC'
    DATA_MIGRATION = 'M', "Data Migration"


class JobScriptStatus(models.TextChoices):
    CREATED = 'C', 'Created'
    SUBMITTED = 'S', 'Submitted'
    FINISHED = 'F', 'Finished'


class DataGeneration(models.TextChoices):
    HISEQ = 'H', 'HiSeq'
    MISEQ = 'M', 'MiSeq'


class PairedEnd(models.TextChoices):
    R1 = "R1"
    R2 = "R2"


class SequencerRead(models.TextChoices):
    R1 = "R1"
    R2 = "R2"
    I1 = "I1"
    I2 = "I2"


class SeqAutoRunStatus(models.TextChoices):
    CREATED = 'c', "Created"
    SCANNING_FILES = 's', "Scanning Files"
    CREATE_MODELS = 'm', "Create Models"
    SCRIPTS_AND_JOBS = 'j', "Scripts and Jobs"
    FINISHED = 'F', "Finished"
    ERROR = 'E', "Error"

    COMPLETED_STATES = Constant([FINISHED[0], ERROR[0]])


class EnrichmentKitType(models.TextChoices):
    AMPLICON = 'A', 'Amplicon'
    CAPTURE = 'C', 'Capture'
    VIRTUAL = 'V', 'Virtual'
