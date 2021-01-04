from django.db import models

from library.utils import Constant


class CaseState(models.TextChoices):
    OPEN = 'O', 'open'
    NO_TEST = 'N', 'No Test'
    CLOSED_SOLVED = 'S', 'solved'
    CLOSED_UNSOLVED = 'U', 'unsolved'

    CLOSED_STATES = Constant((CLOSED_SOLVED, CLOSED_UNSOLVED))


class CaseWorkflowStatus(models.TextChoices):
    NA = 'NA', "N/A"
    SAMPLE_PROCESSING = 'SP', "Sample Processing"
    LIBRARY_PREP_COMPLETE = "LP", "Library Prep Complete"
    SEQUENCING_COMPLETE = "SC", "Sequencing Complete"
    VCF_READY = "VR", "VCF Ready"
    ANALYSIS_COMPLETE = 'AC', "Analysis Complete"
    REPORTING = 'RP', "Reporting"


class InvestigationType(models.TextChoices):
    SINGLE_SAMPLE = 'S', 'Single Sample'
    TRIO = 'T', 'Trio'
    COHORT = 'C', 'Cohort'


class PathologyTestGeneModificationOutcome(models.TextChoices):
    PENDING = 'P', 'Pending'
    ACCEPTED = 'A', 'Accepted'
    REJECTED = 'R', 'Rejected'


class ClinicalSetting(models.TextChoices):
    DIAGNOSTIC_TEST = 'D', 'Diagnostic Test'
    PREDICTIVE_TEST = 'P', 'Predictive Test'
    CARRIER_TEST = 'C', 'Carrier Test'
    PRENATAL = 'N', 'Prenatal'


class PathologyTestType(models.TextChoices):
    COMMON_MUTATION_SCREEN = 'C', 'Common mutation screen'
    FULL_GENE_MUTATION_ANALYSIS = 'F', 'Full gene mutation analysis'
    KNOWN_FAMILIAL_MUTATIONS = 'K', 'Known familial mutation(s)'
