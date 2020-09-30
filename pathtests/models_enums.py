class CaseState:
    OPEN = 'O'
    NO_TEST = 'N'
    CLOSED_SOLVED = 'S'
    CLOSED_UNSOLVED = 'U'

    CHOICES = (
        (OPEN, 'open'),
        (NO_TEST, 'No Test'),
        (CLOSED_SOLVED, 'solved'),
        (CLOSED_UNSOLVED, 'unsolved')
    )

    CLOSED_STATES = (CLOSED_SOLVED, CLOSED_UNSOLVED)


class CaseWorkflowStatus:
    NA = 'NA'
    SAMPLE_PROCESSING = 'SP'
    LIBRARY_PREP_COMPLETE = "LP"
    SEQUENCING_COMPLETE = "SC"
    VCF_READY = "VR"
    ANALYSIS_COMPLETE = 'AC'
    REPORTING = 'RP'

    CHOICES = (
        (NA, "N/A"),
        (SAMPLE_PROCESSING, "Sample Processing"),
        (LIBRARY_PREP_COMPLETE, "Library Prep Complete"),
        (SEQUENCING_COMPLETE, "Sequencing Complete"),
        (VCF_READY, "VCF Ready"),
        (ANALYSIS_COMPLETE, "Analysis Complete"),
        (REPORTING, "Reporting"),
    )


class InvestigationType:
    SINGLE_SAMPLE = 'S'
    TRIO = 'T'
    COHORT = 'C'

    CHOICES = (
        (SINGLE_SAMPLE, 'Single Sample'),
        (TRIO, 'Trio'),
        (COHORT, 'Cohort')
    )


class PathologyTestGeneModificationOutcome:
    PENDING = 'P'
    ACCEPTED = 'A'
    REJECTED = 'R'
    CHOICES = (
        (PENDING, 'Pending'),
        (ACCEPTED, 'Accepted'),
        (REJECTED, 'Rejected'),
    )


class ClinicalSetting:
    DIAGNOSTIC_TEST = 'D'
    PREDICTIVE_TEST = 'P'
    CARRIER_TEST = 'C'
    PRENATAL = 'N'

    CHOICES = (
        ('', 'n/a'),
        (DIAGNOSTIC_TEST, 'Diagnostic Test'),
        (PREDICTIVE_TEST, 'Predictive Test'),
        (CARRIER_TEST, 'Carrier Test'),
        (PRENATAL, 'Prenatal'),
    )


class PathologyTestType:
    COMMON_MUTATION_SCREEN = 'C'
    FULL_GENE_MUTATION_ANALYSIS = 'F'
    KNOWN_FAMILIAL_MUTATIONS = 'K'

    CHOICES = (
        ('', 'n/a'),
        (COMMON_MUTATION_SCREEN, 'Common mutation screen'),
        (FULL_GENE_MUTATION_ANALYSIS, 'Full gene mutation analysis'),
        (KNOWN_FAMILIAL_MUTATIONS, 'Known familial mutation(s)'),
    )
