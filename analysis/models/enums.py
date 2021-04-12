import operator
from functools import reduce

from django.db import models

from library.utils import Constant
from patients.models_enums import Zygosity


class AnalysisType(models.TextChoices):
    SINGLETON = 'S', 'Singleton'
    COHORT = 'C', 'Cohort'
    TRIO = 'T', 'Trio'
    PEDIGREE = 'P', 'Pedigree'


class AnalysisTemplateType(models.TextChoices):
    TEMPLATE = 'T', 'Template'
    SNAPSHOT = 'S', 'Snapshot'


class SetOperations(models.TextChoices):
    NONE = '0', 'none'
    A_NOT_B = 'l', 'A-B'
    INTERSECTION = '&', '&'
    A_ONLY = 'a', 'A'
    B_NOT_A = 'r', 'B-A'
    SYMMETRIC_DIFFERENCE = '^', "^"  # XOR (exclusive or)
    B_ONLY = 'b', 'B'
    UNION = '|', "|"


class SNPMatrix(models.TextChoices):
    TOTAL_PERCENT = 't', 'total percent'
    ROWS_PERCENT = 'r', 'percent of rows'
    COLS_PERCENT = 'c', 'percent of cols'


class TrioInheritance(models.TextChoices):
    RECESSIVE = 'R', 'Recessive'
    COMPOUND_HET = 'C', 'C. Het'
    DOMINANT = 'D', 'Dominant'
    DENOVO = 'N', "Denovo"
    PROBAND_HET = 'P', 'Proband Het'  # can be fed into compound het (so analysis template only takes Trio variable)
    XLINKED_RECESSIVE = 'X', "X-Linked Recessive"


class GroupOperation(models.TextChoices):
    ALL = "L", "All"  # AND
    ANY = "Y", "Any"  # OR

    @classmethod
    def get_operation(cls, group_operation):
        OPERATIONS = {
            cls.ALL: operator.or_,
            cls.ANY: operator.and_,
        }
        return OPERATIONS[group_operation]

    @staticmethod
    def reduce(sequence, group_operation):
        if group_operation == GroupOperation.ALL:
            return reduce(operator.and_, sequence)
        return reduce(operator.or_, sequence)


class MinimisationResultType(models.TextChoices):
    FULL = 'F', "Full"
    BOOTSTRAPPED = 'B', "Bootstrapped"


class MinimisationStrategy(models.TextChoices):
    LEAST_SQUARES = "LS", "Least Squares"
    LEAST_ABSOLUTE = "LA", "Least Absolute Values"


class NodeErrorSource(models.TextChoices):
    ANALYSIS = 'A', 'Analysis'
    CONFIGURATION = 'C', 'Configuration'
    INTERNAL_ERROR = 'I', 'Internal Error'
    PARENT = 'P', 'Parent'


class NodeStatus(models.TextChoices):
    ERROR = 'E', 'Error'
    ERROR_CONFIGURATION = 'C', 'Configuration Error'
    ERROR_WITH_PARENT = 'P', 'Parent Error'
    CANCELLED = 'N', 'Cancelled'
    READY = 'R', 'Ready'
    DIRTY = 'D', 'Dirty'
    QUEUED = 'Q', 'Queued'
    LOADING_CACHE = 'H', 'Loading Cache'
    LOADING = 'L', 'Loading'

    LOADING_STATUSES = Constant([e[0] for e in (DIRTY, QUEUED, LOADING_CACHE, LOADING)])
    ERROR_STATUSES = Constant([e[0] for e in (ERROR, ERROR_CONFIGURATION, ERROR_WITH_PARENT, CANCELLED)])
    READY_STATUSES = Constant([e[0] for e in (ERROR, ERROR_CONFIGURATION, ERROR_WITH_PARENT, CANCELLED, READY)])

    @staticmethod
    def is_loading(status):
        return status in NodeStatus.LOADING_STATUSES

    @staticmethod
    def is_ready(status):
        return status in NodeStatus.READY_STATUSES

    @staticmethod
    def is_error(status):
        return status in NodeStatus.ERROR_STATUSES

    @classmethod
    def get_summary_state(cls, status):
        summary_state = {cls.READY: "Finished"}
        for s in [cls.DIRTY, cls.QUEUED]:
            summary_state[s] = "Queued"

        for s in [cls.LOADING_CACHE, cls.LOADING]:
            summary_state[s] = "Running"

        for s in cls.ERROR_STATUSES:
            summary_state[s] = "User Error"
        summary_state[cls.ERROR] = "Error"

        return summary_state.get(status)


class NodeColors:
    VALID = None
    WARNING = '#FFA500'
    ERROR = '#ee0000'


class ZygosityNodeZygosity(Zygosity):
    MULTIPLE_HIT = 'M'  # >=2 hits, regardless of zygosity
    CHOICES = Zygosity.CHOICES + [(MULTIPLE_HIT, "Multiple hits in gene")]


class TrioSample(models.TextChoices):
    MOTHER = 'M', 'Mother'
    FATHER = 'F', 'Father'
    PROBAND = 'P', 'Proband'


class TagNodeMode(models.TextChoices):
    PARENT = 'P', 'Parent'
    THIS_ANALYSIS = 'T', 'This Analysis'
    ALL_TAGS = 'L', 'All Tags'


class TagLocation(models.TextChoices):
    ANALYSIS = 'A', 'Analysis'
    EXTERNAL = 'E', 'External Import'
    GENE_PAGE = 'G', 'Gene Page'
    VARIANT_DETAILS = 'V', 'Variant Details'
