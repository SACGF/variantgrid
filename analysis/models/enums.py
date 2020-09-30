import operator
from functools import reduce

from patients.models_enums import Zygosity


class AnalysisType:
    SINGLETON = 'S'
    COHORT = 'C'
    TRIO = 'T'
    PEDIGREE = 'P'

    CHOICES = (
        (SINGLETON, 'Singleton'),
        (COHORT, 'Cohort'),
        (TRIO, 'Trio'),
        (PEDIGREE, 'Pedigree'),
    )


class AnalysisTemplateType:
    TEMPLATE = 'T'
    SNAPSHOT = 'S'

    CHOICES = (
        (TEMPLATE, 'Template'),
        (SNAPSHOT, 'Snapshot'),
    )


class SetOperations:
    NONE = '0'
    A_NOT_B = 'l'
    INTERSECTION = '&'
    A_ONLY = 'a'
    B_NOT_A = 'r'
    SYMMETRIC_DIFFERENCE = '^'
    B_ONLY = 'b'
    UNION = '|'
    CHOICES = (
        (NONE, 'none'),
        (A_NOT_B, 'A-B'),
        (INTERSECTION, '&'),
        (A_ONLY, 'A'),
        (B_NOT_A, 'B-A'),
        (SYMMETRIC_DIFFERENCE, "^"),  # XOR (exclusive or)
        (B_ONLY, 'B'),
        (UNION, "|"),
    )


class SNPMatrix:
    TOTAL_PERCENT = 't'
    ROWS_PERCENT = 'r'
    COLS_PERCENT = 'c'
    CHOICES = (
        (TOTAL_PERCENT, 'total percent'),
        (ROWS_PERCENT, 'percent of rows'),
        (COLS_PERCENT, 'percent of cols'),
    )


class TrioInheritance:
    RECESSIVE = 'R'
    COMPOUND_HET = 'C'
    DOMINANT = 'D'
    DENOVO = 'N'
    PROBAND_HET = 'P'  # Proband HET can be fed into compound het (so analysis template only takes Trio variable)
    XLINKED_RECESSIVE = 'X'

    CHOICES = (
        (RECESSIVE, 'Recessive'),
        (COMPOUND_HET, 'C. Het'),
        (DOMINANT, 'Dominant'),
        (DENOVO, "Denovo"),
        (PROBAND_HET, 'Proband Het'),
        (XLINKED_RECESSIVE, "X-Linked Recessive"),
    )


class GroupOperation:
    ALL = "L"
    ANY = "Y"
    CHOICES = [
        (ALL, "All"),  # AND
        (ANY, "Any"),  # OR
    ]
    OPERATIONS = {
        ALL: operator.or_,
        ANY: operator.and_,
    }

    @staticmethod
    def reduce(sequence, group_operation):
        if group_operation == GroupOperation.ALL:
            return reduce(operator.and_, sequence)
        return reduce(operator.or_, sequence)


class MinimisationResultType:
    FULL = 'F'
    BOOTSTRAPPED = 'B'
    CHOICES = [
        (FULL, "Full"),
        (BOOTSTRAPPED, "Bootstrapped"),
    ]


class MinimisationStrategy:
    LEAST_SQUARES = "LS"
    LEAST_ABSOLUTE = "LA"

    CHOICES = [
        (LEAST_SQUARES, "Least Squares"),
        (LEAST_ABSOLUTE, "Least Absolute Values"),
    ]


class NodeErrorSource:
    ANALYSIS = 'A'
    CONFIGURATION = 'C'
    INTERNAL_ERROR = 'I'
    PARENT = 'P'

    CHOICES = (
        (ANALYSIS, 'Analysis'),
        (CONFIGURATION, 'Configuration'),
        (INTERNAL_ERROR, 'Internal Error'),
        (PARENT, 'Parent'),
    )


class NodeStatus:
    ERROR = 'E'
    ERROR_CONFIGURATION = 'C'
    ERROR_WITH_PARENT = 'P'
    CANCELLED = 'N'
    READY = 'R'
    DIRTY = 'D'
    QUEUED = 'Q'
    LOADING_CACHE = 'H'
    LOADING = 'L'

    CHOICES = (
        (ERROR, 'Error'),
        (ERROR_CONFIGURATION, 'Configuration Error'),
        (ERROR_WITH_PARENT, 'Parent Error'),
        (CANCELLED, 'Cancelled'),
        (READY, 'Ready'),
        (DIRTY, 'Dirty'),
        (QUEUED, 'Queued'),
        (LOADING_CACHE, 'Loading Cache'),
        (LOADING, 'Loading'),
    )
    LOADING_STATUSES = (DIRTY, QUEUED, LOADING_CACHE, LOADING)
    ERROR_STATUSES = (ERROR, ERROR_CONFIGURATION, ERROR_WITH_PARENT, CANCELLED)
    READY_STATUSES = (READY,) + ERROR_STATUSES  # Ready = Node load

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
