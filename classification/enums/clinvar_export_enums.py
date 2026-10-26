"""
ClinVar export vocabulary: ClinVarExportStatus (set by the converter,
classification/models/clinvar_export_convertor.py, which loads before the export models), batch,
submission and request statuses, the ClinVar API response outcomes, and the per-ClinVarKey options
(assertion methods, citation modes, exclude pattern modes).
"""
from enum import Enum, auto

from django.db import models
from django.db.models import TextChoices


class ClinVarExportStatus(models.TextChoices):
    """
    The status of an Export as compared to the last batch.
    Note doesn't involve itself in the status as data progresses
    """

    NEW_SUBMISSION = "N", "New Submission"  # new submission and changes pending often work the same, but might be useful to see at a glance, useful if we do approvals
    CHANGES_PENDING = "C", "Changes Pending"
    UP_TO_DATE = "D", "Up to Date"
    IN_ERROR = "E", "Error"
    EXCLUDE = "X", "Exclude"


class ClinVarKeyExcludePatternMode(TextChoices):
    EXCLUDE_IF_MATCH = 'X'
    EXCLUDE_IF_NO_MATCH = 'N'


class ClinVarExportBatchStatus(models.TextChoices):
    """
    Overall status of a submission batch, talking to ClinVar takes a few steps
    """
    AWAITING_UPLOAD = "W", "Awaiting Upload"  # new submission and changes pending often work the same, but might be useful to see at a glance, useful if we do approvals
    UPLOADING = "P", "Processing"
    SUBMITTED = "S", "Submitted"
    REJECTED = "R", "Rejected"


class ClinVarExportSubmissionStatus(TextChoices):
    WAITING = "W", "Waiting"
    SUCCESS = "S", "Success"
    ERROR = "E", "Error"


class ClinVarExportRequestType(TextChoices):
    INITIAL_SUBMISSION = "S", "Initial Submission"
    POLLING_SUBMISSION = "P", "Polling"
    RESPONSE_FILES = "F", "Response Files"


class ClinVarRequestExceptionType(Enum):
    OUR_DATA_ISSUE = auto()
    API_KEY_ISSUE = auto()
    TOO_MANY_REQUESTS = auto()
    CLINVAR_SERVER_ISSUE = auto()
    RESPONSE_MISSING_KEY_DATA = auto()
    NOT_SUPPORTED_YET = auto()


class ClinVarResponseOutcome(Enum):
    ASK_AGAIN_NOW = auto()  # if we can immediately do the next thing in the queue for this batch
    ASK_AGAIN_LATER = auto()
    COMPLETE = auto()


class ClinVarAssertionMethods(TextChoices):
    # "yes", "no", "unknown", "not provided", "not applicable"
    yes = "yes"
    no = "no"
    unknown = "unknown"
    not_provided = "not provided"
    not_applicable = "not applicable"


class ClinVarCitationsModes(TextChoices):
    all = "all"
    interpretation_summary_only = "interpret"
