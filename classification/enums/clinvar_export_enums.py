"""
ClinVarExportStatus, kept apart from the ClinVar export models so the converter that sets it
(classification/models/clinvar_export_convertor.py) can load before them.
"""
from django.db import models


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
