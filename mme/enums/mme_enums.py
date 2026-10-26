"""
MMESubmissionStatus: where an MMESubmission is in its life (draft, submitted to the remote node, or errored).
"""
from django.db import models


class MMESubmissionStatus(models.TextChoices):
    DRAFT = "D", "Draft"
    SUBMITTED = "S", "Submitted"
    ERROR = "E", "Error"
