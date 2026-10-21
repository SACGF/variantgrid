"""
Statuses for classification imports: ClassificationImportRunStatus (an API import run) and
UploadedClassificationsUnmappedStatus (an uploaded file going through mapping and import).
"""
from django.db import models


class ClassificationImportRunStatus(models.TextChoices):
    ONGOING = "O", 'Ongoing'
    COMPLETED = "C", 'Completed'
    UNFINISHED = "U", 'Unfinished'


class UploadedClassificationsUnmappedStatus(models.TextChoices):
    Manual = "MA", "Manual Review Pending"
    Pending = 'P', 'Pending'

    Mapping = "M", "Mapping"
    Importing = "I", "Importing"

    Validated = 'V', 'Validated'
    Error = 'E', 'Error'
    Processed = 'MP', 'Processed'
