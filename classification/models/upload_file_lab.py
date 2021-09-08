from django.contrib.auth.models import User
from django.db import models
from model_utils.models import TimeStampedModel

from snpdb.models import Lab


class UploadedFileLabStatus(models.TextChoices):
    Pending = 'P', 'Pending'
    AutoProcessed = 'AP', 'Automatically-Processed'
    Processed = 'MP', 'Processed'
    Error = 'E', 'Error'


class UploadedFileLab(TimeStampedModel):
    class Meta:
        verbose_name = "Lab Classification File"

    url = models.TextField()
    filename = models.TextField()
    lab = models.ForeignKey(Lab, on_delete=models.CASCADE)
    user = models.ForeignKey(User, on_delete=models.PROTECT)
    comment = models.TextField(default="")
    status = models.CharField(max_length=2, choices=UploadedFileLabStatus.choices, default=UploadedFileLabStatus.Pending)
