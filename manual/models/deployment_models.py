from model_utils.models import TimeStampedModel
from django.db import models


class Deployment(TimeStampedModel):
    git_hash = models.TextField(null=True, blank=True)
    release_note = models.TextField(null=True, blank=True)
    app = models.TextField(null=True, blank=True)