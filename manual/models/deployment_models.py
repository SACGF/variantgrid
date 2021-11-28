from django.db import models
from model_utils.models import TimeStampedModel


class Deployment(TimeStampedModel):
    git_hash = models.TextField(null=True, blank=True)
    release_note = models.TextField(null=True, blank=True)
    app = models.TextField(null=True, blank=True)
