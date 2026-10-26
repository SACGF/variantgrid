"""
IntegrationActivityStatus: the outcome of one integration attempt passed to IntegrationActivity.record()
(success, no change, error).
"""
from django.db import models


class IntegrationActivityStatus(models.TextChoices):
    SUCCESS = "S", "Success"
    NO_CHANGE = "N", "No change"
    ERROR = "E", "Error"
