from typing import Optional

from model_utils.models import TimeStampedModel
from django.db import models

from classification.models import ClassificationModification


class ReportNames:
    DEFAULT_REPORT = "default_report"


class ClassificationReportTemplate(TimeStampedModel):
    name = models.TextField(primary_key=True)
    template = models.TextField(null=False, blank=True, default="")

    @staticmethod
    def preferred_template_for(c: ClassificationModification) -> Optional['ClassificationReportTemplate']:
        return ClassificationReportTemplate.objects.filter(name=ReportNames.DEFAULT_REPORT).first()
