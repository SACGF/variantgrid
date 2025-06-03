from typing import Optional

from django.db import models
from model_utils.models import TimeStampedModel

from classification.models import ClassificationModification


class ReportNames:
    DEFAULT_REPORT = "default_report"


class ClassificationReportTemplate(TimeStampedModel):
    name = models.TextField(primary_key=True)
    template = models.TextField(null=False, blank=True, default="")

    @staticmethod
    def preferred_template_for(c: ClassificationModification) -> Optional['ClassificationReportTemplate']:
        """
        :param c: The classification that we're generating a report template for, in future the allele origin bucket of the
        classification will most likely help define the template
        :return: A report template
        """
        return ClassificationReportTemplate.objects.filter(name=ReportNames.DEFAULT_REPORT).first()
