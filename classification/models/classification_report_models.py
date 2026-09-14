"""
Report templates and the case reports built from them.

ClassificationReportTemplate is a lab maintained document design: `template` is the single record
HTML/Vue report, `case_template` / `json_template` the multi-variant case report, which is rendered
server side so the HTML preview, the PDF and the Word file cannot disagree.

CaseReport is one run of a case_template over a case's classifications, pinned to the
ClassificationModifications it rendered, so a later edit changes nothing about a report already
issued. @see classification/report/case_report_context.py for what the templates are handed.
"""
import os
from typing import Optional

from django.conf import settings
from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied, ValidationError
from django.db import models
from django.db.models import CASCADE, PROTECT, SET_NULL, Prefetch, Q, QuerySet
from model_utils.models import TimeStampedModel

from classification.enums import AlleleOriginBucket
from classification.models.classification import ClassificationModification
from classification.report.template_validation import validate_case_template, validate_json_template
from patients.models import Extraction, Patient, Specimen
from patients.models_enums import SampleSourceLevel
from snpdb.models import Lab, Sample


class ReportNames:
    DEFAULT_REPORT = "default_report"
    DEFAULT_CASE_REPORT = "default_case_report"


class ClassificationReportTemplate(TimeStampedModel):
    name = models.TextField(primary_key=True)
    template = models.TextField(null=False, blank=True, default="")
    # The case report - one Django template rendered to HTML, with the PDF and DOCX derived from it
    case_template = models.TextField(null=False, blank=True, default="")
    # Blank renders the canonical context dump rather than nothing
    json_template = models.TextField(null=False, blank=True, default="")
    # Case level inputs the build form asks for, so a deployment adds them without a schema change:
    # [{"key", "label", "type": "text"|"bool"|"choice", "options": [...], "default", "group",
    #   "prefill_key": the evidence key the form starts the field from}]
    case_fields = models.JSONField(default=list, blank=True)
    # Which cases this template is offered for - null is every case
    allele_origin_bucket = models.CharField(max_length=1, choices=AlleleOriginBucket.choices,
                                            null=True, blank=True)

    @property
    def has_case_report(self) -> bool:
        """ A template can build a case report once it has the document design for one """
        return bool(self.case_template)

    @property
    def case_field_groups(self) -> list[dict]:
        """ case_fields laid out for the build form - the fields sharing a `group` are one row of
            checkboxes, and the context exposes that group as a dict """
        groups = []
        by_name = {}
        for field in self.case_fields or []:
            name = field.get("group") or ""
            if name not in by_name:
                by_name[name] = {"name": name, "fields": []}
                groups.append(by_name[name])
            by_name[name]["fields"].append(field)
        return groups

    def case_values_from_form(self, posted) -> dict:
        """ The build form's answers to case_fields, read back by key. A bool field is the presence
            of its checkbox, and the fields sharing a `group` come back as one dict - which is the
            shape a JSON template reads (`case_values.assay_success`) """
        values = {}
        groups: dict[str, dict] = {}
        for field in self.case_fields or []:
            if not (key := field.get("key")):
                continue
            if field.get("type") == "bool":
                value = f"case_field_{key}" in posted
            else:
                value = posted.get(f"case_field_{key}") or field.get("default") or ""
            if group := field.get("group"):
                groups.setdefault(group, {})[key] = value
            else:
                values[key] = value
        values.update(groups)
        return values

    def clean(self):
        super().clean()
        errors = {}
        if self.case_template:
            if message := validate_case_template(self.case_template):
                errors["case_template"] = message
        if self.json_template:
            if message := validate_json_template(self.json_template):
                errors["json_template"] = message
        if errors:
            raise ValidationError(errors)

    @staticmethod
    def preferred_template_for(c: ClassificationModification) -> Optional['ClassificationReportTemplate']:
        """
        :param c: The classification that we're generating a report template for, in future the allele origin bucket of hte
        classification will most likely help define the template
        :return: A report template
        """
        return ClassificationReportTemplate.objects.filter(name=ReportNames.DEFAULT_REPORT).first()

    @staticmethod
    def case_templates_for_bucket(allele_origin_bucket: Optional[str] = None):
        """ The templates a case can be reported with - a template with no bucket suits every case """
        qs = ClassificationReportTemplate.objects.exclude(case_template="")
        if allele_origin_bucket:
            qs = qs.filter(Q(allele_origin_bucket__isnull=True) |
                           Q(allele_origin_bucket=allele_origin_bucket))
        return qs.order_by("name")

    def __str__(self):
        return self.name


class CaseReportStatus(models.TextChoices):
    DRAFT = 'D', 'Draft'
    FINAL = 'F', 'Final'
    SUPERSEDED = 'S', 'Superseded'


def case_report_upload_path(instance: 'CaseReport', filename: str) -> str:
    return os.path.join("case_reports", str(instance.lab_id), str(instance.pk), filename)


# Which FK each level of the hierarchy uses - the level says which, so one table covers creating a
# report, finding a case's reports and resolving the one it is about
SOURCE_LEVEL_FIELDS = {
    SampleSourceLevel.PATIENT: "patient",
    SampleSourceLevel.SPECIMEN: "specimen",
    SampleSourceLevel.EXTRACTION: "extraction",
    SampleSourceLevel.SAMPLE: "sample",
}


class CaseReport(TimeStampedModel):
    """ One rendering of a case_template over a case's classifications.

        The case is a level of the Patient -> Specimen -> Extraction -> Sample hierarchy, spelt the
        way SampleNode spells it - source_level plus the one FK that level uses - so a TSO 500 case
        with a DNA and an RNA arm is reported from the specimen the measures hang off.

        Everything the templates saw is kept in context_snapshot, so the documents can be re-rendered
        after a template fix and the report's numbers stay answerable even though a SpecimenMeasure
        resend replaces the row they came from. """
    template = models.ForeignKey(ClassificationReportTemplate, on_delete=PROTECT)
    lab = models.ForeignKey(Lab, on_delete=CASCADE)
    user = models.ForeignKey(User, on_delete=PROTECT)
    status = models.CharField(max_length=1, choices=CaseReportStatus.choices, default=CaseReportStatus.DRAFT)
    supersedes = models.ForeignKey('self', null=True, blank=True, on_delete=SET_NULL)

    source_level = models.CharField(max_length=1, choices=SampleSourceLevel.choices)
    patient = models.ForeignKey(Patient, null=True, blank=True, on_delete=SET_NULL)
    specimen = models.ForeignKey(Specimen, null=True, blank=True, on_delete=SET_NULL)
    extraction = models.ForeignKey(Extraction, null=True, blank=True, on_delete=SET_NULL)
    sample = models.ForeignKey(Sample, null=True, blank=True, on_delete=SET_NULL)

    summary = models.TextField(blank=True)  # the case level "Summary Interpretation"
    case_values = models.JSONField(default=dict, blank=True)  # answers to template.case_fields
    external_report_id = models.TextField(null=True, blank=True)  # the LIS's report ID, once issued
    report_date = models.DateField(null=True, blank=True)

    context_snapshot = models.JSONField(default=dict)
    html = models.TextField(blank=True)  # the rendered document - PDF and DOCX are derived from it
    pdf_file = models.FileField(null=True, blank=True, upload_to=case_report_upload_path)
    docx_file = models.FileField(null=True, blank=True, upload_to=case_report_upload_path)
    json_output = models.JSONField(null=True, blank=True)

    class Meta:
        indexes = [
            models.Index(fields=["specimen"]),
            models.Index(fields=["patient"]),
            models.Index(fields=["sample"]),
        ]
        constraints = [
            # Exactly one source - the level says which, and a second FK would make the case ambiguous
            models.CheckConstraint(
                name="case_report_one_source",
                condition=(Q(patient__isnull=False, specimen__isnull=True, extraction__isnull=True, sample__isnull=True) |
                           Q(patient__isnull=True, specimen__isnull=False, extraction__isnull=True, sample__isnull=True) |
                           Q(patient__isnull=True, specimen__isnull=True, extraction__isnull=False, sample__isnull=True) |
                           Q(patient__isnull=True, specimen__isnull=True, extraction__isnull=True, sample__isnull=False))),
        ]

    @property
    def source(self):
        """ The object the case is - whichever level source_level names """
        return getattr(self, SOURCE_LEVEL_FIELDS[self.source_level])

    @property
    def is_editable(self) -> bool:
        """ A FINAL report is the copy that went out with the case - a change makes a new version """
        return self.status == CaseReportStatus.DRAFT

    @staticmethod
    def source_kwargs(source_level: str, source) -> dict:
        """ The level and the one FK it uses, for creating a report or finding a case's """
        return {"source_level": source_level, SOURCE_LEVEL_FIELDS[source_level]: source}

    @staticmethod
    def for_case(source_level: str, source) -> QuerySet['CaseReport']:
        """ A case's reports, newest first - whoever in the lab built them. The Reports card lists
            every report's pinned classifications, so they come along rather than one query a row """
        pinned = CaseReportClassification.objects.select_related(
            "classification_modification__classification__lab")
        return CaseReport.objects.filter(**CaseReport.source_kwargs(source_level, source)) \
            .select_related("template", "user", "lab") \
            .prefetch_related(Prefetch("casereportclassification_set", queryset=pinned)) \
            .order_by("-pk")

    def can_write(self, user: User) -> bool:
        """ Building, finalising and deleting a report is the owning lab's - it is their document """
        if user.is_superuser:
            return True
        return Lab.valid_labs_qs(user).filter(pk=self.lab_id).exists()

    def can_view(self, user: User) -> bool:
        """ The lab's own users, plus anyone who could already see the case it is about """
        if self.can_write(user):
            return True
        if source := self.source:
            return source.can_view(user)
        return False

    def check_can_view(self, user: User):
        if not self.can_view(user):
            raise PermissionDenied(f"You do not have permission to view case report {self.pk}")

    def check_can_write(self, user: User):
        if not self.can_write(user):
            raise PermissionDenied(f"Only {self.lab}'s users can change case report {self.pk}")

    @property
    def rows(self) -> QuerySet['CaseReportClassification']:
        """ The pinned classifications in report order - plain, so for_case's prefetch is used """
        return self.casereportclassification_set.all()

    def get_media_dir(self) -> str:
        return os.path.join(settings.MEDIA_ROOT, "case_reports", str(self.lab_id), str(self.pk))

    def __str__(self):
        return f"{self.get_status_display()} report for {self.source} ({self.template_id})"


class CaseReportClassification(models.Model):
    """ One classification in a report, pinned to the version that was rendered """
    case_report = models.ForeignKey(CaseReport, on_delete=CASCADE)
    classification_modification = models.ForeignKey(ClassificationModification, on_delete=PROTECT)
    order = models.IntegerField()  # position in the report, as computed at build time
    # "Report": Y/N - in the document, vs listed as seen only
    reported = models.BooleanField(default=True)

    class Meta:
        unique_together = ("case_report", "classification_modification")
        ordering = ("order",)

    def __str__(self):
        return f"{self.case_report_id}#{self.order} {self.classification_modification}"
