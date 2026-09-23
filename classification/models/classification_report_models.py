"""
Report templates and the case reports built from them.

ClassificationReportTemplate is a lab maintained document design: `template` is the single record
HTML/Vue report and `case_template` the multi-variant case report, which is rendered server side so
the HTML preview, the PDF and the Word file cannot disagree. The report's JSON is not a template -
@see classification/report/renderers.py:render_json.

CaseReport is one run of a case_template over a case's classifications, pinned to the
ClassificationModifications it rendered, so a later edit changes nothing about a report already
issued. @see classification/report/case_report_context.py for what the templates are handed.
"""
import logging
import os
from dataclasses import dataclass
from typing import Optional

import django.dispatch
from django.conf import settings
from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied, ValidationError
from django.db import models
from django.db.models import CASCADE, PROTECT, SET_NULL, Prefetch, Q, QuerySet
from model_utils.models import TimeStampedModel

from classification.enums import AlleleOriginBucket
from classification.models.classification import ClassificationModification
from classification.report.template_validation import validate_case_template
from library.case_report_delivery import CaseReportDelivery
from patients.models import Extraction, Patient, Specimen, measure_value_description
from patients.models_enums import (
    MEASURE_CONTEXT_KEYS,
    SampleSourceLevel,
)
from seqauto.models import LibraryQC
from seqauto.models.models_enums import CVO_MEASURE_CONTEXT_KEYS, LIBRARY_QC_CONTEXT_KEYS
from snpdb.models import Lab, Sample

# A report has gone out with the case, so a deployment that files it somewhere else can now do so.
# Only the DRAFT -> FINAL step sends it - re-finalising is how a LIS id entered later reaches the records
case_report_finalised_signal = django.dispatch.Signal()  # args: "case_report", "user"
# What each of those deliveries came to, for the Reports card - a receiver answers with
# dict[case report pk, list[CaseReportDelivery]]
case_report_deliveries_signal = django.dispatch.Signal()  # args: "case_reports"
# A report's JSON, built in Python by whichever app owns that template's shape. The JSON is an
# interface to another system rather than a document, so the app that has to keep it in step with
# that system writes it - see classification/report/renderers.py:render_json
case_report_json_signal = django.dispatch.Signal()  # args: "report_template", "context"; returns dict


def get_case_report_deliveries(case_reports: list['CaseReport']) -> dict[int, list[CaseReportDelivery]]:
    """ The deliveries every app has for these reports, keyed by report pk. send_robust, the way
        library/integration_status.py collects its statuses - an app that can't answer costs the
        card its column, not the tab """
    deliveries: dict[int, list[CaseReportDelivery]] = {}
    if not case_reports:
        return deliveries
    for caller, result in case_report_deliveries_signal.send_robust(sender=CaseReport, case_reports=case_reports):
        if isinstance(result, Exception):
            logging.error("Exception getting case report deliveries from %s: %s", caller, result)
        elif result:
            for case_report_id, report_deliveries in result.items():
                deliveries.setdefault(case_report_id, []).extend(report_deliveries)
    return deliveries


@dataclass(frozen=True)
class Measure:
    """ One number a case report prints and a build form ticks from, whatever it came off - a sequencing
        analysis' TMB or the pathologist's tumour content (@see case_report_context.case_measures) """
    value: Optional[float]
    unit: Optional[str]
    call: Optional[str]
    threshold: Optional[str]         # the policy in words, describe_bands(...)
    threshold_source: Optional[str]  # whose policy - the settings that set it
    method: str

    @property
    def value_description(self) -> str:
        return measure_value_description(self.value, self.unit, self.call)


# Every key a case_field's `measure` can name
MEASURE_KEYS = frozenset(CVO_MEASURE_CONTEXT_KEYS) | frozenset(MEASURE_CONTEXT_KEYS.values())


# What a bool case_field's `tick_when` can say about the measure it names. The form starts the tick
# from the rule and the scientist adjusts it - the answer the report goes out with is still theirs.
# A list of rules ticks when any one holds (a purity caveat off the pathologist's call or the number)
TICK_WHEN_CALLED = "called"
TICK_WHEN_CALL_IN = "call_in"
TICK_WHEN_VALUE_BELOW = "value_below"
TICK_WHEN_RULES = (TICK_WHEN_CALLED, TICK_WHEN_CALL_IN, TICK_WHEN_VALUE_BELOW)

# What a bool case_field's `tick_when` can say about the library QC category it names instead -
# 'the assay succeeded for amplifications' is DRAGEN's CNV library QC, not a measure
TICK_WHEN_PASSED = "passed"
TICK_WHEN_COMPLETED = "completed"
TICK_WHEN_QC_RULES = (TICK_WHEN_PASSED, TICK_WHEN_COMPLETED)


def _tick_when_rules(tick_when) -> list[dict]:
    return tick_when if isinstance(tick_when, list) else [tick_when]


def _rule_holds(rule: dict, measure: Measure) -> Optional[bool]:
    if (called := rule.get(TICK_WHEN_CALLED)) is not None:
        return bool(measure.call) == bool(called)
    if (call_in := rule.get(TICK_WHEN_CALL_IN)) is not None:
        return measure.call in call_in
    if (value_below := rule.get(TICK_WHEN_VALUE_BELOW)) is not None:
        return measure.value is not None and measure.value < value_below
    return None


def measure_tick(tick_when, measure: Optional[Measure]) -> Optional[bool]:
    """ Whether the rule a case_field states holds for the case's measure - None where the case has
        no such measure, so the field's own default stands and the form says there is none """
    if measure is None:
        return None
    outcomes = [_rule_holds(rule, measure) for rule in _tick_when_rules(tick_when)]
    if all(outcome is None for outcome in outcomes):
        return None
    return any(outcomes)


def _qc_rule_holds(rule: dict, library_qc: LibraryQC) -> Optional[bool]:
    """ None where the row cannot answer the rule - a category with no QC leaves the field's own
        default standing rather than reading as a failure """
    if (passed := rule.get(TICK_WHEN_PASSED)) is not None:
        return None if library_qc.passed is None else library_qc.passed == bool(passed)
    if (completed := rule.get(TICK_WHEN_COMPLETED)) is not None:
        return None if library_qc.completed is None else library_qc.completed == bool(completed)
    return None


def library_qc_tick(tick_when, library_qc: Optional[LibraryQC]) -> Optional[bool]:
    """ Whether the rule a case_field states holds for the case's library QC - None where the case
        has no QC for that category, so the field's own default stands """
    if library_qc is None:
        return None
    outcomes = [_qc_rule_holds(rule, library_qc) for rule in _tick_when_rules(tick_when)]
    if all(outcome is None for outcome in outcomes):
        return None
    return any(outcomes)


def tick_for(field: dict, measures: dict, library_qc: dict) -> Optional[bool]:
    """ Where a bool case_field's tick starts - a field names either a measure or a QC category,
        and the one it names says which of the case's rows the rule is applied to """
    tick_when = field.get("tick_when")
    if not tick_when:
        return None
    if qc_key := field.get("qc"):
        return library_qc_tick(tick_when, library_qc.get(qc_key))
    return measure_tick(tick_when, measures.get(field.get("measure")))


def _describe_rule(rule: dict, unit: str) -> str:
    if (called := rule.get(TICK_WHEN_CALLED)) is not None:
        return "the measure has a call" if called else "the measure has no call"
    if (call_in := rule.get(TICK_WHEN_CALL_IN)) is not None:
        return "the call is " + " or ".join(str(call) for call in call_in)
    if (value_below := rule.get(TICK_WHEN_VALUE_BELOW)) is not None:
        return f"the value is below {value_below}{unit}"
    if (passed := rule.get(TICK_WHEN_PASSED)) is not None:
        return "the library passed QC" if passed else "the library failed QC"
    if (completed := rule.get(TICK_WHEN_COMPLETED)) is not None:
        return "the run completed for the library" if completed else "the run did not complete for the library"
    return ""


def describe_tick_when(tick_when, unit: Optional[str] = None) -> str:
    """ The rule in words, shown beside the checkbox so a scientist who disagrees with a starting
        tick can see the policy behind it and ask for it to change - 'ticked when the value is below 20%' """
    unit = unit or ""
    if unit not in ("", "%"):
        unit = f" {unit}"
    described = [text for text in (_describe_rule(rule, unit) for rule in _tick_when_rules(tick_when)) if text]
    return "ticked when " + " or ".join(described) if described else ""


def validate_case_fields(case_fields: list) -> Optional[str]:
    """ What a template's JSON has to get right for a measure to reach the form - the keys are hand
        written in admin, so a typo says so at save time rather than silently ticking nothing """
    measure_keys = MEASURE_KEYS
    qc_keys = set(LIBRARY_QC_CONTEXT_KEYS.values())
    for field in case_fields or []:
        if not isinstance(field, dict):
            return f"Each case field must be an object, not '{field}'"
        key = field.get("key") or "(no key)"
        measure = field.get("measure")
        qc = field.get("qc")
        if measure is not None and measure not in measure_keys:
            return f"'{key}' measures '{measure}' - one of {', '.join(sorted(measure_keys))} was expected"
        if qc is not None and qc not in qc_keys:
            return f"'{key}' names library QC '{qc}' - one of {', '.join(sorted(qc_keys))} was expected"
        if measure is not None and qc is not None:
            # The rules are different and the row they judge is different, so a field is one or the other
            return f"'{key}' names both a measure and a library QC category - it can carry one"
        rules = TICK_WHEN_QC_RULES if qc else TICK_WHEN_RULES
        if (tick_when := field.get("tick_when")) is not None:
            for rule in _tick_when_rules(tick_when):
                if not isinstance(rule, dict) or not rule:
                    return f"'{key}' tick_when must be one of {', '.join(rules)}"
                if unknown := set(rule) - set(rules):
                    return f"'{key}' tick_when has no rule '{', '.join(sorted(unknown))}'"
            if measure is None and qc is None:
                return f"'{key}' has a tick_when and names no measure or library QC to apply it to"
    return None


class ReportNames:
    DEFAULT_REPORT = "default_report"
    DEFAULT_CASE_REPORT = "default_case_report"


class ClassificationReportTemplate(TimeStampedModel):
    name = models.TextField(primary_key=True)
    template = models.TextField(null=False, blank=True, default="")
    # The case report - one Django template rendered to HTML, with the PDF and DOCX derived from it
    case_template = models.TextField(null=False, blank=True, default="")
    # Case level inputs the build form asks for, so a deployment adds them without a schema change:
    # [{"key", "label", "type": "text"|"bool"|"choice", "options": [...], "default", "group",
    #   "prefill_key": the evidence key the form starts the field from,
    #   "measure": the Measure shown beside a bool field (a MEASURE_KEYS value),
    #   "qc": the LibraryQC category shown beside a bool field instead (a LIBRARY_QC_CONTEXT_KEYS value),
    #   "tick_when": the rule (or list of rules, any of which) that field's tick starts from - @see tick_for}]
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
        if message := validate_case_fields(self.case_fields):
            errors["case_fields"] = message
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
        after a template fix and the report's numbers stay answerable even though the pathologist's
        tumour content can be re-sent over the row it came from. """
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
