""" The case report half of the Classify & Report tab (#444).

    Ordering, amp_tier and the renderers are covered over fake modifications in
    classification/tests/report/; what needs a real case is here: which samples a specimen or
    extraction case is, who may act on a built report, and what finalising writes.
"""
from datetime import timedelta

from django.contrib.auth.models import User
from django.core.exceptions import ValidationError
from django.test import override_settings
from django.urls import reverse
from django.utils import timezone

from analysis.classify_report import ClassifyReportCase, ReportCandidate
from analysis.tests.test_classify_report import READY_EVIDENCE, ClassifyReportTestCase
from analysis.views.views_classify_report import _case_values_for_form
from classification.enums import SpecialEKeys, SubmissionSource
from classification.models import (
    CaseReport,
    CaseReportStatus,
    Classification,
    ClassificationModification,
    ClassificationReportTemplate,
)
from classification.models.classification_report_models import (
    case_report_deliveries_signal,
    case_report_finalised_signal,
)
from classification.report.case_report_builder import build_case_report, finalise_case_report
from classification.report.case_report_context import specimen_library_qc
from classification.report.default_templates import generic_case_template
from library.case_report_delivery import CaseReportDelivery
from library.guardian_utils import assign_permission_to_user_and_groups
from patients.models import Extraction, Patient, Specimen, SpecimenMeasure
from patients.models_enums import (
    NucleicAcid,
    SampleSourceLevel,
    SpecimenMeasureType,
)
from patients.sample_grouping import get_sample_group
from seqauto.models import LibraryQC
from seqauto.models.models_enums import LibraryQCCategory
from snpdb.models import Sample


@override_settings(CELERY_TASK_ALWAYS_EAGER=True, LIFTOVER_CLASSIFICATIONS=False,
                   CLINGEN_ALLELE_REGISTRY_LOGIN=None)
class SpecimenAndExtractionCaseTest(ClassifyReportTestCase):
    """ A specimen case is its extractions' samples - the same resolution the analysis grouping node
        does, so the tab and the node never disagree about which of a case's arms are in play """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.patient = Patient.objects.create(first_name="Case", last_name="Report")
        assign_permission_to_user_and_groups(cls.user, cls.patient)
        cls.specimen = Specimen.objects.create(reference_id="2600000003", patient=cls.patient)
        cls.dna = Extraction.objects.create(specimen=cls.specimen, reference_id="2600000003C",
                                            nucleic_acid_source=NucleicAcid.DNA)
        cls.rna = Extraction.objects.create(specimen=cls.specimen, reference_id="2600000003R",
                                            nucleic_acid_source=NucleicAcid.RNA)
        Sample.objects.filter(pk=cls.proband.pk).update(extraction=cls.dna)
        Sample.objects.filter(pk=cls.mother.pk).update(extraction=cls.rna)

    def test_a_specimen_case_is_its_extractions_samples(self):
        case = ClassifyReportCase.for_specimen(self.user, self.specimen)

        expected = get_sample_group(self.user, SampleSourceLevel.SPECIMEN, self.specimen).samples
        self.assertEqual(case.samples, expected)
        self.assertEqual(case.source_level, SampleSourceLevel.SPECIMEN)
        self.assertEqual(case.patients, [self.patient])

    def test_an_extraction_case_is_that_arms_samples_alone(self):
        case = ClassifyReportCase.for_extraction(self.user, self.dna)

        expected = get_sample_group(self.user, SampleSourceLevel.EXTRACTION, self.dna).samples
        self.assertEqual(case.samples, expected)
        self.assertEqual([s.pk for s in case.samples], [self.proband.pk])


@override_settings(CELERY_TASK_ALWAYS_EAGER=True, LIFTOVER_CLASSIFICATIONS=False,
                   CLINGEN_ALLELE_REGISTRY_LOGIN=None)
class CaseReportFinaliseTest(ClassifyReportTestCase):
    """ Finalising stamps the report onto its classifications and re-publishes them """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.template = ClassificationReportTemplate.objects.create(
            name="case report test template", case_template=generic_case_template())

    def _build(self) -> tuple[CaseReport, Classification]:
        classification = self._classify(self.proband, data={SpecialEKeys.GENE_SYMBOL: {"value": "RUNX1"}})
        classification.publish_latest(self.user)
        case_report = build_case_report(self.user, self.template, self.lab, SampleSourceLevel.SAMPLE,
                                        self.proband, [classification.last_published_version])
        return case_report, classification

    def test_finalise_stamps_the_report_and_repins_the_published_version(self):
        case_report, classification = self._build()

        result = finalise_case_report(case_report, self.user)

        case_report.refresh_from_db()
        self.assertEqual(case_report.status, CaseReportStatus.FINAL)
        self.assertIsNotNone(case_report.report_date)
        self.assertEqual(result.stamped, [classification])
        classification.refresh_from_db()
        self.assertEqual(classification.get(SpecialEKeys.REPORT_DATE),
                         case_report.report_date.isoformat())
        self.assertEqual(classification.get(SpecialEKeys.VARIANT_REPORTED), "primary_finding")
        # The pinned version is the one that carries the stamp, so the tab's stale check stays quiet
        self.assertEqual(case_report.rows.get().classification_modification,
                         classification.last_published_version)

    def test_finalise_is_idempotent(self):
        case_report, classification = self._build()
        finalise_case_report(case_report, self.user)
        versions = ClassificationModification.objects.filter(classification=classification).count()

        result = finalise_case_report(case_report, self.user)

        self.assertEqual(result.stamped, [])
        self.assertEqual(ClassificationModification.objects.filter(
            classification=classification).count(), versions)

    def test_a_record_with_unsubmitted_edits_is_left_alone(self):
        """ Publishing it would push out edits nobody has submitted, so it is named instead """
        case_report, classification = self._build()
        classification.patch_value({SpecialEKeys.INTERPRETATION_SUMMARY: {"value": "still working"}},
                                   user=self.user, source=SubmissionSource.FORM, save=True)
        versions = ClassificationModification.objects.filter(classification=classification).count()

        result = finalise_case_report(case_report, self.user)

        self.assertEqual(result.unsubmitted, [classification])
        self.assertEqual(result.stamped, [])
        self.assertEqual(ClassificationModification.objects.filter(
            classification=classification).count(), versions)
        classification.refresh_from_db()
        self.assertIsNone(classification.get(SpecialEKeys.REPORT_DATE))


@override_settings(CELERY_TASK_ALWAYS_EAGER=True, LIFTOVER_CLASSIFICATIONS=False,
                   CLINGEN_ALLELE_REGISTRY_LOGIN=None)
class CaseReportDeploymentHooksTest(ClassifyReportTestCase):
    """ The two hooks a deployment specific app files a report through - SA Path sends finalised
        TSO 500 reports to Mocha off these """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.template = ClassificationReportTemplate.objects.create(
            name="case report hooks template", case_template=generic_case_template())
        classification = Classification.create(user=cls.user, lab=cls.lab, sample=cls.proband,
                                               source=SubmissionSource.VARIANT_GRID, variant=cls.variant)
        classification.publish_latest(cls.user)
        cls.case_report = build_case_report(cls.user, cls.template, cls.lab, SampleSourceLevel.SAMPLE,
                                            cls.proband, [classification.last_published_version])

    def test_finalised_signal_fires_only_on_the_draft_to_final_step(self):
        """ Re-finalising is how a LIS id entered later reaches the records - it sends nothing out again """
        finalised = []

        def receiver(sender, case_report, user, **kwargs):
            finalised.append((case_report, user))

        case_report_finalised_signal.connect(receiver)
        try:
            finalise_case_report(self.case_report, self.user)
            self.assertEqual(finalised, [(self.case_report, self.user)])
            finalise_case_report(self.case_report, self.user)
            self.assertEqual(len(finalised), 1)
        finally:
            case_report_finalised_signal.disconnect(receiver)

    def test_the_card_shows_what_a_delivery_receiver_returns(self):
        def receiver(sender, case_reports, **kwargs):
            return {cr.pk: [CaseReportDelivery(label="Mocha", status="warning", text="Sent, pending",
                                               action_url="/sapath/mocha/send", action_label="Send to Mocha")]
                    for cr in case_reports}

        case_report_deliveries_signal.connect(receiver)
        self.client.force_login(self.user)
        try:
            response = self.client.get(reverse("sample_classify_report_tab",
                                               kwargs={"sample_id": self.proband.pk}))
        finally:
            case_report_deliveries_signal.disconnect(receiver)

        self.assertContains(response, "Sent, pending")
        self.assertContains(response, "Send to Mocha")


@override_settings(CELERY_TASK_ALWAYS_EAGER=True, LIFTOVER_CLASSIFICATIONS=False,
                   CLINGEN_ALLELE_REGISTRY_LOGIN=None)
class CaseReportPermissionTest(ClassifyReportTestCase):
    """ A report is its lab's document, and the case's patient's - neither is open to a passer-by """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.template = ClassificationReportTemplate.objects.create(
            name="case report permission template", case_template=generic_case_template())
        classification = Classification.create(user=cls.user, lab=cls.lab, sample=cls.proband,
                                                source=SubmissionSource.VARIANT_GRID,
                                                variant=cls.variant)
        classification.publish_latest(cls.user)
        cls.case_report = build_case_report(cls.user, cls.template, cls.lab,
                                            SampleSourceLevel.SAMPLE, cls.proband,
                                            [classification.last_published_version])
        cls.outsider = User.objects.create_user(username="case_report_outsider")

    def _urls(self) -> dict[str, str]:
        return {
            "download": reverse("case_report_download",
                                kwargs={"case_report_id": self.case_report.pk, "document_format": "pdf"}),
            "finalise": reverse("case_report_finalise", kwargs={"case_report_id": self.case_report.pk}),
        }

    def test_a_user_outside_the_lab_cannot_download_or_finalise(self):
        self.client.force_login(self.outsider)
        urls = self._urls()

        self.assertEqual(self.client.get(urls["download"]).status_code, 403)
        self.assertEqual(self.client.post(urls["finalise"]).status_code, 403)

    def test_the_labs_own_user_can_download_and_finalise(self):
        self.client.force_login(self.user)
        urls = self._urls()

        self.assertEqual(self.client.get(urls["download"]).status_code, 200)
        self.assertEqual(self.client.post(urls["finalise"]).status_code, 200)

    def test_the_preview_can_be_framed_by_the_built_report_modal(self):
        """ case_report_built.html shows it in an iframe, which the site-wide DENY blanks """
        self.client.force_login(self.user)
        response = self.client.get(reverse("view_case_report", kwargs={"case_report_id": self.case_report.pk}))

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response["X-Frame-Options"], "SAMEORIGIN")
        # HtmlMinifyMiddleware turns &nbsp; into the character, which reads as "Â " without a charset
        self.assertIn("charset=utf-8", response["Content-Type"])


@override_settings(CELERY_TASK_ALWAYS_EAGER=True, LIFTOVER_CLASSIFICATIONS=False,
                   CLINGEN_ALLELE_REGISTRY_LOGIN=None)
class CaseReportBuildFormTest(ClassifyReportTestCase):
    """ The Build report form: the ticked classifications, their Report Y/N, and the template's own
        case_fields - which are how a deployment adds report level inputs without a schema change """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.template = ClassificationReportTemplate.objects.create(
            name="case report form template", case_template=generic_case_template(),
            case_fields=[{"key": "panel", "label": "Panel", "type": "text"},
                         {"key": "tmb", "label": "TMB", "type": "bool", "group": "assay_success"},
                         {"key": "msi", "label": "MSI", "type": "bool", "group": "assay_success"}])
        cls.reported = Classification.create(user=cls.user, lab=cls.lab, sample=cls.proband,
                                              source=SubmissionSource.VARIANT_GRID, variant=cls.variant,
                                              data={**READY_EVIDENCE, SpecialEKeys.GENE_SYMBOL: {"value": "RUNX1"}})
        cls.reported.publish_latest(cls.user)
        cls.unreported = Classification.create(user=cls.user, lab=cls.lab, sample=cls.proband,
                                                source=SubmissionSource.VARIANT_GRID,
                                                variant=cls.shared_variant,
                                                data={**READY_EVIDENCE, SpecialEKeys.GENE_SYMBOL: {"value": "TP53"}})
        cls.unreported.publish_latest(cls.user)

    def setUp(self):
        super().setUp()
        self.client.force_login(self.user)

    def _post_data(self) -> dict:
        reported_pk = self.reported.last_published_version.pk
        return {
            "report_template": self.template.pk,
            "lab": self.lab.pk,
            "classification_modification_id": [reported_pk, self.unreported.last_published_version.pk],
            f"reported_{reported_pk}": "on",  # the other is left off the document
            "case_field_panel": "TSO500",
            "case_field_tmb": "on",
            "summary": "Summary interpretation.",
        }

    def _url(self, name: str) -> str:
        return reverse(name, kwargs={"case_type": "sample", "case_id": self.proband.pk})

    def test_the_dialog_lists_the_ticked_classifications_and_the_templates_fields(self):
        response = self.client.post(self._url("case_report_build_dialog"), self._post_data())

        self.assertContains(response, "RUNX1")
        self.assertContains(response, "TP53")
        self.assertContains(response, "case_field_panel")
        self.assertContains(response, "assay_success")

    def test_a_record_with_errors_cannot_be_put_in_a_report(self):
        """ Creating a record publishes it to the lab with its mandatory keys still empty, so being
            published is not being ready - its tick box is disabled, and posting it anyway drops it """
        incomplete = Classification.create(user=self.user, lab=self.lab, sample=self.mother,
                                           source=SubmissionSource.VARIANT_GRID, variant=self.variant,
                                           data={SpecialEKeys.GENE_SYMBOL: {"value": "BRCA2"}})
        incomplete.publish_latest(self.user)
        candidate = ReportCandidate(incomplete.last_published_version)
        self.assertFalse(candidate.ready)
        self.assertGreater(candidate.error_count, 0)

        url = reverse("case_report_build_dialog", kwargs={"case_type": "sample", "case_id": self.mother.pk})
        response = self.client.post(url, {"classification_modification_id": [candidate.modification.pk]})
        self.assertEqual(response.status_code, 404)

    def test_a_record_with_unsubmitted_edits_is_not_ready(self):
        candidate = ReportCandidate(self.reported.last_published_version)
        self.assertTrue(candidate.ready)

        self.reported.patch_value({SpecialEKeys.INTERPRETATION_SUMMARY: {"value": "still working"}},
                                  user=self.user, source=SubmissionSource.FORM, save=True)

        candidate = ReportCandidate(self.reported.last_published_version)
        self.assertTrue(candidate.unsubmitted)
        self.assertFalse(candidate.ready)

    def test_the_dialog_starts_a_field_from_its_default_and_its_prefill_key(self):
        """ A case level value the records already carry is not worth retyping, and a flag that is
            normally on starts on - both only reach the form through the field definition """
        self.template.case_fields = [
            {"key": "panel", "label": "Panel", "type": "text", "default": "TSO500"},
            {"key": "gene", "label": "Gene", "type": "text", "prefill_key": SpecialEKeys.GENE_SYMBOL},
        ]
        self.template.save()

        response = self.client.post(self._url("case_report_build_dialog"), self._post_data())

        content = response.content.decode()
        self.assertRegex(content, r'name="case_field_panel"[^>]*>TSO500</textarea>')
        self.assertRegex(content, r'name="case_field_gene"[^>]*>(RUNX1|TP53)</textarea>')

    def test_building_pins_the_report_flags_and_the_case_fields(self):
        response = self.client.post(self._url("create_case_report"), self._post_data())
        self.assertEqual(response.status_code, 200)

        case_report = CaseReport.objects.get(sample=self.proband)
        self.assertEqual(case_report.summary, "Summary interpretation.")
        # A group comes back as one dict, which is the shape a JSON template reads
        self.assertEqual(case_report.case_values,
                         {"panel": "TSO500", "assay_success": {"tmb": True, "msi": False}})
        reported_by_classification = {row.classification_modification.classification_id: row.reported
                                      for row in case_report.rows}
        self.assertTrue(reported_by_classification[self.reported.pk])
        self.assertFalse(reported_by_classification[self.unreported.pk])
        self.assertTrue(case_report.pdf_file)
        self.assertTrue(case_report.docx_file)

    def test_preview_renders_the_document_without_saving_one(self):
        response = self.client.post(self._url("preview_case_report"), self._post_data())

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "RUNX1")
        self.assertFalse(CaseReport.objects.exists())


@override_settings(CELERY_TASK_ALWAYS_EAGER=True, LIFTOVER_CLASSIFICATIONS=False,
                   CLINGEN_ALLELE_REGISTRY_LOGIN=None)
class CaseReportMeasureTickTest(ClassifyReportTestCase):
    """ The assay flags the scientist used to answer by hand: a bool case_field naming a measure and
        a tick_when rule starts ticked from what the case was measured at (sapath#454) """

    MEASURE_FIELDS = [
        {"key": "assay_success_msi", "label": "MSI", "type": "bool", "default": True,
         "group": "assay_success", "measure": "msi", "tick_when": {"called": True}},
        {"key": "assay_success_tmb", "label": "TMB", "type": "bool", "default": True,
         "group": "assay_success", "measure": "tmb", "tick_when": {"called": True}},
        {"key": "caveat_purity", "label": "Purity", "type": "bool", "default": False,
         "group": "caveats", "measure": "tumour_fraction",
         "tick_when": [{"call_in": ["Insufficient", "No tumour"]}, {"value_below": 20}]},
        {"key": "caveat_low_gis", "label": "GIS", "type": "bool", "default": False,
         "group": "caveats", "measure": "gis", "tick_when": {"value_below": 42}},
    ]

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.patient = Patient.objects.create(first_name="Measure", last_name="Tick")
        assign_permission_to_user_and_groups(cls.user, cls.patient)
        cls.specimen = Specimen.objects.create(reference_id="2600000004", patient=cls.patient)
        cls.dna = Extraction.objects.create(specimen=cls.specimen, reference_id="2600000004C",
                                            nucleic_acid_source=NucleicAcid.DNA)
        Sample.objects.filter(pk=cls.proband.pk).update(extraction=cls.dna)
        cls.template = ClassificationReportTemplate.objects.create(
            name="measure tick template", case_template=generic_case_template(),
            case_fields=cls.MEASURE_FIELDS)
        cls.classification = Classification.create(
            user=cls.user, lab=cls.lab, sample=cls.proband, source=SubmissionSource.VARIANT_GRID,
            variant=cls.variant,
            data={**READY_EVIDENCE, SpecialEKeys.GENE_SYMBOL: {"value": "RUNX1"}})
        cls.classification.publish_latest(cls.user)

    @staticmethod
    def _measures(**by_key) -> dict:
        return {key: SpecimenMeasure(value=value, unit=unit, call=call)
                for key, (value, unit, call) in by_key.items()}

    def _values(self, measures: dict, draft=None) -> dict:
        return _case_values_for_form(self.template, [], draft, measures)

    def test_each_rule_starts_its_tick_from_the_measure(self):
        measures = self._measures(msi=(2.48, "%", "Stable"), tmb=(7.1, "mut/Mb", None),
                                  tumour_fraction=(0.05, None, "Insufficient"), gis=(31.0, None, None))

        values = self._values(measures)

        self.assertTrue(values["assay_success_msi"])   # called
        self.assertFalse(values["assay_success_tmb"])  # the import left it uncalled
        self.assertTrue(values["caveat_purity"])       # call_in
        self.assertTrue(values["caveat_low_gis"])      # value_below

    def test_a_list_of_rules_ticks_when_any_holds(self):
        """ Purity is the pathologist's call from Mocha or the number off the pipeline, whichever the
            case has - a call with no number, a number under 20% with no call, and a number over it """
        self.assertTrue(self._values(self._measures(tumour_fraction=(None, "%", "No tumour")))["caveat_purity"])
        self.assertTrue(self._values(self._measures(tumour_fraction=(15.0, "%", None)))["caveat_purity"])
        self.assertFalse(self._values(self._measures(tumour_fraction=(60.0, "%", None)))["caveat_purity"])

    def test_a_measure_the_case_lacks_leaves_the_fields_default(self):
        values = self._values({})

        self.assertTrue(values["assay_success_msi"])
        self.assertFalse(values["caveat_purity"])

    def test_a_drafts_own_answer_wins_over_the_rule(self):
        """ The scientist's adjustment is the answer - a rule that disagrees does not undo it """
        draft = CaseReport.objects.create(template=self.template, lab=self.lab, user=self.user,
                                          case_values={"assay_success": {"assay_success_msi": False}},
                                          **CaseReport.source_kwargs(SampleSourceLevel.SPECIMEN,
                                                                     self.specimen))

        values = self._values(self._measures(msi=(2.48, "%", "Stable")), draft=draft)

        self.assertFalse(values["assay_success_msi"])

    def test_the_dialog_shows_each_measure_beside_its_checkbox(self):
        SpecimenMeasure.objects.create(specimen=self.specimen, measure_type=SpecimenMeasureType.MSI,
                                       value=2.48, unit="%", call="MSS",
                                       threshold="MSI-High >= 30%, MSI-Low >= 10%, MSS < 10% unstable sites",
                                       threshold_source="settings.TSO500_MSI_CALL_BANDS")
        self.client.force_login(self.user)
        url = reverse("case_report_build_dialog",
                      kwargs={"case_type": "specimen", "case_id": self.specimen.pk})

        response = self.client.post(url, {
            "report_template": self.template.pk,
            "classification_modification_id": [self.classification.last_published_version.pk],
        })

        content = response.content.decode()
        self.assertIn("2.48% (MSS)", content)
        self.assertIn("no measure", content)  # the case has no TMB, tumour fraction or GIS
        self.assertRegex(content, r'checked[^>]*id="case_field_assay_success_msi"')
        # The policy behind the tick is on the form, so a scientist can ask for it to change
        self.assertIn("ticked when the measure has a call", content)
        self.assertIn("policy MSI-High &gt;= 30%, MSI-Low &gt;= 10%, MSS &lt; 10% unstable sites (settings.TSO500_MSI_CALL_BANDS)", content)
        self.assertIn("ticked when the call is Insufficient or No tumour or the value is below 20", content)

    def test_a_field_naming_an_unknown_measure_or_rule_is_not_saved(self):
        """ The JSON is hand written in admin, so a typo says so rather than ticking nothing """
        for case_fields in ([{"key": "flag", "type": "bool", "measure": "msi_status"}],
                            [{"key": "flag", "type": "bool", "measure": "msi",
                              "tick_when": {"call_is": "Stable"}}],
                            [{"key": "flag", "type": "bool", "measure": "msi",
                              "tick_when": [{"called": True}, {}]}],
                            [{"key": "flag", "type": "bool", "tick_when": {"called": True}}]):
            with self.subTest(case_fields=case_fields):
                template = ClassificationReportTemplate(name="bad", case_fields=case_fields)
                with self.assertRaises(ValidationError):
                    template.clean()


@override_settings(CELERY_TASK_ALWAYS_EAGER=True, LIFTOVER_CLASSIFICATIONS=False,
                   CLINGEN_ALLELE_REGISTRY_LOGIN=None)
class CaseReportLibraryQCTickTest(ClassifyReportTestCase):
    """ The other five assay flags: a bool case_field naming a library QC category starts ticked from
        what DRAGEN's own QC said about that library (sapath#455) """

    QC_FIELDS = [
        {"key": "assay_success_amplifications", "label": "Amplifications", "type": "bool",
         "default": True, "group": "assay_success", "qc": "cnv", "tick_when": {"passed": True}},
        {"key": "assay_success_fusions", "label": "Fusions", "type": "bool", "default": True,
         "group": "assay_success", "qc": "rna", "tick_when": {"passed": True}},
        {"key": "caveat_fail", "label": "Fail", "type": "bool", "default": False,
         "group": "caveats", "qc": "dna",
         "tick_when": [{"passed": False}, {"completed": False}]},
    ]

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.patient = Patient.objects.create(first_name="Library", last_name="QC")
        assign_permission_to_user_and_groups(cls.user, cls.patient)
        cls.specimen = Specimen.objects.create(reference_id="2600000005", patient=cls.patient)
        cls.dna = Extraction.objects.create(specimen=cls.specimen, reference_id="2600000005C",
                                            nucleic_acid_source=NucleicAcid.DNA)
        Sample.objects.filter(pk=cls.proband.pk).update(extraction=cls.dna)
        cls.template = ClassificationReportTemplate.objects.create(
            name="library qc template", case_template=generic_case_template(),
            case_fields=cls.QC_FIELDS)
        cls.classification = Classification.create(
            user=cls.user, lab=cls.lab, sample=cls.proband, source=SubmissionSource.VARIANT_GRID,
            variant=cls.variant,
            data={**READY_EVIDENCE, SpecialEKeys.GENE_SYMBOL: {"value": "RUNX1"}})
        cls.classification.publish_latest(cls.user)

    def _library_qc(self, sequencing_run_name: str, category: str, passed=True, completed=True,
                    measured_date=None, metrics=None) -> LibraryQC:
        return LibraryQC.objects.create(sequencing_run_name=sequencing_run_name,
                                        pair_id="5_C0000005_ABCD_2600000005",
                                        specimen_reference="2600000005", specimen=self.specimen,
                                        category=category, passed=passed, completed=completed,
                                        measured_date=measured_date, metrics=metrics or {})

    def _values(self, library_qc: dict) -> dict:
        return _case_values_for_form(self.template, [], None, {}, library_qc)

    def test_the_passed_rule_starts_the_tick_from_the_librarys_qc(self):
        library_qc = {"cnv": LibraryQC(passed=True, completed=True),
                      "rna": LibraryQC(passed=False, completed=True),
                      "dna": LibraryQC(passed=True, completed=True)}

        values = self._values(library_qc)

        self.assertTrue(values["assay_success_amplifications"])
        self.assertFalse(values["assay_success_fusions"])
        self.assertFalse(values["caveat_fail"])

    def test_a_run_that_did_not_complete_ticks_the_fail_caveat(self):
        """ The caveat is either half - the library failed QC, or the run never finished for it """
        self.assertTrue(self._values({"dna": LibraryQC(passed=True, completed=False)})["caveat_fail"])
        self.assertTrue(self._values({"dna": LibraryQC(passed=False, completed=True)})["caveat_fail"])

    def test_a_category_the_case_has_no_qc_for_leaves_the_fields_default(self):
        values = self._values({})

        self.assertTrue(values["assay_success_amplifications"])
        self.assertFalse(values["caveat_fail"])

    def test_the_newest_library_is_the_specimens_qc(self):
        """ A repeat sequencing is a new library with its own QC, and it supersedes the one it replaced """
        self._library_qc("RUN_1", LibraryQCCategory.CNV, passed=False,
                         measured_date=timezone.now() - timedelta(days=7))
        latest = self._library_qc("RUN_2", LibraryQCCategory.CNV, passed=True,
                                  measured_date=timezone.now())

        library_qc = specimen_library_qc(self.specimen)

        self.assertEqual(latest, library_qc["cnv"])
        self.assertTrue(self._values(library_qc)["assay_success_amplifications"])

    def test_the_dialog_shows_each_categorys_metrics_against_its_guideline(self):
        self._library_qc("RUN_1", LibraryQCCategory.CNV, passed=True, metrics={
            "GENE_SCALED_MAD": {"value": 0.059, "unit": "Count", "lsl": 0, "usl": 0.134, "passed": True},
            "MEDIAN_BIN_COUNT_CNV_TARGET": {"value": 6.4, "unit": "Count", "lsl": 1, "usl": None, "passed": True},
        })
        self.client.force_login(self.user)
        url = reverse("case_report_build_dialog",
                      kwargs={"case_type": "specimen", "case_id": self.specimen.pk})

        response = self.client.post(url, {
            "report_template": self.template.pk,
            "classification_modification_id": [self.classification.last_published_version.pk],
        })

        content = response.content.decode()
        self.assertIn("GENE_SCALED_MAD 0.059 (&lt;= 0.134), MEDIAN_BIN_COUNT_CNV_TARGET 6.4 (&gt;= 1)", content)
        self.assertIn("no QC", content)  # the case has no RNA or DNA library QC
        self.assertRegex(content, r'checked[^>]*id="case_field_assay_success_amplifications"')
        self.assertIn("ticked when the library passed QC", content)
        self.assertIn("ticked when the library failed QC or the run did not complete for the library", content)

    def test_a_field_naming_an_unknown_category_a_measure_rule_or_both_kinds_is_not_saved(self):
        """ The JSON is hand written in admin, so a typo says so rather than ticking nothing """
        for case_fields in ([{"key": "flag", "type": "bool", "qc": "amplifications"}],
                            [{"key": "flag", "type": "bool", "qc": "cnv",
                              "tick_when": {"called": True}}],
                            [{"key": "flag", "type": "bool", "qc": "cnv", "measure": "msi",
                              "tick_when": {"passed": True}}]):
            with self.subTest(case_fields=case_fields):
                template = ClassificationReportTemplate(name="bad", case_fields=case_fields)
                with self.assertRaises(ValidationError):
                    template.clean()
