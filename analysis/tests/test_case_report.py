""" The case report half of the Classify & Report tab (#444).

    Ordering, amp_tier and the renderers are covered over fake modifications in
    classification/tests/report/; what needs a real case is here: which samples a specimen or
    extraction case is, who may act on a built report, and what finalising writes.
"""
from django.contrib.auth.models import User
from django.test import override_settings
from django.urls import reverse

from analysis.classify_report import ClassifyReportCase, ReportCandidate
from analysis.tests.test_classify_report import READY_EVIDENCE, ClassifyReportTestCase
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
from classification.report.default_templates import generic_case_template
from library.case_report_delivery import CaseReportDelivery
from library.guardian_utils import assign_permission_to_user_and_groups
from patients.models import Extraction, Patient, Specimen
from patients.models_enums import NucleicAcid, SampleSourceLevel
from patients.sample_grouping import get_sample_group
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
