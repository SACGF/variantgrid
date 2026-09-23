"""
The specimen page's measures - the pathologist's tumour content off the specimen, and each sequencing
analysis' own record (#1904) and library QC (sapath#455)
"""
from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse

from library.guardian_utils import assign_permission_to_user_and_groups
from patients.models import Patient, Specimen, SpecimenMeasure
from patients.models_enums import NucleicAcid, SpecimenMeasureType
from seqauto.models import DragenTSO500CombinedVariantOutput, LibraryQC
from seqauto.models.models_enums import LibraryQCCategory


class SpecimenPageMeasuresTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user(username="measure_user", password="x")
        cls.patient = Patient.objects.create(first_name="MEASURE", last_name="PATIENT")
        assign_permission_to_user_and_groups(cls.user, cls.patient)
        cls.specimen = Specimen.objects.create(patient=cls.patient, reference_id="2600000001")

    def _specimen_page(self) -> str:
        self.client.force_login(self.user)
        response = self.client.get(reverse("view_specimen", kwargs={"specimen_id": self.specimen.pk}))
        self.assertEqual(response.status_code, 200)
        return response.content.decode()

    def test_pathology_and_sequencing_analyses_show_on_the_specimen_page(self):
        SpecimenMeasure.objects.create(specimen=self.specimen, measure_type=SpecimenMeasureType.TUMOUR_FRACTION,
                                       value=60.0, unit="%", method="Mocha")
        DragenTSO500CombinedVariantOutput.objects.create(sequencing_run_name="260101_M02027_0001_000000000-TSO500",
                                                         pair_id="5_C0000001_FCUP_2600000001",
                                                         specimen_reference="2600000001", specimen=self.specimen,
                                                         module_version="2.1.1", total_tmb=12.3)

        html = self._specimen_page()
        self.assertIn("Tumour content (pathology)", html)
        self.assertIn("60.0", html)
        self.assertIn("12.3 mut/Mb", html)
        self.assertIn("DRAGEN TSO500 CombinedVariantOutput 2.1.1", html)

    def test_library_qc_shows_on_the_specimen_page(self):
        """ What the caller's own QC said about each library sequenced off this specimen (sapath#455) """
        LibraryQC.objects.create(pair_id="5_C0000001_FCUP_2600000001",
                                 sequencing_run_name="260101_M02027_0001_000000000-TSO500",
                                 specimen_reference="2600000001", specimen=self.specimen,
                                 category=LibraryQCCategory.CNV, nucleic_acid=NucleicAcid.DNA,
                                 passed=False, completed=True,
                                 method="DRAGEN TSO500 MetricsOutput 2.6.2.4",
                                 metrics={"GENE_SCALED_MAD": {"value": 0.9, "unit": "Count",
                                                              "lsl": 0, "usl": 0.134, "passed": False}})

        html = self._specimen_page()
        self.assertIn("Library QC", html)
        self.assertIn("failed", html)
        self.assertIn("260101_M02027_0001_000000000-TSO500", html)
        self.assertIn("5_C0000001_FCUP_2600000001", html)
        self.assertIn("CNV", html)
        # The metrics themselves are on the pair's page, which the summary links to
        self.assertIn(reverse("view_tso500_pair",
                              kwargs={"sequencing_run_name": "260101_M02027_0001_000000000-TSO500",
                                      "pair_id": "5_C0000001_FCUP_2600000001"}), html)
