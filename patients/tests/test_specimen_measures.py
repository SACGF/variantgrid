"""
#1559 - specimen-level measures (TMB, MSI, GIS), transcribed by the lab client out of vendor output
"""
from django.contrib.auth.models import User
from django.urls import reverse
from rest_framework import status
from rest_framework.test import APITestCase

from library.guardian_utils import assign_permission_to_user_and_groups
from patients.models import Extraction, Patient, Specimen, SpecimenMeasure
from patients.models_enums import NucleicAcid, SpecimenMeasureType

TMB_PAYLOAD = {"measure_type": SpecimenMeasureType.TMB, "value": 12.3, "unit": "mut/Mb",
               "call": "High", "threshold": ">=10 mut/Mb", "threshold_source": "lab policy 2026",
               "method": "DRAGEN TSO500 v2.6.2",
               "source_payload": {"file": "MetricsOutput.tsv", "TMB": "12.3"}}
MSI_PAYLOAD = {"measure_type": SpecimenMeasureType.MSI, "value": 0.05, "unit": "%", "call": "Stable"}
GIS_PAYLOAD = {"measure_type": SpecimenMeasureType.GIS, "value": 42.0, "method": "PhenoHRD"}


class SpecimenMeasureAPITest(APITestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user(username="measure_user", password="x")
        cls.other_user = User.objects.create_user(username="measure_other", password="x")
        cls.patient = Patient.objects.create(first_name="MEASURE", last_name="PATIENT")
        assign_permission_to_user_and_groups(cls.user, cls.patient)
        cls.specimen = Specimen.objects.create(patient=cls.patient, reference_id="2600000001")
        cls.dna = Extraction.objects.create(specimen=cls.specimen, reference_id="2600000001C",
                                            nucleic_acid_source=NucleicAcid.DNA)

    def setUp(self):
        self.client.force_authenticate(user=self.user)

    def _bulk_post(self, measures, specimen="2600000001"):
        return self.client.post(reverse("api_specimen_measure_bulk_create"),
                                {"specimen": specimen, "measures": measures}, format="json")

    def test_bulk_post_of_a_runs_measures(self):
        response = self._bulk_post([TMB_PAYLOAD, MSI_PAYLOAD, GIS_PAYLOAD])
        self.assertEqual(response.status_code, status.HTTP_201_CREATED, response.data)

        measures = SpecimenMeasure.objects.filter(specimen=self.specimen)
        self.assertEqual(measures.count(), 3)
        tmb = measures.get(measure_type=SpecimenMeasureType.TMB)
        self.assertEqual(tmb.value, 12.3)
        self.assertEqual(tmb.call, "High")
        self.assertEqual(tmb.threshold_source, "lab policy 2026")
        self.assertEqual(tmb.user, self.user)

    def test_source_payload_round_trips(self):
        self._bulk_post([TMB_PAYLOAD])
        tmb = SpecimenMeasure.objects.get(specimen=self.specimen, measure_type=SpecimenMeasureType.TMB)
        self.assertEqual(tmb.source_payload, TMB_PAYLOAD["source_payload"])

    def test_re_posting_a_measure_replaces_it(self):
        """ The report wants one current MSI rather than a history to choose from """
        self._bulk_post([MSI_PAYLOAD])
        original = SpecimenMeasure.objects.get(specimen=self.specimen,
                                               measure_type=SpecimenMeasureType.MSI)

        self._bulk_post([{**MSI_PAYLOAD, "value": 21.5, "call": "High"}])
        measures = SpecimenMeasure.objects.filter(specimen=self.specimen,
                                                  measure_type=SpecimenMeasureType.MSI)
        self.assertEqual(measures.count(), 1)
        updated = measures.get()
        self.assertEqual(updated.pk, original.pk)
        self.assertEqual(updated.value, 21.5)
        self.assertEqual(updated.call, "High")
        self.assertGreater(updated.modified, original.modified)

    def test_measure_naming_an_unknown_specimen_is_400(self):
        response = self._bulk_post([TMB_PAYLOAD], specimen="NO-SUCH-SPECIMEN")
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertIn("NO-SUCH-SPECIMEN", str(response.data))
        self.assertFalse(SpecimenMeasure.objects.exists())

    def test_measure_can_name_the_arm_it_came_off(self):
        self._bulk_post([{**TMB_PAYLOAD, "extraction": "2600000001C"}])
        tmb = SpecimenMeasure.objects.get(specimen=self.specimen, measure_type=SpecimenMeasureType.TMB)
        self.assertEqual(tmb.extraction, self.dna)

    def test_single_measure_endpoint(self):
        response = self.client.post(reverse("api_specimen_measure-list"),
                                    {"specimen": "2600000001", **GIS_PAYLOAD}, format="json")
        self.assertEqual(response.status_code, status.HTTP_201_CREATED, response.data)
        self.assertEqual(SpecimenMeasure.objects.get().measure_type, SpecimenMeasureType.GIS)

    def test_another_users_specimen_is_not_reachable(self):
        self.client.force_authenticate(user=self.other_user)
        response = self._bulk_post([TMB_PAYLOAD])
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertFalse(SpecimenMeasure.objects.exists())

    def test_measures_are_listed_for_their_owner_only(self):
        self._bulk_post([TMB_PAYLOAD])
        response = self.client.get(reverse("api_specimen_measure-list"))
        self.assertEqual(len(response.data), 1)

        self.client.force_authenticate(user=self.other_user)
        self.assertEqual(self.client.get(reverse("api_specimen_measure-list")).data, [])

    def test_measures_show_on_the_specimen_page(self):
        self._bulk_post([TMB_PAYLOAD, MSI_PAYLOAD])
        self.client.force_authenticate(user=None)
        self.client.force_login(self.user)

        response = self.client.get(reverse("view_specimen", kwargs={"specimen_id": self.specimen.pk}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        html = response.content.decode()
        self.assertIn("Tumour mutational burden", html)
        self.assertIn("12.3", html)
        self.assertIn("Stable", html)
