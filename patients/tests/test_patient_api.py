"""
#1707 - accessioning patients, specimens and extractions over the API
"""
from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import override_settings
from django.urls import reverse
from rest_framework import status
from rest_framework.test import APITestCase

from library.guardian_utils import assign_permission_to_user_and_groups
from patients.models import Extraction, ExternalModelManager, Patient, Specimen
from patients.models_enums import NucleicAcid, TissueStatus

RECONCILE_TASK = "patients.views_rest.reconcile_pending_extractions.delay"


class PatientAPITest(APITestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user(username="patient_api_user", password="x")
        cls.other_user = User.objects.create_user(username="patient_api_other", password="x")
        cls.admin_user = User.objects.create_superuser(username="patient_api_admin")
        cls.external_manager = ExternalModelManager.objects.create(name="HELIX", can_modify=True)

    def setUp(self):
        self.client.force_authenticate(user=self.user)

    def _post(self, url_name, data, user=None):
        if user:
            self.client.force_authenticate(user=user)
        with patch(RECONCILE_TASK):
            return self.client.post(reverse(url_name), data, format="json")

    def _create_patient(self, user=None, **overrides):
        data = {"patient_code": "2600000001P", "first_name": "TSO", "last_name": "PATIENT", **overrides}
        return self._post("api_patient-list", data, user=user)

    def _create_specimen(self, user=None, **overrides):
        data = {"patient": "2600000001P", "reference_id": "2600000001",
                "tissue_status": TissueStatus.AFFECTED, **overrides}
        return self._post("api_specimen-list", data, user=user)

    def _create_extraction(self, **overrides):
        data = {"specimen": "2600000001", "reference_id": "2600000001C",
                "nucleic_acid_source": NucleicAcid.DNA, **overrides}
        return self._post("api_extraction-list", data)

    def _create_run(self):
        self._create_patient()
        self._create_specimen()
        return self._create_extraction()

    def test_creates_patient_specimen_extraction(self):
        response = self._create_run()
        self.assertEqual(response.status_code, status.HTTP_201_CREATED, response.data)

        extraction = Extraction.objects.get(reference_id="2600000001C")
        self.assertEqual(extraction.specimen.reference_id, "2600000001")
        self.assertEqual(extraction.specimen.patient.patient_code, "2600000001P")
        self.assertEqual(extraction.nucleic_acid_source, NucleicAcid.DNA)

    def test_re_posting_the_same_payload_yields_one_row_of_each(self):
        self._create_run()
        self._create_run()

        self.assertEqual(Patient.objects.filter(patient_code="2600000001P").count(), 1)
        self.assertEqual(Specimen.objects.filter(reference_id="2600000001").count(), 1)
        self.assertEqual(Extraction.objects.filter(reference_id="2600000001C").count(), 1)

    def test_both_arms_post_against_one_specimen(self):
        """ nucleic_acid_source is on Extraction, so the DNA and RNA arms share their block """
        self._create_patient()
        self._create_specimen()
        self._create_extraction()
        self._create_extraction(reference_id="2600000001B", nucleic_acid_source=NucleicAcid.RNA)

        specimen = Specimen.objects.get(reference_id="2600000001")
        self.assertEqual(specimen.extraction_set.count(), 2)

    def test_extraction_create_fires_reconcile(self):
        self._create_patient()
        self._create_specimen()
        with patch(RECONCILE_TASK) as mock_delay:
            self.client.post(reverse("api_extraction-list"),
                             {"specimen": "2600000001", "reference_id": "2600000001C"}, format="json")
        mock_delay.assert_called_once()

    def test_specimen_naming_an_unknown_patient_is_400(self):
        response = self._create_specimen(patient="NOBODY")
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertIn("NOBODY", str(response.data))

    def test_specimen_by_external_pk_reference(self):
        self._create_patient(external_pk={"code": "H12345", "external_type": "HelixID",
                                          "external_manager": "HELIX"})
        response = self._create_specimen(patient={"code": "H12345", "external_type": "HelixID"})
        self.assertEqual(response.status_code, status.HTTP_201_CREATED, response.data)
        self.assertEqual(Specimen.objects.get(reference_id="2600000001").patient.patient_code,
                         "2600000001P")

    def test_external_pk_keys_the_upsert(self):
        """ A client that names an ExternalPK gets the same patient back whatever else it sends """
        external_pk = {"code": "H12345", "external_type": "HelixID", "external_manager": "HELIX"}
        self._create_patient(external_pk=external_pk)
        self._create_patient(patient_code="DIFFERENT-CODE", external_pk=external_pk)

        self.assertEqual(Patient.objects.filter(external_pk__code="H12345").count(), 1)
        self.assertEqual(Patient.objects.get(external_pk__code="H12345").patient_code, "DIFFERENT-CODE")

    def test_tissue_status_round_trips(self):
        self._create_patient()
        self._create_specimen(tissue_status=TissueStatus.REFERENCE)

        specimen = Specimen.objects.get(reference_id="2600000001")
        self.assertEqual(specimen.tissue_status, TissueStatus.REFERENCE)

        response = self.client.get(reverse("api_specimen-detail", kwargs={"pk": specimen.pk}))
        self.assertEqual(response.data["tissue_status"], TissueStatus.REFERENCE)

    def test_specimen_echoes_its_patient_identifiers(self):
        self._create_patient()
        self._create_specimen()

        specimen = Specimen.objects.get(reference_id="2600000001")
        response = self.client.get(reverse("api_specimen-detail", kwargs={"pk": specimen.pk}))
        self.assertEqual(response.data["patient"], {"reference_id": "2600000001P"})

    def test_another_users_patient_is_not_visible(self):
        self._create_patient()
        self.client.force_authenticate(user=self.other_user)

        response = self.client.get(reverse("api_patient-list"))
        self.assertEqual(response.data, [])

    def test_another_users_patient_cannot_be_attached_to(self):
        self._create_patient()
        response = self._create_specimen(user=self.other_user)
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)

    def test_unknown_external_manager_is_400_naming_the_known_ones(self):
        response = self._create_patient(external_pk={"code": "H12345", "external_type": "HelixID",
                                                     "external_manager": "HELIXX"})
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertIn("HELIX", str(response.data))
        self.assertFalse(ExternalModelManager.objects.filter(name="HELIXX").exists())

    def test_superuser_creates_an_external_manager(self):
        response = self._create_patient(user=self.admin_user,
                                        external_pk={"code": "H12345", "external_type": "HelixID",
                                                     "external_manager": "NEW_LIMS"})
        self.assertEqual(response.status_code, status.HTTP_201_CREATED, response.data)
        self.assertTrue(ExternalModelManager.objects.filter(name="NEW_LIMS").exists())

    @override_settings(PATIENTS_API_EXTERNAL_MANAGER_CREATE_ADMIN_ONLY=False)
    def test_anyone_creates_an_external_manager_where_configured(self):
        response = self._create_patient(external_pk={"code": "H12345", "external_type": "HelixID",
                                                     "external_manager": "PUBLIC_LIMS"})
        self.assertEqual(response.status_code, status.HTTP_201_CREATED, response.data)
        self.assertTrue(ExternalModelManager.objects.filter(name="PUBLIC_LIMS").exists())

    def test_code_without_external_type_is_400(self):
        self._create_patient()
        response = self._create_specimen(patient={"code": "H12345"})
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertIn("external_type", str(response.data))

    def test_anonymous_request_is_rejected(self):
        self.client.force_authenticate(user=None)
        response = self.client.get(reverse("api_patient-list"))
        self.assertIn(response.status_code, (status.HTTP_401_UNAUTHORIZED, status.HTTP_403_FORBIDDEN))


class PatientAPIPermissionTest(APITestCase):
    """ Guardian gives all three filter_for_user, so a client only reaches its own records """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.owner = User.objects.create_user(username="api_perm_owner", password="x")
        cls.other_user = User.objects.create_user(username="api_perm_other", password="x")
        cls.patient = Patient.objects.create(first_name="PERM", last_name="PATIENT",
                                             patient_code="PERM-1")
        assign_permission_to_user_and_groups(cls.owner, cls.patient)
        cls.specimen = Specimen.objects.create(patient=cls.patient, reference_id="PERMSPEC")
        cls.extraction = Extraction.objects.create(specimen=cls.specimen, reference_id="PERMSPECC")

    def test_owner_lists_their_own(self):
        self.client.force_authenticate(user=self.owner)
        for url_name in ("api_patient-list", "api_specimen-list", "api_extraction-list"):
            response = self.client.get(reverse(url_name))
            self.assertEqual(len(response.data), 1, url_name)

    def test_other_user_lists_nothing(self):
        self.client.force_authenticate(user=self.other_user)
        for url_name in ("api_patient-list", "api_specimen-list", "api_extraction-list"):
            response = self.client.get(reverse(url_name))
            self.assertEqual(response.data, [], url_name)

    def test_other_user_cannot_fetch_by_pk(self):
        self.client.force_authenticate(user=self.other_user)
        response = self.client.get(reverse("api_specimen-detail", kwargs={"pk": self.specimen.pk}))
        self.assertEqual(response.status_code, status.HTTP_404_NOT_FOUND)
