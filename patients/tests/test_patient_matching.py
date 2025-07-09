from django.test import TestCase

from patients.models import Patient, Sex
from patients.models_enums import PatientRecordMatchType


class TestPatientMatching(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        # Delete any existing stuff
        patients = [
            # Test UniSex matching
            Patient(first_name="unrelated", last_name="rando", date_of_birth=None),
            Patient(first_name="JESSIE", last_name="SMITH", sex=Sex.UNKNOWN, date_of_birth=None),
            Patient(first_name="BILLIE", last_name="JEAN", sex=Sex.FEMALE),
            Patient(first_name="CHARLIE", last_name="JONES", sex=Sex.MALE),
            # Test DOB provided
            Patient(first_name="BOB", last_name="DOBALINA", date_of_birth=None),
            Patient(first_name="DAVID", last_name="LAWRENCE", date_of_birth='1980-03-05'),
            Patient(first_name="DAVID", last_name="LAWRENCE", date_of_birth='2020-01-01'),
        ]

        Patient.objects.bulk_create(patients)

    def test_patient_match_sex(self):
        jessie_smith = Patient.objects.get(first_name="JESSIE", last_name="SMITH", sex=Sex.UNKNOWN)
        billie_jean = Patient.objects.get(first_name="BILLIE", last_name="JEAN", sex=Sex.FEMALE)
        charlie_jones = Patient.objects.get(first_name="CHARLIE", last_name="JONES", sex=Sex.MALE)
        # Test DOB provided
        bob_dobalina = Patient.objects.get(first_name="BOB", last_name="DOBALINA", date_of_birth=None)
        dl_me = Patient.objects.get(first_name="DAVID", last_name="LAWRENCE", date_of_birth='1980-03-05')
        dl_other = Patient.objects.get(first_name="DAVID", last_name="LAWRENCE", date_of_birth='2020-01-01')

        # Test matches against UNKNOWN
        self.assertEqual(Patient.match(first_name="JESSIE", last_name="SMITH"),
                         (jessie_smith, PatientRecordMatchType.EXACT))
        self.assertEqual(Patient.match(first_name="JESSIE", last_name="SMITH", sex=Sex.UNKNOWN),
                         (jessie_smith, PatientRecordMatchType.EXACT))
        self.assertEqual(Patient.match(first_name="JESSIE", last_name="SMITH", sex=Sex.MALE),
                         (jessie_smith, PatientRecordMatchType.PARTIAL))
        self.assertEqual(Patient.match(first_name="JESSIE", last_name="SMITH", sex=Sex.FEMALE),
                         (jessie_smith, PatientRecordMatchType.PARTIAL))

        # Test matches against FEMALE
        self.assertEqual(Patient.match(first_name="BILLIE", last_name="JEAN"),
                         (billie_jean, PatientRecordMatchType.EXACT))
        self.assertEqual(Patient.match(first_name="BILLIE", last_name="JEAN", sex=Sex.UNKNOWN),
                         (billie_jean, PatientRecordMatchType.EXACT))
        self.assertEqual(Patient.match(first_name="BILLIE", last_name="JEAN", sex=Sex.MALE),
                         (None, None))  # No match
        self.assertEqual(Patient.match(first_name="BILLIE", last_name="JEAN", sex=Sex.FEMALE),
                         (billie_jean, PatientRecordMatchType.EXACT))

        # Test matches against MALE
        self.assertEqual(Patient.match(first_name="CHARLIE", last_name="JONES"),
                         (charlie_jones, PatientRecordMatchType.EXACT))
        self.assertEqual(Patient.match(first_name="CHARLIE", last_name="JONES", sex=Sex.UNKNOWN),
                         (charlie_jones, PatientRecordMatchType.EXACT))
        self.assertEqual(Patient.match(first_name="CHARLIE", last_name="JONES", sex=Sex.MALE),
                         (charlie_jones, PatientRecordMatchType.EXACT))
        self.assertEqual(Patient.match(first_name="CHARLIE", last_name="JONES", sex=Sex.FEMALE),
                         (None, None))  # No match

        # Test matches against unknown person
        self.assertEqual(Patient.match(first_name="BOB", last_name="DOBALINA"),
                         (bob_dobalina, PatientRecordMatchType.EXACT))
        self.assertEqual(Patient.match(first_name="BOB", last_name="DOBALINA", date_of_birth=None),
                         (bob_dobalina, PatientRecordMatchType.EXACT))
        self.assertEqual(Patient.match(first_name="BOB", last_name="DOBALINA", date_of_birth='2025-01-01'),
                         (bob_dobalina, PatientRecordMatchType.PARTIAL))

        # Test dupes
        with self.assertRaises(Patient.MultipleObjectsReturned):
            Patient.match(first_name="DAVID", last_name="LAWRENCE")

        with self.assertRaises(Patient.MultipleObjectsReturned):
            Patient.match(first_name="DAVID", last_name="LAWRENCE", date_of_birth=None)

        # Test matches
        self.assertEqual(Patient.match(first_name="DAVID", last_name="LAWRENCE", date_of_birth='1930-01-01'),
                         (None, None))  # Nobody
        self.assertEqual(Patient.match(first_name="DAVID", last_name="LAWRENCE", date_of_birth='1980-03-05'),
                     (dl_me, PatientRecordMatchType.EXACT))
        self.assertEqual(Patient.match(first_name="DAVID", last_name="LAWRENCE", date_of_birth='2020-01-01'),
                         (dl_other, PatientRecordMatchType.EXACT))
