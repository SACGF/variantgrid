from django.contrib.auth.models import User
from django.test import TestCase

from patients.models import Patient, Sex


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
        self.assertEqual(Patient.match(first_name="JESSIE", last_name="SMITH"), jessie_smith)
        self.assertEqual(Patient.match(first_name="JESSIE", last_name="SMITH", sex=Sex.UNKNOWN), jessie_smith)
        self.assertEqual(Patient.match(first_name="JESSIE", last_name="SMITH", sex=Sex.MALE), jessie_smith)
        self.assertEqual(Patient.match(first_name="JESSIE", last_name="SMITH", sex=Sex.FEMALE), jessie_smith)

        # Test matches against FEMALE
        self.assertEqual(Patient.match(first_name="BILLIE", last_name="JEAN"), billie_jean)
        self.assertEqual(Patient.match(first_name="BILLIE", last_name="JEAN", sex=Sex.UNKNOWN), billie_jean)
        self.assertEqual(Patient.match(first_name="BILLIE", last_name="JEAN", sex=Sex.MALE), None)  # No match
        self.assertEqual(Patient.match(first_name="BILLIE", last_name="JEAN", sex=Sex.FEMALE), billie_jean)

        # Test matches against MALE
        self.assertEqual(Patient.match(first_name="CHARLIE", last_name="JONES"), charlie_jones)
        self.assertEqual(Patient.match(first_name="CHARLIE", last_name="JONES", sex=Sex.UNKNOWN), charlie_jones)
        self.assertEqual(Patient.match(first_name="CHARLIE", last_name="JONES", sex=Sex.MALE), charlie_jones)
        self.assertEqual(Patient.match(first_name="CHARLIE", last_name="JONES", sex=Sex.FEMALE), None)  # No match

        # Test matches against unknown person
        self.assertEqual(Patient.match(first_name="BOB", last_name="DOBALINA"), bob_dobalina)
        self.assertEqual(Patient.match(first_name="BOB", last_name="DOBALINA", date_of_birth=None), bob_dobalina)
        self.assertEqual(Patient.match(first_name="BOB", last_name="DOBALINA", date_of_birth='2025-01-01'), bob_dobalina)

        # Test dupes
        with self.assertRaises(Patient.MultipleObjectsReturned):
            Patient.match(first_name="DAVID", last_name="LAWRENCE")

        with self.assertRaises(Patient.MultipleObjectsReturned):
            Patient.match(first_name="DAVID", last_name="LAWRENCE", date_of_birth=None)

        # Test matches
        self.assertEqual(Patient.match(first_name="DAVID", last_name="LAWRENCE", date_of_birth='1930-01-01'), None)  # Nobody
        self.assertEqual(Patient.match(first_name="DAVID", last_name="LAWRENCE", date_of_birth='1980-03-05'), dl_me)
        self.assertEqual(Patient.match(first_name="DAVID", last_name="LAWRENCE", date_of_birth='2020-01-01'), dl_other)



