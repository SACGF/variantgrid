"""
Patient.merge - folding one patient into another, and what blocks it
"""
from datetime import date

from django.test import TestCase

from patients.models import Patient, PatientModification, PatientRecordOriginType, Specimen
from patients.models_enums import Sex


class PatientMergeTest(TestCase):

    def test_merge(self):
        patient = Patient.objects.create(last_name="Smith", date_of_birth=date(1980, 1, 2))
        other = Patient.objects.create(patient_code="C1", sex=Sex.FEMALE)
        Specimen.objects.create(patient=other, reference_id="S1")
        PatientModification.objects.create(patient=other, origin=PatientRecordOriginType.INTERNAL_VG)

        self.assertEqual(patient.get_merge_conflicts(other), [])  # sex unset (default) on patient
        patient.merge(other, "test", PatientRecordOriginType.INTERNAL_VG)

        patient.refresh_from_db()
        self.assertEqual((patient.patient_code, patient.sex, patient.last_name), ("C1", Sex.FEMALE, "Smith"))
        self.assertTrue(patient.specimen_set.filter(reference_id="S1").exists())
        self.assertEqual(patient.patientmodification_set.count(), 2)
        self.assertFalse(Patient.objects.filter(pk=other.pk).exists())

    def test_conflicts_block_merge(self):
        patient = Patient.objects.create(last_name="Smith")
        other = Patient.objects.create(last_name="Jones")
        Specimen.objects.create(patient=patient, reference_id="S1")
        Specimen.objects.create(patient=other, reference_id="S1")

        self.assertEqual(len(patient.get_merge_conflicts(other)), 2)
        with self.assertRaises(ValueError):
            patient.merge(other, "test", PatientRecordOriginType.INTERNAL_VG)
        self.assertTrue(Patient.objects.filter(pk=other.pk).exists())
