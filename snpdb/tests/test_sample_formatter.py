"""Tests for Sample._get_sample_formatter_params (issue #230)."""
from django.test import TestCase

from patients.models import Patient
from snpdb.models import Sample


class TestSampleFormatterPatientCode(TestCase):
    def test_patient_code_param_uses_patient_code_field(self):
        patient = Patient(pk=42, first_name="REAL", last_name="SURNAME",
                          patient_code="DEID-42")
        sample = Sample(pk=99, name="sample", patient=patient)
        params = sample._get_sample_formatter_params()
        self.assertEqual(params["patient_code"], "DEID-42")
        self.assertNotEqual(params["patient_code"], patient.last_name)

    def test_patient_code_param_blank_when_unset(self):
        patient = Patient(pk=42, first_name="REAL", last_name="SURNAME")
        sample = Sample(pk=99, name="sample", patient=patient)
        params = sample._get_sample_formatter_params()
        self.assertEqual(params["patient_code"], "")
