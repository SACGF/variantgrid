"""
Sequencing run page Patients tab - sample sheet rows grouped patient -> extraction
"""
from django.contrib.auth.models import User
from django.test import TestCase
from django.utils import timezone

from patients.models import Extraction, Patient, Specimen
from patients.models_enums import MatchStatus
from seqauto.models import SampleFromSequencingSample
from seqauto.tests.test_extraction_link import make_sample_sheet, make_sequencing_run
from seqauto.views import _get_sequencing_run_patients
from snpdb.models import VCF, Sample


class SequencingRunPatientsTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_superuser(username="run_patients_user")
        sequencing_run = make_sequencing_run("RUN_PATIENTS")
        cls.sample_sheet, (cls.ss_claimed, cls.ss_hand_set, cls.ss_unmatched) = make_sample_sheet(
            sequencing_run, ["CLAIMED", "HAND_SET", "UNMATCHED"])

        cls.patient = Patient.objects.create(first_name="RUN", last_name="PATIENT")
        specimen = Specimen.objects.create(patient=cls.patient, reference_id="S1")
        cls.extraction = Extraction.objects.create(specimen=specimen, reference_id="S1A")

        # Extraction only claimed by the sequencing sample, not yet carried down to the Sample
        cls.ss_claimed.extraction = cls.extraction
        cls.ss_claimed.save()
        cls._make_sample(cls.ss_claimed)
        # Patient set by hand on the Sample, no extraction
        cls._make_sample(cls.ss_hand_set, patient=cls.patient)
        cls.ss_unmatched.extraction_match_status = MatchStatus.PENDING
        cls.ss_unmatched.save()

    @classmethod
    def _make_sample(cls, sequencing_sample, patient=None):
        vcf = VCF.objects.create(name=sequencing_sample.sample_name, date=timezone.now(), user=cls.user,
                                 genotype_samples=1)
        sample = Sample.objects.create(vcf=vcf, name=sequencing_sample.sample_name,
                                       vcf_sample_name=sequencing_sample.sample_name, patient=patient)
        SampleFromSequencingSample.objects.create(sample=sample, sequencing_sample=sequencing_sample)
        return sample

    def test_grouping(self):
        data = _get_sequencing_run_patients(self.sample_sheet, self.user)
        (run_patient,) = data["run_patients"]
        self.assertEqual(run_patient["patient"], self.patient)
        rows_by_extraction = dict(run_patient["extractions"])
        self.assertEqual([r.sequencing_sample for r in rows_by_extraction[self.extraction]], [self.ss_claimed])
        self.assertEqual([r.sequencing_sample for r in rows_by_extraction[None]], [self.ss_hand_set])

        (unmatched,) = data["unmatched_run_samples"]
        self.assertIsNone(unmatched.sample)
        self.assertEqual(unmatched.match_record.extraction_match_status, MatchStatus.PENDING)
