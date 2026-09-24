from django.contrib.auth.models import User
from django.test import TestCase
from django.utils import timezone

from annotation.phenotype_matcher import PhenotypeMatcher
from library.guardian_utils import assign_permission_to_user_and_groups
from ontology.tests.test_data_ontology import (
    create_ontology_test_data,
    create_test_ontology_version,
)
from patients.models import Patient
from patients.signals.patient_search import (
    patient_preview_phenotype_extra,
    sample_preview_patient_extra,
)
from snpdb.models import VCF, GenomeBuild, ImportStatus, Sample


class TestSamplePreviewPatientExtra(TestCase):

    @classmethod
    def setUpTestData(cls):
        create_ontology_test_data()
        create_test_ontology_version()
        cls.user = User.objects.get_or_create(username='testuser_sample_preview')[0]
        cls.other_user = User.objects.get_or_create(username='testuser_sample_preview_other')[0]

        cls.patient = Patient(patient_code="PT-1", phenotype="Raised TSH")
        cls.patient.save(phenotype_matcher=PhenotypeMatcher())
        assign_permission_to_user_and_groups(cls.user, cls.patient)

        vcf = VCF.objects.create(name="sample_preview_vcf", genotype_samples=1,
                                 genome_build=GenomeBuild.get_name_or_alias("GRCh37"),
                                 import_status=ImportStatus.SUCCESS, user=cls.user, date=timezone.now())
        cls.sample = Sample.objects.create(name="sample_preview", vcf=vcf, patient=cls.patient)

    def _extras(self, user) -> dict:
        extras = sample_preview_patient_extra(sender=Sample, user=user, obj=self.sample) or []
        return {kv.key: kv.value for kv in extras}

    def test_patient_and_phenotype(self):
        self.assertEqual(self._extras(self.user), {"Patient": "PT-1", "HPO": "Increased thyroid-stimulating hormone level"})

    def test_patient_the_user_cannot_view_is_left_out(self):
        self.assertEqual(self._extras(self.other_user), {})

    def test_patient_preview_phenotype(self):
        extras = patient_preview_phenotype_extra(sender=Patient, user=self.user, obj=self.patient)
        self.assertEqual({kv.key: kv.value for kv in extras}, {"HPO": "Increased thyroid-stimulating hormone level"})
