"""
    Sex detected from chrX genotype stats, which the trio wizard compares against the
    patient record - prenatal cases can be entered under the mother's record (sapath#439)
"""
from django.contrib.auth.models import User
from django.test import TestCase

from annotation.fake_annotation import get_fake_annotation_version
from patients.models import Patient
from patients.models_enums import Sex
from snpdb.models import (
    CohortGenotypeCollection,
    CohortGenotypeStats,
    GenomeBuild,
    SampleStatsCodeVersion,
)
from snpdb.models.models_cohort_stats import MIN_CHRX_VARIANTS_FOR_SEX_GUESS
from snpdb.tests.utils.fake_cohort_data import create_fake_cohort


class TestChrXSexGuess(TestCase):
    @staticmethod
    def _stats(x_hom_count, x_het_count) -> CohortGenotypeStats:
        return CohortGenotypeStats(x_hom_count=x_hom_count, x_het_count=x_het_count)

    def test_female_ratio(self):
        self.assertEqual(Sex.FEMALE, self._stats(10, 200).chrx_sex_guess)

    def test_male_ratio(self):
        self.assertEqual(Sex.MALE, self._stats(200, 10).chrx_sex_guess)

    def test_ambiguous_ratio_unknown(self):
        self.assertEqual(Sex.UNKNOWN, self._stats(100, 200).chrx_sex_guess)

    def test_below_minimum_count_unknown(self):
        """ Targeted panels have too few chrX calls for the ratio to mean anything """
        low = MIN_CHRX_VARIANTS_FOR_SEX_GUESS - 2
        self.assertEqual(Sex.UNKNOWN, self._stats(low, 1).chrx_sex_guess)

    def test_no_hets_unknown(self):
        self.assertEqual(Sex.UNKNOWN, self._stats(500, 0).chrx_sex_guess)


class TestSampleDetectedSex(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser_genetic_sex')[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.genome_build)
        cls.cohort = create_fake_cohort(cls.user, cls.genome_build)
        cls.sample = cls.cohort.cohortsample_set.get(sample__name='proband').sample
        cls.cgc = CohortGenotypeCollection.objects.get(cohort=cls.cohort)
        cls.code_version = SampleStatsCodeVersion.objects.create(name="SampleStats", version=1, code_git_hash="test")

    def setUp(self):
        self.sample.patient = None
        self.sample.save()
        CohortGenotypeStats.objects.filter(sample=self.sample).delete()
        self.sample = type(self.sample).objects.get(pk=self.sample.pk)  # clear cached_property

    def _set_patient_sex(self, sex):
        patient = Patient.objects.create(family_code="fam", sex=sex)
        self.sample.patient = patient
        self.sample.save()

    def _set_genotype_stats(self, x_hom_count, x_het_count):
        CohortGenotypeStats.objects.create(cohort_genotype_collection=self.cgc, sample=self.sample,
                                           code_version=self.code_version,
                                           x_hom_count=x_hom_count, x_het_count=x_het_count)

    def test_no_stats_detects_nothing(self):
        self._set_patient_sex(Sex.FEMALE)
        self.assertEqual(Sex.UNKNOWN, self.sample.detected_sex)
        self.assertEqual(Sex.FEMALE, self.sample.patient_sex)

    def test_detected_sex_can_disagree_with_patient_record(self):
        """ Male fetus imported under the mother's record - the wizard flags this """
        self._set_patient_sex(Sex.FEMALE)
        self._set_genotype_stats(x_hom_count=200, x_het_count=10)
        self.assertEqual(Sex.MALE, self.sample.detected_sex)
        self.assertEqual(Sex.FEMALE, self.sample.patient_sex)

    def test_no_patient_is_unknown(self):
        self.assertEqual(Sex.UNKNOWN, self.sample.patient_sex)
