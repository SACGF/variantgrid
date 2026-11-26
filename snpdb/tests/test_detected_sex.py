"""
    Sex detected from a ploidy caller's sex karyotype, else chrX genotype stats, which the trio wizard
    compares against the patient record - prenatal cases can be entered under the mother's record (sapath#439)
"""
import os

from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase
from django.utils import timezone

from annotation.fake_data import get_fake_annotation_version
from library.genomics.ploidy_vcf_header import PloidyHeader, read_ploidy_header
from patients.enums import Sex
from patients.models import Extraction, Patient, Specimen
from snpdb.enums import ImportStatus
from snpdb.fake_data import create_fake_cohort
from snpdb.models import (
    VCF,
    CohortGenotypeCollection,
    CohortGenotypeStats,
    GenomeBuild,
    Sample,
    SampleStatsCodeVersion,
)
from snpdb.models.models_cohort_stats import MIN_CHRX_VARIANTS_FOR_SEX_GUESS
from snpdb.models.models_vcf import SamplePloidy


class TestChrXSexGuess(TestCase):
    @staticmethod
    def _stats(x_hom_count, x_het_count) -> CohortGenotypeStats:
        return CohortGenotypeStats(x_hom_count=x_hom_count, x_het_count=x_het_count)

    def test_female_ratio(self):
        self.assertEqual(Sex.FEMALE, self._stats(200, 200).chrx_sex_guess)

    def test_male_ratio(self):
        self.assertEqual(Sex.MALE, self._stats(200, 10).chrx_sex_guess)

    def test_ambiguous_ratio_unknown(self):
        self.assertEqual(Sex.UNKNOWN, self._stats(300, 100).chrx_sex_guess)

    def test_below_minimum_count_unknown(self):
        """ Targeted panels have too few chrX calls for the ratio to mean anything """
        low = MIN_CHRX_VARIANTS_FOR_SEX_GUESS - 2
        self.assertEqual(Sex.UNKNOWN, self._stats(low, 1).chrx_sex_guess)

    def test_no_hets_male(self):
        """ A caller that knows chrX is haploid in males makes no hets there at all """
        self.assertEqual(Sex.MALE, self._stats(500, 0).chrx_sex_guess)


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
        cls.code_version = SampleStatsCodeVersion.objects.create(name="SampleStats", version=1)

    def setUp(self):
        self.sample.patient = None
        self.sample.save()
        CohortGenotypeStats.objects.filter(sample=self.sample).delete()
        self.sample = type(self.sample).objects.get(pk=self.sample.pk)  # clear cached_property

    def _set_patient_sex(self, sex):
        patient = Patient.objects.create(patient_code="PROBAND", family_code="fam", sex=sex)
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


class TestPloidyHeader(TestCase):
    def test_dragen_ploidy_header(self):
        filename = os.path.join(settings.BASE_DIR, "upload", "test_data", "vcf", "dragen_ploidy",
                                "trisomy21.ploidy.vcf")
        with open(filename, encoding="utf-8") as f:
            header = "".join(line for line in f if line.startswith("##"))
        self.assertEqual(PloidyHeader(sex_karyotype="XY", reference_sex_karyotype="XY",
                                      autosome_depth_of_coverage=35.1771), read_ploidy_header(header))

    def test_no_karyotype_is_no_ploidy(self):
        self.assertIsNone(read_ploidy_header("##fileformat=VCFv4.2\n##autosomeDepthOfCoverage=35.1\n"))

    def test_karyotype_sex(self):
        expected = {"XY": Sex.MALE, "XXY": Sex.MALE, "XYY": Sex.MALE,
                    "XX": Sex.FEMALE, "XXX": Sex.FEMALE, "X0": Sex.FEMALE,
                    "UNKNOWN": Sex.UNKNOWN}
        for karyotype, sex in expected.items():
            with self.subTest(karyotype=karyotype):
                self.assertEqual(sex, SamplePloidy(sex_karyotype=karyotype).sex)


class TestSampleDetectedSexFromKaryotype(TestCase):
    """ A DRAGEN ploidy VCF is a separate upload, linked to the small variant VCF's sample by extraction """
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser_sex_karyotype')[0]
        genome_build = GenomeBuild.get_name_or_alias("GRCh38")
        patient = Patient.objects.create(first_name="KARYOTYPE", last_name="PATIENT")
        specimen = Specimen.objects.create(patient=patient, reference_id="SPEC1")
        cls.extraction = Extraction.objects.create(specimen=specimen, reference_id="EXT1")
        other_specimen = Specimen.objects.create(patient=patient, reference_id="SPEC2")
        cls.other_extraction = Extraction.objects.create(specimen=other_specimen, reference_id="EXT2")

        def _sample(name):
            vcf = VCF.objects.create(name=name, genotype_samples=1, genome_build=genome_build,
                                     import_status=ImportStatus.SUCCESS, user=cls.user, date=timezone.now())
            return Sample.objects.create(name=name, vcf=vcf)

        cls.small_variant_sample = _sample("SMP5008")
        cls.ploidy_sample = _sample("SMP5008_ploidy")
        SamplePloidy.objects.create(sample=cls.ploidy_sample, sex_karyotype="XXY", reference_sex_karyotype="XY")

    def _reloaded(self, sample) -> Sample:
        return Sample.objects.get(pk=sample.pk)  # clear cached_property

    def _link(self, sample, extraction):
        sample.extraction = extraction
        sample.save()

    def test_own_karyotype(self):
        sex, detail = self._reloaded(self.ploidy_sample).get_detected_sex()
        self.assertEqual(Sex.MALE, sex)
        self.assertIn("XXY", detail)
        self.assertNotIn("shares this sample's extraction", detail)

    def test_karyotype_from_same_extraction(self):
        self._link(self.ploidy_sample, self.extraction)
        self._link(self.small_variant_sample, self.extraction)
        sample = self._reloaded(self.small_variant_sample)
        self.assertEqual(Sex.MALE, sample.detected_sex)
        self.assertIn("SMP5008_ploidy", sample.get_detected_sex()[1])

    def test_no_karyotype_from_the_patients_other_extraction(self):
        self._link(self.ploidy_sample, self.extraction)
        self._link(self.small_variant_sample, self.other_extraction)
        self.assertEqual(Sex.UNKNOWN, self._reloaded(self.small_variant_sample).detected_sex)

    def test_no_extraction_shares_nothing(self):
        """ A null extraction isn't a shared one """
        self.assertEqual(Sex.UNKNOWN, self._reloaded(self.small_variant_sample).detected_sex)
