"""
Somalier: what we hand to `somalier extract` (issue #183 - the exported depths decide how somalier
genotypes), what `relate` is run with (--unknown, and the sites VCF of brentp/somalier#163), how ancestry
is batched across VCFs (#1850), and that a failing stage never takes the import down.
"""
import os
import shutil
import tempfile
from pathlib import Path
from unittest.mock import patch

import cyvcf2
from django.conf import settings
from django.contrib.auth.models import User
from django.core.management import call_command
from django.core.management.base import CommandError
from django.test import TestCase, override_settings
from django.urls import reverse
from django.utils import timezone

from library.guardian_utils import assign_permission_to_user_and_groups
from patients.models import Patient
from snpdb.enums import ImportStatus, ProcessingStatus, SuperPopulationCode
from snpdb.fake_data import create_fake_cohort, make_cohort_genotype
from snpdb.management.commands.somalier_existing_vcfs import Command as SomalierExistingVCFsCommand
from snpdb.management.commands.somalier_existing_vcfs import batch_by_sample_count
from snpdb.models import (
    VCF,
    Cohort,
    CohortGenotype,
    CohortSample,
    GenomeBuild,
    Sample,
    SomalierAllSamplesRelate,
    SomalierAncestry,
    SomalierAncestryBatch,
    SomalierAncestryRun,
    SomalierCohortRelate,
    SomalierConfig,
    SomalierPatientRelate,
    SomalierPatientRelatePair,
    SomalierRelatePairs,
    SomalierSampleExtract,
    SomalierVCFExtract,
    Zygosity,
)
from snpdb.models.models_somalier import get_same_individual_relatedness
from snpdb.tasks.somalier_tasks import (
    ANCESTRY_TSV_FILENAME,
    SOMALIER_PAIRS_FILENAME,
    _create_somalier_pairs,
    _load_ancestry_tsv,
    _load_somalier_pairs,
    _read_somalier_pairs,
    somalier_ancestry_batch,
    somalier_patient_relate,
    somalier_vcf_id,
)
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant
from snpdb.variants_to_vcf import _allele_depths, vcf_export_to_file

MISSING = -1  # CohortGenotype.MISSING_NUMBER_VALUE
MISSING_AF = -1.0  # the AF array is floats, so its missing value must be too


class SomalierVCFExportTest(TestCase):
    """ Samples are (proband, mother, father) in packed order """

    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.get_or_create(username='somalier_export_user')[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.cohort = create_fake_cohort(cls.user, cls.genome_build)
        cls.vcf = cls.cohort.vcf
        # the fixture leaves vcf_sample_name blank, and that's what the export writes as the column name
        for sample in cls.vcf.sample_set.all():
            sample.vcf_sample_name = sample.name
            sample.save()
        cgc = cls.cohort.cohort_genotype_collection

        # het / hom-ref / hom-alt, all depths present
        make_cohort_genotype(cgc, slowly_create_test_variant("3", 1000, "A", "T", cls.genome_build),
                             "ERO", allele_depth=[10, 0, 25], allele_frequency=[0.45, 0.0, 1.0])
        # no read depth: ref comes from AF, or is undecidable when AF is 0
        make_cohort_genotype(cgc, slowly_create_test_variant("3", 2000, "A", "T", cls.genome_build),
                             "EEU", allele_depth=[10, 10, MISSING], allele_frequency=[0.5, 0.0, MISSING_AF])
        # a missing alt depth has no honest AD at all
        make_cohort_genotype(cgc, slowly_create_test_variant("3", 3000, "A", "T", cls.genome_build),
                             "EEE", allele_depth=[MISSING, 5, 5], allele_frequency=[0.0, 0.2, 0.2])

    def _export(self) -> tuple[list[str], list[cyvcf2.Variant]]:
        with tempfile.TemporaryDirectory() as tmp_dir:
            filename = os.path.join(tmp_dir, "somalier.vcf.bgz")
            vcf_export_to_file(self.vcf, filename)
            reader = cyvcf2.Reader(filename)
            return reader.raw_header.split("\n"), list(reader)

    def _calls(self, record) -> list[str]:
        """ The sample columns exactly as written """
        return str(record).rstrip("\n").split("\t")[9:]

    def test_gt_ad_export_derives_ref_depth(self):
        self.vcf.read_depth_field = "DP"
        self.vcf.save()  # make_cohort_genotype writes samples_read_depth of 30

        header, records = self._export()
        self.assertTrue(any(line.startswith("##FORMAT=<ID=AD,") for line in header))
        self.assertTrue(any(line.startswith("##contig=<ID=3,") for line in header),
                        "contig IDs must be names, matching the CHROM we write and the sites file")

        # ref = DP - AD
        self.assertEqual(["0/1:20,10", "0/0:30,0", "1/1:5,25"], self._calls(records[0]))
        # unknown zygosity is not called at all
        self.assertEqual(["0/1:20,10", "0/1:20,10", "./.:."], self._calls(records[1]))
        # a -1 alt depth is missing, not a depth of -1
        self.assertEqual(["0/1:.", "0/1:25,5", "0/1:25,5"], self._calls(records[2]))

    def test_ref_depth_from_allele_frequency_when_no_read_depth(self):
        self.vcf.read_depth_field = None
        self.vcf.allele_frequency_field = "AF"
        self.vcf.save()
        # a VCF with no DP stored none - the format decision is per VCF, the derivation per call
        CohortGenotype.objects.filter(collection=self.cohort.cohort_genotype_collection) \
            .update(samples_read_depth=[MISSING] * 3)

        _, records = self._export()
        # AF 0.45 of 10 alt reads -> 12 ref; a hom-ref call with AF 0 has no derivable ref depth
        self.assertEqual(["0/1:12,10", "0/0:.", "1/1:0,25"], self._calls(records[0]))
        self.assertEqual(["0/1:10,10", "0/1:.", "./.:."], self._calls(records[1]))

    def test_gt_only_export_when_no_depth_fields(self):
        """ Declaring AD we can't fill in would zero out every sample, so don't declare it """
        self.vcf.allele_depth_field = None
        self.vcf.read_depth_field = None
        self.vcf.allele_frequency_field = None
        self.vcf.save()

        header, records = self._export()
        self.assertFalse(any(line.startswith("##FORMAT=<ID=AD,") for line in header))
        self.assertEqual(["0/1", "0/0", "1/1"], self._calls(records[0]))
        self.assertEqual(["0/1", "0/1", "./."], self._calls(records[1]))

    def test_zygosity_counts_returned_per_sample(self):
        counts = vcf_export_to_file(self.vcf, tempfile.mkstemp(suffix=".vcf.bgz")[1])
        proband = self.cohort.cohortsample_set.get(sample__name="proband").sample
        father = self.cohort.cohortsample_set.get(sample__name="father").sample
        self.assertEqual(3, counts[proband][Zygosity.HET])
        self.assertEqual(1, counts[father][Zygosity.HOM_ALT])
        self.assertEqual(1, counts[father][Zygosity.UNKNOWN_ZYGOSITY])

    def test_allele_frequency_percent(self):
        self.vcf.allele_frequency_percent = True
        # 25% of the reads are alt, so 10 alt reads means 30 ref
        self.assertEqual((30, 10), _allele_depths(self.vcf, 10, None, 25.0))

    def _flipped_site_records(self, compensate_allele_order=True):
        """ A variant whose ALT sorts before its REF, which a pre-0.3.5 somalier reads against the
            other allele - the deployment setting decides, so never the box the tests run on """
        variant = slowly_create_test_variant("3", 4000, "T", "A", self.genome_build)
        make_cohort_genotype(self.cohort.cohort_genotype_collection, variant, "ROE",
                             allele_depth=[10, 0, 25], allele_frequency=[0.45, 0.0, 1.0])
        with override_settings(SOMALIER={**settings.SOMALIER,
                                         "compensate_allele_order": compensate_allele_order}):
            _, records = self._export()
        return next(r for r in records if r.POS == 4000)

    def test_allele_depths_written_in_somalier_site_order(self):
        """ Only the AD pair compensates - REF, ALT and GT stay as called (brentp/somalier#163) """
        self.vcf.read_depth_field = "DP"
        self.vcf.save()

        flipped = self._flipped_site_records()
        self.assertEqual(("T", ["A"]), (flipped.REF, flipped.ALT), "REF is the reference base")
        self.assertEqual(["0/0:10,20", "1/1:0,30", "0/1:25,5"], self._calls(flipped))

    def test_allele_depths_left_as_called_for_a_fixed_somalier(self):
        """ v0.3.5+ reads the site's REF/ALT out of the sites VCF at relate, so nothing is swapped """
        self.vcf.read_depth_field = "DP"
        self.vcf.save()

        flipped = self._flipped_site_records(compensate_allele_order=False)
        self.assertEqual(["0/0:20,10", "1/1:30,0", "0/1:5,25"], self._calls(flipped))

    def test_genotype_carries_it_when_there_are_no_depths(self):
        """ Nothing else can, so a depth-less VCF flips the genotype instead (brentp/somalier#163) """
        self.vcf.allele_depth_field = None
        self.vcf.read_depth_field = None
        self.vcf.allele_frequency_field = None
        self.vcf.save()

        flipped = self._flipped_site_records()
        self.assertEqual(("T", ["A"]), (flipped.REF, flipped.ALT))
        self.assertEqual(["1/1", "0/0", "0/1"], self._calls(flipped))


class SomalierRelateTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.get_or_create(username='somalier_relate_user')[0]
        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.cohort = create_fake_cohort(cls.user, genome_build)
        cls.vcf_extract = SomalierVCFExtract.objects.create(vcf=cls.cohort.vcf)
        for sample in cls.cohort.vcf.sample_set.all():
            SomalierSampleExtract.objects.create(vcf_extract=cls.vcf_extract, sample=sample,
                                                 ref_count=100, het_count=50, hom_count=50)

    def test_has_hom_ref_calls(self):
        relate = SomalierCohortRelate.objects.create(cohort=self.cohort, cohort_version=self.cohort.version)
        self.assertTrue(relate.has_hom_ref_calls)

        SomalierSampleExtract.objects.filter(sample__name="mother").update(ref_count=0)
        self.assertFalse(relate.has_hom_ref_calls, "A VCF with no 0/0 for a sample needs --unknown")

    def test_all_samples_relate_never_joint_called(self):
        self.assertFalse(SomalierAllSamplesRelate.objects.create().has_hom_ref_calls)

    def test_all_samples_relate_uses_the_build_most_samples_are_in(self):
        """ relate takes one sites file and these samples span builds """
        grch38 = GenomeBuild.get_name_or_alias("GRCh38")
        create_fake_cohort(self.user, grch38)  # 3 GRCh37 samples in setUpTestData, 3 here
        create_fake_cohort(self.user, grch38)
        Sample.objects.all().update(import_status=ImportStatus.SUCCESS)
        self.assertEqual(grch38, SomalierAllSamplesRelate.objects.create().genome_build)


class SomalierRelateSitesTest(TestCase):
    """ relate is handed the sites VCF exactly when we're not compensating in the AD order (#163) """

    def setUp(self):
        self.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        self.annotation_dir = tempfile.mkdtemp()
        self.sites_filename = os.path.join(self.annotation_dir, "sites.GRCh37.vcf.gz")

    def _somalier_settings(self, compensate_allele_order: bool) -> dict:
        return {**settings.SOMALIER, "annotation_base_dir": self.annotation_dir,
                "annotation": {**settings.SOMALIER["annotation"],
                               "sites": {"GRCh37": "sites.GRCh37.vcf.gz"}},
                "compensate_allele_order": compensate_allele_order}

    def test_no_sites_arg_while_compensating(self):
        with override_settings(SOMALIER=self._somalier_settings(True)):
            self.assertEqual([], SomalierConfig().get_relate_sites_args(self.genome_build))

    def test_sites_passed_when_not_compensating(self):
        Path(self.sites_filename).touch()
        with override_settings(SOMALIER=self._somalier_settings(False)):
            self.assertEqual(["--sites", self.sites_filename],
                             SomalierConfig().get_relate_sites_args(self.genome_build))

    def test_missing_sites_file_is_fatal_when_not_compensating(self):
        with override_settings(SOMALIER=self._somalier_settings(False)):
            with self.assertRaises(ValueError):
                SomalierConfig().get_relate_sites_args(self.genome_build)


SOMALIER_TEST_DATA = Path(__file__).parent / "test_data" / "somalier"
FIRST_SAMPLE_PK = 9_000_001  # SOMALIER_TEST_DATA files name samples somalier-style, name_pk, from this pk


def _create_vcf_samples(user: User, first_pk: int = FIRST_SAMPLE_PK,
                        names: tuple[str, ...] = ("proband", "mother", "father")) -> list[Sample]:
    """ Fixed pks, so the samples' somalier names match the SOMALIER_TEST_DATA files """
    vcf = VCF.objects.create(name=f"somalier_{first_pk}", genotype_samples=len(names),
                             genome_build=GenomeBuild.get_name_or_alias("GRCh37"),
                             import_status=ImportStatus.SUCCESS, user=user, date=timezone.now())
    assign_permission_to_user_and_groups(user, vcf)
    samples = []
    for pk, name in enumerate(names, start=first_pk):
        sample = Sample.objects.create(pk=pk, vcf=vcf, vcf_sample_name=name, name=name)
        assign_permission_to_user_and_groups(user, sample)
        samples.append(sample)
    return samples


def _pairs(filename: str) -> Path:
    return SOMALIER_TEST_DATA / filename


def _related_dir(filename: str) -> Path:
    """ A relate output dir holding a SOMALIER_TEST_DATA pairs file """
    related_dir = Path(tempfile.mkdtemp())
    shutil.copy(_pairs(filename), related_dir / SOMALIER_PAIRS_FILENAME)
    return related_dir


class SomalierAllSamplesPairsTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.get_or_create(username='somalier_pairs_user')[0]
        cls.proband, cls.mother, cls.father = _create_vcf_samples(cls.user)

    def test_loads_only_pairs_over_thresholds_for_existing_samples(self):
        """ Also has a pair under min_relatedness, and one with a deleted sample """
        relate = SomalierAllSamplesRelate.objects.create(status=ProcessingStatus.PROCESSING)
        self.assertEqual(1, _load_somalier_pairs(relate, _pairs("pairs_thresholds.tsv")))
        pair = SomalierRelatePairs.objects.get()
        self.assertEqual(self.proband.pk, pair.sample_a_id)
        self.assertEqual(self.mother.pk, pair.sample_b_id)

    def test_loads_somalier_0_3_5_concordance_column(self):
        relate = SomalierAllSamplesRelate.objects.create(status=ProcessingStatus.PROCESSING)
        self.assertEqual(1, _load_somalier_pairs(relate, _pairs("pairs_somalier_0_3_5.tsv")))
        self.assertAlmostEqual(0.9, SomalierRelatePairs.objects.get().hom_concordance)

    def test_previous_runs_pairs_are_replaced(self):
        old_relate = SomalierAllSamplesRelate.objects.create(status=ProcessingStatus.SUCCESS)
        _load_somalier_pairs(old_relate, _pairs("pairs_thresholds.tsv"))

        new_relate = SomalierAllSamplesRelate.objects.create(status=ProcessingStatus.PROCESSING)
        _load_somalier_pairs(new_relate, _pairs("pairs_same_individual.tsv"))

        pair = SomalierRelatePairs.objects.get()
        self.assertEqual(new_relate, pair.relate)
        self.assertAlmostEqual(0.98, pair.relatedness)

    def test_get_for_sample_only_pairs_user_can_view(self):
        """ #1843 - the other sample of a pair is named on the page, so it must be one the user can view.
            somalier's order puts the proband on either side of a pair """
        other_user = User.objects.get_or_create(username='somalier_pairs_other_user')[0]
        _create_vcf_samples(other_user, first_pk=FIRST_SAMPLE_PK + 3, names=("hidden",))
        relate = SomalierAllSamplesRelate.objects.create(status=ProcessingStatus.PROCESSING)
        _load_somalier_pairs(relate, _pairs("pairs_visibility.tsv"))
        pairs = SomalierRelatePairs.get_for_sample(self.proband, self.user)
        self.assertEqual({(self.proband.pk, self.mother.pk), (self.father.pk, self.proband.pk)},
                         set(pairs.values_list("sample_a", "sample_b")))

    def test_same_individual_relatedness(self):
        """ #196 - a patient's samples: a missing pair between samples in the last relate is unrelated, a sample
            that relate couldn't have compared is listed rather than paired """
        vcf_extract = SomalierVCFExtract.objects.create(vcf=self.proband.vcf, status=ProcessingStatus.SUCCESS)
        for sample, het_count in [(self.proband, 5000), (self.mother, 5000), (self.father, 10)]:
            SomalierSampleExtract.objects.create(vcf_extract=vcf_extract, sample=sample,
                                                 het_count=het_count, hom_count=5000)
        relate = SomalierAllSamplesRelate.objects.create(status=ProcessingStatus.SUCCESS)
        # somalier's own order, mother before proband
        _load_somalier_pairs(relate, _pairs("pairs_same_individual.tsv"))

        samples = [self.proband, self.mother, self.father]
        relatedness = get_same_individual_relatedness(samples)
        self.assertEqual([(self.father, "Too few HET/HOM sites")], relatedness.not_checked)
        self.assertEqual([True], [sp.same_individual for sp in relatedness.sample_pairs])
        self.assertFalse(relatedness.mismatch)

        SomalierRelatePairs.objects.all().delete()
        relatedness = get_same_individual_relatedness(samples)
        self.assertIsNone(relatedness.sample_pairs[0].pair)
        self.assertTrue(relatedness.mismatch)

        vcf_extract.save()  # re-extracted since the relate
        relatedness = get_same_individual_relatedness(samples)
        self.assertEqual([], relatedness.sample_pairs)
        self.assertFalse(relatedness.mismatch)

    def test_patient_relate_pair_used_while_newer_than_extracts(self):
        """ A patient relate keeps every pair, so a low one is a different individual rather than missing """
        vcf_extract = SomalierVCFExtract.objects.create(vcf=self.proband.vcf, status=ProcessingStatus.SUCCESS)
        for sample in [self.proband, self.mother]:
            SomalierSampleExtract.objects.create(vcf_extract=vcf_extract, sample=sample,
                                                 het_count=5000, hom_count=5000)
        patient_relate = SomalierPatientRelate.objects.create(patient=Patient.objects.create(patient_code="RELATE"),
                                                              status=ProcessingStatus.SUCCESS)
        _create_somalier_pairs(SomalierPatientRelatePair, patient_relate,
                               _read_somalier_pairs(_pairs("pairs_different_individual.tsv")))

        sample_pair = get_same_individual_relatedness([self.proband, self.mother]).sample_pairs[0]
        self.assertIsInstance(sample_pair.pair, SomalierPatientRelatePair)
        self.assertFalse(sample_pair.same_individual)

        vcf_extract.save()  # re-extracted since the patient relate, but not the nightly one after it
        all_samples_relate = SomalierAllSamplesRelate.objects.create(status=ProcessingStatus.SUCCESS)
        _load_somalier_pairs(all_samples_relate, _pairs("pairs_same_individual.tsv"))
        sample_pair = get_same_individual_relatedness([self.proband, self.mother]).sample_pairs[0]
        self.assertIsInstance(sample_pair.pair, SomalierRelatePairs)
        self.assertTrue(sample_pair.same_individual)


class SomalierVCFTaskTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.get_or_create(username='somalier_task_user')[0]
        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.cohort = create_fake_cohort(cls.user, genome_build)
        cls.vcf = cls.cohort.vcf

    def test_failure_is_recorded_not_raised(self):
        """ #432 - the sites VCF isn't imported here, so the export blows up before somalier runs """
        somalier_vcf_id(self.vcf.pk)
        vcf_extract = SomalierVCFExtract.objects.get(vcf=self.vcf)
        self.assertEqual(ProcessingStatus.ERROR, vcf_extract.status)
        self.assertIn("Traceback", vcf_extract.error_exception)

    def test_concurrent_run_is_skipped(self):
        processing = SomalierVCFExtract.objects.create(vcf=self.vcf, status=ProcessingStatus.PROCESSING)
        somalier_vcf_id(self.vcf.pk)
        vcf_extract = SomalierVCFExtract.objects.get(vcf=self.vcf)
        self.assertEqual(processing.pk, vcf_extract.pk)
        self.assertEqual(ProcessingStatus.PROCESSING, vcf_extract.status)


class SomalierAncestryBatchTest(TestCase):
    """ #1850 - one somalier ancestry call over several VCFs """

    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.get_or_create(username='somalier_ancestry_user')[0]
        cls.vcf_extracts = []
        for i in range(3):
            samples = _create_vcf_samples(cls.user, first_pk=FIRST_SAMPLE_PK + 3 * i)
            vcf_extract = SomalierVCFExtract.objects.create(vcf=samples[0].vcf, status=ProcessingStatus.SUCCESS)
            for sample in samples:
                SomalierSampleExtract.objects.create(vcf_extract=vcf_extract, sample=sample,
                                                     het_count=5000, hom_count=5000)
            cls.vcf_extracts.append(vcf_extract)

    def setUp(self):
        report_dir = tempfile.mkdtemp()
        self.enterContext(override_settings(SOMALIER={**settings.SOMALIER, "report_base_dir": report_dir,
                                                      "vcf_base_dir": report_dir, "ancestry_enabled": True,
                                                      "min_genotyped_sites": 100}))

    @staticmethod
    def _copy_tsv(batch: SomalierAncestryBatch):
        report_dir = Path(batch.get_report_dir())
        report_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(SOMALIER_TEST_DATA / "somalier-ancestry.tsv", report_dir / ANCESTRY_TSV_FILENAME)

    def _fake_ancestry(self, batch: SomalierAncestryBatch, _runs, _processing_dir):
        self._copy_tsv(batch)
        batch.status = ProcessingStatus.SUCCESS
        batch.save()

    def test_load_ancestry_tsv_splits_samples_by_run(self):
        """ somalier-ancestry.tsv has a different ancestry per VCF, and a sample that isn't in any run """
        batch = SomalierAncestryBatch.objects.create()
        runs = [SomalierAncestryRun.objects.create(vcf_extract=ve, batch=batch) for ve in self.vcf_extracts]
        self._copy_tsv(batch)

        _load_ancestry_tsv(batch, runs)
        for run, ancestry in zip(runs, [SuperPopulationCode.EUR, SuperPopulationCode.EAS, SuperPopulationCode.AFR]):
            self.assertEqual({s.pk for s in run.get_samples()},
                             set(run.somalierancestry_set.values_list("sample_extract__sample", flat=True)))
            self.assertEqual({ancestry}, set(run.somalierancestry_set.values_list("predicted_ancestry", flat=True)))
        self.assertEqual(9, SomalierAncestry.objects.count())

    def test_vcf_page_shows_only_its_own_samples_ancestry(self):
        """ The batch's somalier report covers every VCF in it - a VCF's page must not reveal the others """
        batch = SomalierAncestryBatch.objects.create(status=ProcessingStatus.SUCCESS)
        runs = [SomalierAncestryRun.objects.create(vcf_extract=ve, batch=batch, status=ProcessingStatus.SUCCESS)
                for ve in self.vcf_extracts]
        self._copy_tsv(batch)
        _load_ancestry_tsv(batch, runs)

        vcf = self.vcf_extracts[0].vcf
        cohort = Cohort.objects.create(name="ancestry_page", user=self.user, vcf=vcf, genome_build=vcf.genome_build,
                                       import_status=ImportStatus.SUCCESS)
        for i, sample in enumerate(vcf.sample_set.order_by("pk")):
            CohortSample.objects.create(cohort=cohort, sample=sample, cohort_genotype_packed_field_index=i, sort_order=i)
        self.client.force_login(self.user)
        with override_settings(SOMALIER={**settings.SOMALIER, "enabled": True, "admin_only": False}):
            response = self.client.get(reverse("view_vcf", kwargs={"vcf_id": vcf.pk}))
        self.assertEqual(200, response.status_code)
        ancestry_sample_ids = {a.sample_extract.sample_id for a in response.context["somalier_ancestries"]}
        self.assertEqual({s.pk for s in vcf.sample_set.all()}, ancestry_sample_ids)
        self.assertNotContains(response, str(batch.uuid))

    def test_too_few_sites_skipped_while_neighbour_runs(self):
        too_few, runs_ok = self.vcf_extracts[:2]
        too_few.somaliersampleextract_set.update(het_count=10, hom_count=10)
        with patch("snpdb.tasks.somalier_tasks._somalier_ancestry", side_effect=self._fake_ancestry):
            somalier_ancestry_batch([too_few.pk, runs_ok.pk])

        skipped = SomalierAncestryRun.objects.get(vcf_extract=too_few)
        self.assertEqual(ProcessingStatus.SKIPPED, skipped.status)
        self.assertIsNone(skipped.batch)
        ran = SomalierAncestryRun.objects.get(vcf_extract=runs_ok)
        self.assertEqual(ProcessingStatus.SUCCESS, ran.status)
        self.assertEqual(3, ran.somalierancestry_set.count())

    def test_all_skipped_makes_no_batch(self):
        SomalierSampleExtract.objects.update(het_count=10, hom_count=10)
        with patch("snpdb.tasks.somalier_tasks._somalier_ancestry") as mock_ancestry:
            somalier_ancestry_batch([ve.pk for ve in self.vcf_extracts])
        mock_ancestry.assert_not_called()
        self.assertFalse(SomalierAncestryBatch.objects.exists())
        self.assertEqual(3, SomalierAncestryRun.objects.filter(status=ProcessingStatus.SKIPPED).count())

    def test_failure_marks_batch_and_every_run(self):
        """ The command isn't in the (empty) annotation dir, so somalier never runs """
        with override_settings(SOMALIER={**settings.SOMALIER, "annotation_base_dir": tempfile.mkdtemp()}):
            somalier_ancestry_batch([ve.pk for ve in self.vcf_extracts])
        batch = SomalierAncestryBatch.objects.get()
        self.assertEqual(ProcessingStatus.ERROR, batch.status)
        statuses = set(SomalierAncestryRun.objects.filter(batch=batch).values_list("status", flat=True))
        self.assertEqual({ProcessingStatus.ERROR}, statuses)
        self.assertEqual(3, SomalierAncestryRun.objects.count())

    def test_batch_by_sample_count(self):
        self.assertEqual([[1], [2, 3]], batch_by_sample_count([(1, 600), (2, 600), (3, 300)], 1000))
        self.assertEqual([[1], [2]], batch_by_sample_count([(1, 1500), (2, 10)], 1000),
                         "an extract over the limit is a batch of its own")



class SomalierExistingVCFsCommandTest(TestCase):
    """ A backfill whose extracts all failed must not look finished to the upgrader """

    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.get_or_create(username='somalier_backfill_user')[0]
        cls.vcf = _create_vcf_samples(cls.user)[0].vcf

    def setUp(self):
        self.annotation_dir = tempfile.mkdtemp()
        self.enterContext(override_settings(SOMALIER={**settings.SOMALIER, "enabled": True,
                                                      "annotation_base_dir": self.annotation_dir}))

    def test_missing_binary_fails_before_clear(self):
        SomalierVCFExtract.objects.create(vcf=self.vcf, status=ProcessingStatus.SUCCESS)
        with self.assertRaisesRegex(CommandError, "download_annotation_data"):
            call_command("somalier_existing_vcfs", "--clear")
        self.assertTrue(SomalierVCFExtract.objects.filter(vcf=self.vcf).exists())

    def test_failed_extract_is_retried_and_fails_the_command(self):
        somalier_bin = Path(SomalierConfig().get_annotation("command"))
        somalier_bin.touch(mode=0o755)
        SomalierVCFExtract.objects.create(vcf=self.vcf, status=ProcessingStatus.ERROR)

        def fail_extract(vcf_id):
            SomalierVCFExtract.objects.filter(vcf_id=vcf_id).delete()
            SomalierVCFExtract.objects.create(vcf_id=vcf_id, status=ProcessingStatus.ERROR)

        module = "snpdb.management.commands.somalier_existing_vcfs"
        with patch.object(SomalierExistingVCFsCommand, "_require_worker"), \
                patch(f"{module}.somalier_vcf_extract") as mock_extract, \
                patch(f"{module}.somalier_ancestry_batch"), \
                patch(f"{module}.somalier_all_samples"):
            mock_extract.delay.side_effect = fail_extract
            with self.assertRaisesRegex(CommandError, "extracts .* failed"):
                call_command("somalier_existing_vcfs")
        mock_extract.delay.assert_any_call(self.vcf.pk)


class SomalierPatientRelateTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.get_or_create(username='somalier_patient_user')[0]
        cls.patient = Patient.objects.create(patient_code="SOMALIER-1")
        cls.samples = _create_vcf_samples(cls.user)
        vcf_extract = SomalierVCFExtract.objects.create(vcf=cls.samples[0].vcf, status=ProcessingStatus.SUCCESS)
        for sample in cls.samples:
            sample.patient = cls.patient
            sample.save()
            SomalierSampleExtract.objects.create(vcf_extract=vcf_extract, sample=sample,
                                                 het_count=5000, hom_count=5000)

    def test_relates_once_per_change_of_samples(self):
        related_dirs = [_related_dir("pairs_patient_trio.tsv"), _related_dir("pairs_patient_mother_father.tsv")]
        with patch("snpdb.tasks.somalier_tasks._somalier_relate", side_effect=related_dirs) as mock_relate:
            somalier_patient_relate(self.patient.pk)
            somalier_patient_relate(self.patient.pk)
            self.assertEqual(1, mock_relate.call_count, "Already a relate of these samples")
            self.assertEqual(3, SomalierPatientRelatePair.objects.count(), "Unrelated pairs are kept")

            self.samples[0].no_dna_control = True
            self.samples[0].save()
            somalier_patient_relate(self.patient.pk)
            self.assertEqual(2, mock_relate.call_count)

        relate = SomalierPatientRelate.objects.get()  # the earlier one was replaced
        self.assertEqual(ProcessingStatus.SUCCESS, relate.status)
        self.assertEqual(1, relate.somalierpatientrelatepair_set.count())

    @override_settings(SOMALIER={**settings.SOMALIER, "enabled": True})
    def test_saving_an_extracted_sample_queues_its_patient(self):
        with patch("snpdb.signals.somalier_patient_relate.somalier_patient_relate") as mock_task:
            with self.captureOnCommitCallbacks(execute=True):
                self.samples[0].save()
        mock_task.si.assert_called_once_with(self.patient.pk)
