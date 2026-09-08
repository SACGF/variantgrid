"""
Somalier: what we hand to `somalier extract` (issue #183 - the exported depths decide how somalier
genotypes), whether a relate needs --unknown, and that a failing stage never takes the import down.
"""
import os
import tempfile

import cyvcf2
from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.tests.inheritance_node_mixin import make_cohort_genotype
from snpdb.models import (
    CohortGenotype,
    GenomeBuild,
    ProcessingStatus,
    Sample,
    SomalierAllSamplesRelate,
    SomalierCohortRelate,
    SomalierRelatePairs,
    SomalierSampleExtract,
    SomalierVCFExtract,
    Zygosity,
)
from snpdb.tasks.somalier_tasks import _load_somalier_pairs, somalier_vcf_id
from snpdb.tests.utils.fake_cohort_data import create_fake_cohort
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

    def _flipped_site_records(self):
        """ A variant whose ALT sorts before its REF, so somalier reads it against the other allele """
        variant = slowly_create_test_variant("3", 4000, "T", "A", self.genome_build)
        make_cohort_genotype(self.cohort.cohort_genotype_collection, variant, "ROE",
                             allele_depth=[10, 0, 25], allele_frequency=[0.45, 0.0, 1.0])
        _, records = self._export()
        return next(r for r in records if r.POS == 4000)

    def test_allele_depths_written_in_somalier_site_order(self):
        """ Only the AD pair compensates - REF, ALT and GT stay as called (brentp/somalier#163) """
        self.vcf.read_depth_field = "DP"
        self.vcf.save()

        flipped = self._flipped_site_records()
        self.assertEqual(("T", ["A"]), (flipped.REF, flipped.ALT), "REF is the reference base")
        self.assertEqual(["0/0:10,20", "1/1:0,30", "0/1:25,5"], self._calls(flipped))

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


PAIRS_HEADER = "#sample_a\tsample_b\trelatedness\tibs0\tibs2\thom_concordance\thets_a\thets_b\thets_ab\t" \
               "shared_hets\thom_alts_a\thom_alts_b\tshared_hom_alts\tn\tx_ibs0\tx_ibs2\texpected_relatedness"


def _pairs_row(sample_a: str, sample_b: str, relatedness: float, shared_hets: int, shared_hom_alts: int) -> str:
    return f"{sample_a}\t{sample_b}\t{relatedness}\t10\t20\t0.9\t50\t50\t60\t{shared_hets}\t" \
           f"30\t30\t{shared_hom_alts}\t1000\t1\t2\t-1"


class SomalierAllSamplesPairsTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.get_or_create(username='somalier_pairs_user')[0]
        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.cohort = create_fake_cohort(cls.user, genome_build)
        cls.proband = cls.cohort.cohortsample_set.get(sample__name="proband").sample
        cls.mother = cls.cohort.cohortsample_set.get(sample__name="mother").sample

    def _write_pairs(self, rows: list[str]) -> str:
        filename = tempfile.mkstemp(suffix=".pairs.tsv")[1]
        with open(filename, "w") as f:
            f.write("\n".join([PAIRS_HEADER, *rows]) + "\n")
        return filename

    def _name(self, sample: Sample) -> str:
        return f"{sample.name}_{sample.pk}"

    def test_loads_only_pairs_over_thresholds_for_existing_samples(self):
        deleted_id = Sample.objects.order_by("-pk").first().pk + 1
        rows = [
            _pairs_row(self._name(self.proband), self._name(self.mother), 0.5, 2000, 300),
            _pairs_row(self._name(self.proband), f"gone_{deleted_id}", 0.5, 2000, 300),
            _pairs_row(self._name(self.mother), self._name(self.proband), 0.05, 2000, 300),
        ]
        relate = SomalierAllSamplesRelate.objects.create(status=ProcessingStatus.PROCESSING)
        self.assertEqual(1, _load_somalier_pairs(relate, self._write_pairs(rows)))
        pair = SomalierRelatePairs.objects.get()
        self.assertEqual(self.proband.pk, pair.sample_a_id)
        self.assertEqual(self.mother.pk, pair.sample_b_id)

    def test_previous_runs_pairs_are_replaced(self):
        old_relate = SomalierAllSamplesRelate.objects.create(status=ProcessingStatus.SUCCESS)
        _load_somalier_pairs(old_relate, self._write_pairs([
            _pairs_row(self._name(self.proband), self._name(self.mother), 0.5, 2000, 300)]))

        new_relate = SomalierAllSamplesRelate.objects.create(status=ProcessingStatus.PROCESSING)
        _load_somalier_pairs(new_relate, self._write_pairs([
            _pairs_row(self._name(self.proband), self._name(self.mother), 0.9, 2000, 300)]))

        pair = SomalierRelatePairs.objects.get()
        self.assertEqual(new_relate, pair.relate)
        self.assertAlmostEqual(0.9, pair.relatedness)


@override_settings(SOMALIER={"enabled": True, "admin_only": False,
                             "vcf_base_dir": "/tmp/somalier_test", "report_base_dir": "/tmp/somalier_test",
                             "annotation_base_dir": "/tmp/somalier_test",
                             "annotation": {"command": "somalier", "ancestry_labels": "labels.tsv",
                                            "ancestry_somalier_dir": "1kg-somalier",
                                            "sites": {"GRCh37": "sites.GRCh37.vcf.gz"}},
                             "ancestry_enabled": True, "min_genotyped_sites": 100,
                             "all_samples_relate_hour": 2,
                             "relatedness": {"min_relatedness": 0.1, "min_shared_hets": 1000,
                                             "min_shared_hom_alts": 200}})
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
