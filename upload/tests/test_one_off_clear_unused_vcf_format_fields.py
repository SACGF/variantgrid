""" one_off_clear_unused_vcf_format_fields picks its candidates from sample stats alone (#2833) """
from django.contrib.auth.models import User
from django.test import TestCase

from annotation.tasks.calculate_sample_stats import get_sample_stats_code_version
from snpdb.fake_data import create_fake_cohort
from snpdb.models import CohortGenotypeStats, GenomeBuild
from upload.management.commands.one_off_clear_unused_vcf_format_fields import (
    get_no_genotype_call_vcf_ids,
    get_uncounted_unknown_call_vcf_ids,
)


class NoGenotypeCallVCFsTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        user = User.objects.create(username="no_genotype_call_user")
        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        code_version = get_sample_stats_code_version()

        def vcf_with_stats(name, genotype_field="GT", **counts):
            cohort = create_fake_cohort(user, genome_build, name=name)
            vcf = cohort.vcf
            vcf.genotype_field = genotype_field
            vcf.save()
            for sample in vcf.sample_set.all():
                CohortGenotypeStats.objects.create(cohort_genotype_collection=cohort.cohort_genotype_collection,
                                                   sample=sample, code_version=code_version, **counts)
            return vcf

        cls.no_calls_vcf = vcf_with_stats("no_calls", unk_count=9)
        cls.called_vcf = vcf_with_stats("called", het_count=3, unk_count=9)
        cls.no_gt_field_vcf = vcf_with_stats("no_gt_field", genotype_field=None, unk_count=9, variant_count=9)
        cls.no_gt_field_uncounted_vcf = vcf_with_stats("no_gt_field_uncounted", genotype_field=None, unk_count=9)

    def test_only_a_gt_field_with_no_calls_is_a_candidate(self):
        vcf_ids = get_no_genotype_call_vcf_ids()
        self.assertIn(self.no_calls_vcf.pk, vcf_ids)
        self.assertNotIn(self.called_vcf.pk, vcf_ids)
        self.assertNotIn(self.no_gt_field_vcf.pk, vcf_ids)
        self.assertNotIn(self.no_gt_field_uncounted_vcf.pk, vcf_ids)

    def test_no_gt_field_with_uncounted_unknown_calls_is_recounted(self):
        vcf_ids = get_uncounted_unknown_call_vcf_ids()
        self.assertIn(self.no_gt_field_uncounted_vcf.pk, vcf_ids)
        self.assertNotIn(self.no_gt_field_vcf.pk, vcf_ids)
        self.assertNotIn(self.no_calls_vcf.pk, vcf_ids)
