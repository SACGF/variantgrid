""" The copy number call a sample's genotype autopopulates - which key it lands under is the VCF's
    declared field, since CN is an absolute count and SM/FC are ratios against the normal. """
from django.contrib.auth.models import User
from django.test import TestCase

from classification.autopopulate_evidence_keys.evidence_from_sample_and_patient import (
    get_copy_number_evidence,
)
from classification.enums import SpecialEKeys
from snpdb.models import CohortGenotype, GenomeBuild
from snpdb.tests.utils.fake_cohort_data import create_fake_cohort
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


class CopyNumberEvidenceTest(TestCase):

    def setUp(self):
        user = User.objects.get_or_create(username='testuser')[0]
        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        self.cohort = create_fake_cohort(user, genome_build)
        self.sample = self.cohort.cohortsample_set.get(sample__name="proband").sample
        self.variant = slowly_create_test_variant("3", 128198980, 'G', 'C', genome_build)

    def _copy_number_evidence(self, copy_number_field, value):
        vcf = self.sample.vcf
        vcf.copy_number_field = copy_number_field
        vcf.save()
        sample_format = [{copy_number_field: [value]} if copy_number_field else {}, {}, {}]
        cg = CohortGenotype.objects.create(collection=self.cohort.cohort_genotype_collection,
                                           variant=self.variant, samples_zygosity="OOO",
                                           format=sample_format)
        return get_copy_number_evidence(cg.get_sample_genotype(self.sample))

    def test_a_copy_ratio_is_a_fold_change(self):
        self.assertEqual((SpecialEKeys.FOLD_CHANGE, 4.31428), self._copy_number_evidence("SM", 4.31428))

    def test_a_fold_change_field_is_a_fold_change(self):
        self.assertEqual((SpecialEKeys.FOLD_CHANGE, 0.532), self._copy_number_evidence("FC", 0.532))

    def test_a_copy_number_is_an_integer_count(self):
        """ A caller may write the count as 12.0 - the key is an integer one """
        self.assertEqual((SpecialEKeys.COPY_NUMBER, 12), self._copy_number_evidence("CN", 12.0))

    def test_nothing_is_populated_without_a_declared_field(self):
        self.assertIsNone(self._copy_number_evidence(None, None))
