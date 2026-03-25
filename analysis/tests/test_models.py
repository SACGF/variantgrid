import unittest

from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.utils import timezone

from analysis.models import Analysis, AnalysisLock, SampleNode, GeneListNode, ZygosityNode
from annotation.fake_annotation import get_fake_annotation_version
from library.guardian_utils import assign_permission_to_user_and_groups
from snpdb.models import GenomeBuild, VCF, Sample, Cohort, CohortSample, CohortGenotypeCollection, ImportStatus


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class AnalysisModelTestCase(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        owner_username = f"test_user_{__file__}_owner"
        non_owner_username = f"test_user_{__file__}_non_owner"
        admin_username = f"test_user_{__file__}_admin"

        cls.owner_user = User.objects.get_or_create(username=owner_username)[0]
        cls.non_owner_user = User.objects.get_or_create(username=non_owner_username)[0]
        cls.admin_user = User.objects.create_superuser(admin_username)

        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)

    def test_permissions(self):
        analysis = Analysis(genome_build=self.grch37)
        analysis.set_defaults_and_save(self.owner_user)

        self.assertTrue(analysis.can_write(self.owner_user))
        self.assertTrue(analysis.can_write(self.admin_user))
        self.assertFalse(analysis.can_write(self.non_owner_user))

    def test_locking(self):
        analysis = Analysis(genome_build=self.grch37)
        analysis.set_defaults_and_save(self.owner_user)
        analysis.is_locked

        AnalysisLock.objects.create(analysis=analysis, locked=True, user=self.owner_user, date=timezone.now())
        # Bump version to expire cache
        analysis.version += 1
        analysis.save()
        self.assertTrue(analysis.is_locked)

        # Nobody should be able to write if locked
        self.assertFalse(analysis.can_write(self.owner_user))
        self.assertFalse(analysis.can_write(self.admin_user))
        self.assertFalse(analysis.can_write(self.non_owner_user))


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class AncestorSampleNoGenotypeTestCase(TestCase):
    """ Test that ancestor sample nodes work with VCFs that have no genotype samples (issue #1418) """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        cls.user = User.objects.get_or_create(username=f"test_user_{__file__}_no_geno")[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)

        # Create a VCF with no genotype samples (genotype_samples=0)
        cls.vcf = VCF.objects.create(name="no_genotype.vcf", genotype_samples=0,
                                     genome_build=cls.grch37,
                                     import_status=ImportStatus.SUCCESS,
                                     user=cls.user, date=timezone.now())
        cls.sample = Sample.objects.create(name="no_genotype.vcf", vcf=cls.vcf,
                                           import_status=ImportStatus.SUCCESS)
        assign_permission_to_user_and_groups(cls.user, cls.vcf)
        assign_permission_to_user_and_groups(cls.user, cls.sample)

        cohort = Cohort.objects.create(name="no_geno_cohort", user=cls.user, vcf=cls.vcf,
                                       genome_build=cls.grch37,
                                       import_status=ImportStatus.SUCCESS)
        CohortSample.objects.create(cohort=cohort, sample=cls.sample,
                                    cohort_genotype_packed_field_index=0, sort_order=0)
        assign_permission_to_user_and_groups(cls.user, cohort)
        CohortGenotypeCollection.objects.create(cohort=cohort,
                                                cohort_version=cohort.version,
                                                num_samples=1)

        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(cls.user)
        cls.sample_node = SampleNode.objects.create(analysis=cls.analysis, sample=cls.sample)

    def _create_child_node(self, node_class, **kwargs):
        node = node_class.objects.create(analysis=self.analysis, sample=self.sample, **kwargs)
        node.add_parent(self.sample_node)
        node._cached_parents = None  # Clear stale cache from create()'s save()
        node.save()
        return node

    def test_no_genotype_vcf_has_genotype_false(self):
        """ Sanity check: VCF with genotype_samples=0 reports has_genotype=False """
        self.assertFalse(self.vcf.has_genotype)
        self.assertFalse(self.sample.has_genotype)

    def test_ancestor_samples_includes_no_genotype_sample(self):
        """ _get_ancestor_samples should find samples from no-genotype VCFs """
        gene_list_node = self._create_child_node(GeneListNode)
        ancestor_samples = gene_list_node._get_ancestor_samples()
        self.assertIn(self.sample, ancestor_samples)

    def test_gene_list_node_no_config_error_for_no_genotype_sample(self):
        """ GeneListNode should not report a configuration error for a no-genotype ancestor sample """
        gene_list_node = self._create_child_node(GeneListNode)
        config_errors = gene_list_node._get_configuration_errors()
        ancestor_error = [e for e in config_errors if "is not set as a sample in any ancestors" in str(e)]
        self.assertEqual(ancestor_error, [])

    def test_zygosity_node_no_config_error_for_no_genotype_sample(self):
        """ ZygosityNode should not report a configuration error for a no-genotype ancestor sample """
        zygosity_node = self._create_child_node(ZygosityNode)
        config_errors = zygosity_node._get_configuration_errors()
        ancestor_error = [e for e in config_errors if "is not set as a sample in any ancestors" in str(e)]
        self.assertEqual(ancestor_error, [])

    def test_handle_ancestor_input_samples_changed_auto_sets_sample(self):
        """ Auto-sample-selection should work for no-genotype VCF samples """
        gene_list_node = GeneListNode.objects.create(analysis=self.analysis)
        gene_list_node.add_parent(self.sample_node)
        gene_list_node._cached_parents = None  # Clear stale cache from create()'s save()
        gene_list_node.version = 1  # Must be non-zero for auto-set logic
        gene_list_node.save()

        gene_list_node.handle_ancestor_input_samples_changed()
        self.assertEqual(gene_list_node.sample, self.sample)


if __name__ == '__main__':
    unittest.main()
