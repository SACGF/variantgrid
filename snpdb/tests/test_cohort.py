from django.contrib.auth.models import User
from django.test import TestCase

from annotation.fake_annotation import get_fake_annotation_version
from snpdb.models import Cohort, CohortGenotypeCollection, CohortSample, GenomeBuild, Trio
from snpdb.tasks.cohort_genotype_tasks import create_cohort_genotype_and_launch_task
from snpdb.tests.utils.fake_cohort_data import create_fake_trio
from snpdb.views.vcf_cohort_page import _family_groups_by_sample_id


class CohortGenotypeTestCase(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        owner_username = f"test_user_{__file__}_owner"
        cls.user_owner = User.objects.get_or_create(username=owner_username)[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.trio1 = create_fake_trio(cls.user_owner, cls.grch37)
        cls.trio2 = create_fake_trio(cls.user_owner, cls.grch37)
        get_fake_annotation_version(cls.grch37)

    def test_sub_cohort(self):
        """ A cohort is a sub-cohort only while every sample comes from the one parent cohort """
        cohort = Cohort.objects.create(name="sub_cohort", user=self.user_owner, genome_build=self.grch37)
        samples = list(self.trio1.get_samples())
        cohort.add_sample(samples[0].pk)
        cohort.add_sample(samples[1].pk)
        create_cohort_genotype_and_launch_task(cohort, run_async=False)
        self.assertTrue(cohort.is_sub_cohort)

        samples2 = list(self.trio2.get_samples())
        cohort.add_sample(samples2[0].pk)  # Add a sample that is NOT in parent cohort
        self.assertFalse(cohort.is_sub_cohort)

    def test_set_samples_single_version_bump(self):
        """ A batch membership edit costs one version bump, whatever the size of the diff """
        cohort = Cohort.objects.create(name="batch_edit", user=self.user_owner, genome_build=self.grch37)
        samples = list(self.trio1.get_samples())
        samples2 = list(self.trio2.get_samples())
        cohort.set_samples([s.pk for s in samples])

        version = cohort.version
        # 3 removed (all of trio1), 3 added (all of trio2)
        cohort.set_samples([s.pk for s in samples2])
        self.assertEqual(version + 1, cohort.version)
        self.assertEqual([s.pk for s in samples2], list(cohort.get_cohort_samples().values_list("sample_id", flat=True)))

    def test_set_samples_sort_order_matches_payload(self):
        cohort = Cohort.objects.create(name="sorted", user=self.user_owner, genome_build=self.grch37)
        samples = list(self.trio1.get_samples())
        cohort.set_samples([s.pk for s in samples])

        reversed_sample_ids = [s.pk for s in reversed(samples)]
        cohort.set_samples(reversed_sample_ids)
        cohort_samples = cohort.get_cohort_samples()
        self.assertEqual(reversed_sample_ids, [cs.sample_id for cs in cohort_samples])
        self.assertEqual([0, 1, 2], [cs.sort_order for cs in cohort_samples])

    def test_set_samples_sub_cohort_copies_parent_packed_index(self):
        parent_cohort = self.trio1.cohort
        sub_cohort = parent_cohort.create_sub_cohort(self.user_owner, list(parent_cohort.get_samples())[:2])

        third_sample = list(parent_cohort.get_samples())[2]
        sample_ids = list(sub_cohort.get_cohort_samples().values_list("sample_id", flat=True)) + [third_sample.pk]
        sub_cohort.set_samples(sample_ids)

        self.assertTrue(sub_cohort.is_sub_cohort)
        parent_packed_index = dict(parent_cohort.cohortsample_set.values_list(
            "sample_id", "cohort_genotype_packed_field_index"))
        # Packed indexes need to be unique in the sub cohort and valid in the parent's CohortGenotype arrays
        packed_indexes = list(sub_cohort.cohortsample_set.values_list("cohort_genotype_packed_field_index", flat=True))
        self.assertEqual(len(sample_ids), len(set(packed_indexes)))
        self.assertLessEqual(max(packed_indexes), max(parent_packed_index.values()))

    def test_set_samples_detaches_parent_cohort(self):
        parent_cohort = self.trio1.cohort
        sub_cohort = parent_cohort.create_sub_cohort(self.user_owner, list(parent_cohort.get_samples())[:2])
        self.assertTrue(sub_cohort.is_sub_cohort)

        outside_sample = list(self.trio2.get_samples())[0]
        sample_ids = list(sub_cohort.get_cohort_samples().values_list("sample_id", flat=True)) + [outside_sample.pk]
        sub_cohort.set_samples(sample_ids)
        self.assertFalse(sub_cohort.is_sub_cohort)

    def test_set_samples_requires_a_sample(self):
        cohort = Cohort.objects.create(name="empty", user=self.user_owner, genome_build=self.grch37)
        self.assertRaises(ValueError, cohort.set_samples, [])

    def test_set_samples_within_one_vcf_becomes_sub_cohort(self):
        """ A batch edit that lands inside a single VCF reuses that cohort's genotypes instead of rebuilding """
        cohort = Cohort.objects.create(name="converts", user=self.user_owner, genome_build=self.grch37)
        samples = list(self.trio1.get_samples())
        samples2 = list(self.trio2.get_samples())
        cohort.set_samples([samples[0].pk, samples2[0].pk])
        create_cohort_genotype_and_launch_task(cohort, run_async=False)
        self.assertFalse(cohort.is_sub_cohort)

        cohort.set_samples([s.pk for s in samples])
        create_cohort_genotype_and_launch_task(cohort, run_async=False)
        self.assertEqual(self.trio1.cohort, cohort.parent_cohort)
        self.assertFalse(CohortGenotypeCollection.objects.filter(cohort=cohort,
                                                                 cohort_version=cohort.version).exists())

    def test_family_groups_warned_about_before_removal(self):
        """ Duo/Trio/Quad cascade off CohortSample, so the editor names what a removal would delete """
        cohort = Cohort.objects.create(name="has_a_trio", user=self.user_owner, genome_build=self.grch37)
        samples = list(self.trio1.get_samples())
        cohort.set_samples([s.pk for s in samples])
        mother, father, proband = cohort.get_cohort_samples()
        trio = Trio.objects.create(name="fam01", user=self.user_owner, cohort=cohort,
                                   mother=mother, father=father, proband=proband)

        family_groups = _family_groups_by_sample_id(cohort)
        self.assertEqual({s.pk for s in samples}, set(family_groups))
        for labels in family_groups.values():
            self.assertEqual([f"Trio '{trio.name}'"], labels)

    def test_cohort_genotype_packed_field_index(self):
        """ Add/Remove CohortSamples - ensure cohort_genotype_packed_field_index stays in range """
        # Add samples from multiple parent cohorts - to ensure it's not a sub cohort
        cohort = Cohort.objects.create(name="not_a_sub_cohort", user=self.user_owner, genome_build=self.grch37)
        samples = list(self.trio1.get_samples())
        samples2 = list(self.trio2.get_samples())
        cohort.add_sample(samples[0].pk)
        cohort.add_sample(samples[1].pk)
        cohort.add_sample(samples2[0].pk)

        # Delete existing cohort sample
        CohortSample.objects.get(cohort=cohort, sample=samples[1]).delete()

        create_cohort_genotype_and_launch_task(cohort, run_async=False)
        cs = cohort.cohortsample_set.order_by("-cohort_genotype_packed_field_index").first()
        self.assertLess(cs.cohort_genotype_packed_field_index, cohort.cohort_genotype_collection.num_samples)
