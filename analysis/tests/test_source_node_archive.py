"""
Source-node archive tolerance: when a source's VCF is archived, the node surfaces
a configuration error and the cohort_genotype_collection mirror returns None.

@see claude/issue_1536_data_archive_plan.md §3
"""

from django.contrib.auth.models import User
from django.test import TestCase
from django.utils import timezone

from analysis.models.nodes.cohort_mixin import CohortMixin
from snpdb.archive import DataArchivedError
from snpdb.models import GenomeBuild
from snpdb.tests.utils.fake_cohort_data import create_fake_cohort, create_fake_trio, create_fake_quad


class SourceNodeArchiveToleranceTests(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.create_user(username="src_archive_user")
        cls.gb = GenomeBuild.grch38()

    def _archive(self, vcf):
        vcf.data_archived_date = timezone.now()
        vcf.data_archived_by = self.user
        vcf.data_archive_reason = "test reason"
        vcf.save()

    def test_chokepoint_raises_for_archived_vcf_cohort(self):
        cohort = create_fake_cohort(self.user, self.gb)
        self._archive(cohort.vcf)
        cohort = type(cohort).objects.get(pk=cohort.pk)
        with self.assertRaises(DataArchivedError):
            _ = cohort.cohort_genotype_collection

    def test_chokepoint_raises_for_archived_vcf_trio(self):
        trio = create_fake_trio(self.user, self.gb)
        self._archive(trio.cohort.vcf)
        cohort = type(trio.cohort).objects.get(pk=trio.cohort.pk)
        with self.assertRaises(DataArchivedError):
            _ = cohort.cohort_genotype_collection

    def test_chokepoint_raises_for_archived_vcf_quad(self):
        quad = create_fake_quad(self.user, self.gb)
        self._archive(quad.cohort.vcf)
        cohort = type(quad.cohort).objects.get(pk=quad.cohort.pk)
        with self.assertRaises(DataArchivedError):
            _ = cohort.cohort_genotype_collection

    def test_cohort_mixin_returns_none_when_archived(self):
        """ Analysis-side accessor catches DataArchivedError and returns None
            so node-internal 'no source' code paths still work. """
        cohort = create_fake_cohort(self.user, self.gb)
        self._archive(cohort.vcf)

        class _StubNode(CohortMixin):
            def __init__(self, cohort):
                self._cohort = cohort

            def _get_cohort(self):
                return self._cohort

        node = _StubNode(type(cohort).objects.get(pk=cohort.pk))
        self.assertIsNone(node.cohort_genotype_collection)
