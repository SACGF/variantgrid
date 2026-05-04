"""
Tests for DataArchiveMixin and the chokepoint behaviour added in #1536.

@see claude/issue_1536_data_archive_plan.md
"""

from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase
from django.utils import timezone

from annotation.models import VariantAnnotationVersion
from genes.models import GeneCoverageCollection
from snpdb.archive import archive_vcf, ArchivePreconditionError, DataArchivedError
from snpdb.models import GenomeBuild, VCF
from snpdb.tests.utils.fake_cohort_data import create_fake_cohort


class DataArchiveMixinTests(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.create_user(username="archive_user")
        cls.genome_build = GenomeBuild.grch38()
        cls.cohort = create_fake_cohort(cls.user, cls.genome_build)
        cls.vcf: VCF = cls.cohort.vcf

    def test_mixin_fields_present(self):
        """ Mixin fields are exposed on VCF, VAV, GCC. """
        for field in ("data_archived_date", "data_archived_by", "data_archive_reason", "data_restorable_from"):
            self.assertTrue(hasattr(self.vcf, field), f"VCF missing {field}")
            self.assertIn(field, [f.name for f in VariantAnnotationVersion._meta.get_fields()])
            self.assertIn(field, [f.name for f in GeneCoverageCollection._meta.get_fields()])

    def test_data_archived_property_reflects_date(self):
        self.assertFalse(self.vcf.data_archived)
        self.vcf.data_archived_date = timezone.now()
        self.assertTrue(self.vcf.data_archived)

    def test_chokepoint_raises_on_archived_vcf(self):
        """ Cohort.cohort_genotype_collection raises DataArchivedError when its VCF is archived. """
        self.vcf.data_archived_date = timezone.now()
        self.vcf.data_archived_by = self.user
        self.vcf.save()
        # cached_property: invalidate by re-fetching cohort
        cohort = type(self.cohort).objects.get(pk=self.cohort.pk)
        with self.assertRaises(DataArchivedError):
            _ = cohort.cohort_genotype_collection

    def test_chokepoint_passes_when_not_archived(self):
        cohort = type(self.cohort).objects.get(pk=self.cohort.pk)
        self.assertIsNotNone(cohort.cohort_genotype_collection)

    def test_data_archived_error_includes_metadata(self):
        date = timezone.now()
        self.vcf.data_archived_date = date
        self.vcf.data_archived_by = self.user
        self.vcf.data_archive_reason = "Testing archive flow"
        self.vcf.save()
        cohort = type(self.cohort).objects.get(pk=self.cohort.pk)
        try:
            _ = cohort.cohort_genotype_collection
            self.fail("Expected DataArchivedError")
        except DataArchivedError as e:
            msg = str(e)
            self.assertIn("Testing archive flow", msg)
            self.assertIn(date.strftime("%Y-%m-%d"), msg)


class ArchiveVCFHelperTests(TestCase):
    """ archive_vcf preconditions + idempotence. The full delete_internal_data path
        depends on partition tables that aren't created in unit tests, so we patch it. """

    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.create_user(username="archive_helper_user")
        cls.genome_build = GenomeBuild.grch38()
        cls.cohort = create_fake_cohort(cls.user, cls.genome_build)
        cls.vcf: VCF = cls.cohort.vcf

    def test_archive_idempotent(self):
        """ archive_vcf no-ops when already archived. """
        self.vcf.data_archived_date = timezone.now()
        self.vcf.save()
        with patch.object(VCF, "delete_internal_data") as mock_delete:
            archive_vcf(self.vcf, self.user, reason="should noop")
        mock_delete.assert_not_called()

    def test_archive_raises_when_uploaded_file_missing(self):
        """ archive_vcf raises ArchivePreconditionError when no UploadedVCF exists. """
        # The fake VCF has no uploadedvcf record at all → AttributeError caught → raised.
        with self.assertRaises(ArchivePreconditionError):
            archive_vcf(self.vcf, self.user, reason="no source")

    def test_filter_for_user_include_archived_default(self):
        """ Default include_archived=True keeps archived VCFs visible. """
        self.vcf.data_archived_date = timezone.now()
        self.vcf.save()
        qs = VCF.filter_for_user(self.user)
        self.assertIn(self.vcf, qs)

    def test_filter_for_user_excludes_archived_when_requested(self):
        self.vcf.data_archived_date = timezone.now()
        self.vcf.save()
        qs = VCF.filter_for_user(self.user, include_archived=False)
        self.assertNotIn(self.vcf, qs)
