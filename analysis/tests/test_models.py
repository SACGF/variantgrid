import unittest

from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.utils import timezone

from analysis.models import Analysis, AnalysisLock
from annotation.fake_annotation import get_fake_annotation_version
from snpdb.models import GenomeBuild


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
        analysis.is_locked()

        AnalysisLock.objects.create(analysis=analysis, locked=True, user=self.owner_user, date=timezone.now())
        # Bump version to expire cache
        analysis.version += 1
        analysis.save()
        self.assertTrue(analysis.is_locked())

        # Nobody should be able to write if locked
        self.assertFalse(analysis.can_write(self.owner_user))
        self.assertFalse(analysis.can_write(self.admin_user))
        self.assertFalse(analysis.can_write(self.non_owner_user))


if __name__ == '__main__':
    unittest.main()
