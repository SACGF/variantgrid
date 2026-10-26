from unittest import mock

from django.contrib.auth.models import User
from django.test import TestCase

from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.models import GenomeBuild, UserSettings


class GenomeBuildManagerTest(TestCase):

    def test_cached_build_not_returned_for_another_user(self):
        """ The cached build is the current user's - asking for another user's resolves theirs """
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        grch38 = GenomeBuild.get_name_or_alias("GRCh38")
        other_user = User.objects.get_or_create(username='genome_build_manager_other_user')[0]

        GenomeBuildManager.set_current_genome_build(grch38)
        self.addCleanup(GenomeBuildManager.set_current_genome_build, None)

        with mock.patch.object(UserSettings, "get_genome_build_or_default", return_value=grch37):
            self.assertEqual(grch37, GenomeBuildManager.get_current_genome_build(other_user))
        self.assertEqual(grch38, GenomeBuildManager.get_current_genome_build(), "Cache left as it was")
