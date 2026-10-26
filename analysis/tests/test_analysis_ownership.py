""" Taking over the analysis of a deactivated owner (analysis/views/views_analysis_settings.py:analysis_take_ownership) """
from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse
from guardian.shortcuts import assign_perm

from analysis.models import Analysis
from annotation.fake_data import get_fake_annotation_version
from library.django_utils.unittest_utils import prevent_request_warnings
from library.guardian_utils import DjangoPermission
from snpdb.models import GenomeBuild


@patch("analysis.views.views_analysis_settings.node_utils.reload_analysis_nodes")
class TestAnalysisTakeOwnership(TestCase):

    @classmethod
    def setUpTestData(cls):
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(grch37)
        cls.owner = User.objects.create_user(username="ownership_owner")
        cls.editor = User.objects.create_user(username="ownership_editor")
        cls.viewer = User.objects.create_user(username="ownership_viewer")

        cls.analysis = Analysis(genome_build=grch37)
        cls.analysis.set_defaults_and_save(cls.owner)
        for permission in (DjangoPermission.READ, DjangoPermission.WRITE):
            assign_perm(DjangoPermission.perm(cls.analysis, permission), cls.editor, cls.analysis)
        assign_perm(DjangoPermission.perm(cls.analysis, DjangoPermission.READ), cls.viewer, cls.analysis)
        cls.url = reverse("analysis_take_ownership", kwargs={"analysis_id": cls.analysis.pk})

    def _post_as(self, user):
        self.client.force_login(user)
        return self.client.post(self.url)

    @prevent_request_warnings
    def test_only_from_a_deactivated_owner(self, mock_reload):
        self.assertEqual(self._post_as(self.editor).status_code, 403)

        User.objects.filter(pk=self.owner.pk).update(is_active=False)
        self.assertEqual(self._post_as(self.viewer).status_code, 403)  # needs write
        self.assertEqual(self._post_as(self.editor).status_code, 302)

        self.analysis.refresh_from_db()
        self.assertEqual(self.analysis.user, self.editor)
        mock_reload.assert_called_once_with(self.analysis.pk)
