from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.urls import reverse

from analysis.models import AllVariantsNode, Analysis
from annotation.fake_annotation import get_fake_annotation_version
from snpdb.models import GenomeBuild


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class HorizontalModeTestCase(TestCase):
    """ The horizontal layout (full width canvas, grid docked along the bottom) is rendered into the
        page, so both orientations need to keep rendering - see analysis.html """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_superuser(f"test_user_{__file__}")
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)

    def _get_analysis_page(self, horizontal_mode: bool) -> str:
        analysis = Analysis(genome_build=self.grch37, name="horizontal mode test")
        analysis.set_defaults_and_save(self.user)
        analysis.analysis_horizontal_mode = horizontal_mode
        analysis.save()
        AllVariantsNode.objects.create(analysis=analysis, x=10, y=20)

        self.client.force_login(self.user)
        response = self.client.get(reverse("analysis", kwargs={"analysis_id": analysis.pk}))
        self.assertEqual(response.status_code, 200)

        editor_and_grid = self.client.get(reverse("analysis_editor_and_grid", kwargs={"analysis_id": analysis.pk}))
        self.assertEqual(editor_and_grid.status_code, 200)
        return response.content.decode() + editor_and_grid.content.decode()

    def test_horizontal_mode(self):
        page = self._get_analysis_page(True)
        self.assertIn("ANALYSIS_HORIZONTAL_MODE = true", page)
        self.assertIn("node-editor-drawer", page)
        self.assertIn("bottom-pane-tabs", page)

    def test_vertical_mode(self):
        page = self._get_analysis_page(False)
        self.assertIn("ANALYSIS_HORIZONTAL_MODE = false", page)
        self.assertNotIn("node-editor-drawer", page)
        self.assertNotIn("bottom-pane-tabs", page)
