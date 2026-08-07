from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse

from analysis.models import Analysis, AllVariantsNode
from analysis.models.nodes.analysis_node import NodeCount
from annotation.fake_annotation import get_fake_annotation_version
from snpdb.models import GenomeBuild
from snpdb.models.models_enums import BuiltInFilters
from snpdb.models.models_user_settings import GlobalSettings


class NodeGridAutoLoadViewTest(TestCase):
    """ Large nodes (count >= threshold) defer their grid behind a placeholder;
        small nodes auto-load exactly as before. """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username="test_node_grid_auto_load")[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(cls.user)
        # The threshold is read from the cascaded settings (Global is the base), so set it there.
        global_settings = GlobalSettings.objects.get()
        global_settings.node_grid_auto_load_max_variants = 100
        global_settings.save()

    def setUp(self):
        self.client.force_login(self.user)

    def _get_grid(self, node, extra_filters="default"):
        kwargs = {
            "analysis_id": self.analysis.pk,
            "analysis_version": self.analysis.version,
            "node_id": node.pk,
            "node_version": node.version,
            "extra_filters": extra_filters,
        }
        return self.client.get(reverse("node_data_grid", kwargs=kwargs))

    def test_large_node_defers_grid(self):
        node = AllVariantsNode.objects.create(analysis=self.analysis, count=500)
        response = self._get_grid(node)
        self.assertEqual(response.status_code, 200)
        self.assertFalse(response.context["grid_auto_load"])
        # The placeholder block (with its Load/CSV/VCF icons) only renders for deferred nodes
        self.assertContains(response, f'id="grid-placeholder-{node.pk}"')

    def test_small_node_auto_loads(self):
        node = AllVariantsNode.objects.create(analysis=self.analysis, count=10)
        response = self._get_grid(node)
        self.assertEqual(response.status_code, 200)
        self.assertTrue(response.context["grid_auto_load"])
        self.assertNotContains(response, f'id="grid-placeholder-{node.pk}"')

    def test_uncomputed_count_defers_grid(self):
        node = AllVariantsNode.objects.create(analysis=self.analysis, count=None)
        response = self._get_grid(node)
        self.assertEqual(response.status_code, 200)
        self.assertFalse(response.context["grid_auto_load"])
        self.assertContains(response, f'id="grid-placeholder-{node.pk}"')

    def test_extra_filter_uses_filtered_count(self):
        """ A large node deferred on its total count auto-loads when viewing an extra filter
            (eg ClinVar) whose filtered count is below the threshold. """
        node = AllVariantsNode.objects.create(analysis=self.analysis, count=500)
        NodeCount.objects.create(node_version=node.node_version, label=BuiltInFilters.CLINVAR, count=10)

        # Default (unfiltered) view still defers - 500 >= 100
        default_response = self._get_grid(node)
        self.assertFalse(default_response.context["grid_auto_load"])

        # Filtered ClinVar view auto-loads - the filtered count (10) is below the threshold
        clinvar_response = self._get_grid(node, extra_filters=BuiltInFilters.CLINVAR)
        self.assertEqual(clinvar_response.status_code, 200)
        self.assertEqual(clinvar_response.context["grid_row_count"], 10)
        self.assertTrue(clinvar_response.context["grid_auto_load"])
        self.assertNotContains(clinvar_response, f'id="grid-placeholder-{node.pk}"')
