from django.contrib.auth.models import User
from django.test import TestCase
from django.test.client import Client
from django.urls.base import reverse

from analysis.models.enums import ClassificationsNodeInput
from analysis.models.nodes.filters.damage_node import DamageNode
from analysis.models.nodes.filters.intersection_node import IntersectionNode
from analysis.models.nodes.sources.classifications_node import ClassificationsNode
from analysis.tests.utils import AnalysisSetupMixin
from library.genomics.vcf_enums import VariantClass
from library.guardian_utils import assign_permission_to_user_and_groups


class NodeEditorRenderTest(AnalysisSetupMixin, TestCase):
    """ The editors for the nodes carrying the new controls - a template that can't render the
        form (a renamed field, a form helper that doesn't return what the loop expects) is a 500 """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get(username=f"test_{cls.__name__}")
        assign_permission_to_user_and_groups(cls.user, cls.analysis)

    def _get_editor(self, node):
        client = Client()
        client.force_login(self.user)
        url = reverse("node_view", kwargs={"analysis_id": self.analysis.pk,
                                           "analysis_version": self.analysis.version,
                                           "node_id": node.pk,
                                           "node_version": node.version,
                                           "extra_filters": "default"})
        response = client.get(url)
        self.assertEqual(200, response.status_code)
        return response

    def test_damage_node_variant_class_groups(self):
        node = DamageNode.objects.create(analysis=self.analysis, variant_class=[VariantClass.SNV])
        response = self._get_editor(node)
        content = response.content.decode()
        for group_name in ("SNV", "Indel", "Copy number", "Rearrangement", "Fusion", "Other"):
            self.assertIn(group_name, content)
        # The checkbox for the node's stored value comes back ticked
        self.assertRegex(content, r'<input checked[^>]*id="id_variant_class_0"[^>]*'
                                  rf'value="{VariantClass.SNV.value}"')

    def test_intersection_node_variant_text(self):
        node = IntersectionNode.objects.create(analysis=self.analysis,
                                               accordion_panel=IntersectionNode.VARIANTS,
                                               variant_text="1:100-200", variant_regions=["1:100-200"])
        content = self._get_editor(node).content.decode()
        self.assertIn("1:100-200", content)

    def test_classifications_node_clinvar_section(self):
        node = ClassificationsNode.objects.create(analysis=self.analysis,
                                                  node_input=ClassificationsNodeInput.PARENT_NOT_MATCHING,
                                                  clinvar_benign=True, clinvar_variation_ids=[12345])
        content = self._get_editor(node).content.decode()
        self.assertIn("id_clinvar_benign", content)
        self.assertIn("12345", content)
