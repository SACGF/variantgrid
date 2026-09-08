from html.parser import HTMLParser

from django.contrib.auth.models import User
from django.test import TestCase
from django.test.client import Client
from django.urls.base import reverse

from analysis.models import VariantTag
from analysis.models.enums import NodeMatchInput
from analysis.models.nodes.filters.classifications_node import ClassificationsNode
from analysis.models.nodes.filters.clinvar_node import ClinVarNode
from analysis.models.nodes.filters.damage_node import DamageNode
from analysis.models.nodes.filters.intersection_node import IntersectionNode
from analysis.models.nodes.filters.tag_node import TagNode
from analysis.tests.utils import AnalysisSetupMixin
from analysis.views.views_node import NODE_DISPATCHER
from annotation.fake_annotation import create_fake_variants
from library.genomics.vcf_enums import VariantClass
from library.guardian_utils import assign_permission_to_user_and_groups
from snpdb.models import Tag, Variant
from snpdb.tests.utils.tag_testing_utils import create_classify_queue_tag


class _FormSubmitDataParser(HTMLParser):
    """ Collects what a browser would POST for a form - the successful controls only """

    def __init__(self, form_id):
        super().__init__()
        self.form_id = form_id
        self.data = {}
        self._in_form = False
        self._select_name = None

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "form":
            self._in_form = attrs.get("id") == self.form_id
        elif not self._in_form:
            return
        elif tag == "input":
            name = attrs.get("name")
            if name and (attrs.get("type") not in ("checkbox", "radio") or "checked" in attrs):
                self.data[name] = attrs.get("value", "")
        elif tag == "select":
            self._select_name = attrs.get("name")
        elif tag == "option" and self._select_name and "selected" in attrs:
            self.data[self._select_name] = attrs.get("value", "")

    def handle_endtag(self, tag):
        if tag == "form":
            self._in_form = False
        elif tag == "select":
            self._select_name = None


def form_submit_data(html, form_id) -> dict:
    parser = _FormSubmitDataParser(form_id)
    parser.feed(html)
    return parser.data


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

    def test_every_node_editor_renders(self):
        # A freshly added node has nothing configured yet - every editor has to cope with that
        for node_class in NODE_DISPATCHER:
            with self.subTest(node_class=node_class.__name__):
                node = node_class.objects.create(analysis=self.analysis)
                self._get_editor(node)

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

    def test_clinvar_node_editor(self):
        node = ClinVarNode.objects.create(analysis=self.analysis,
                                          node_input=NodeMatchInput.PARENT_NOT_MATCHING,
                                          germline_benign=False, variation_ids=[12345])
        content = self._get_editor(node).content.decode()
        self.assertIn("id_germline_benign", content)
        self.assertIn("12345", content)

    def _test_editor_posts_back_valid(self, node, form_id):
        """ The editor's own HTML has to survive a round trip through its form - a widget that submits
            nothing for a required field (a radio group with no option checked) blocks every save """
        response = self._get_editor(node)
        data = form_submit_data(response.content.decode(), form_id)
        data["node_input"] = NodeMatchInput.PARENT_MATCHING

        client = Client()
        client.force_login(self.user)
        post_response = client.post(response.request["PATH_INFO"], data)
        # JSON back means saved - an invalid form comes back as the re-rendered editor HTML
        self.assertEqual({}, post_response.json())
        node.refresh_from_db()
        self.assertEqual(NodeMatchInput.PARENT_MATCHING, node.node_input)

    def test_classifications_node_editor_posts_back_valid(self):
        node = ClassificationsNode.objects.create(analysis=self.analysis)
        self._test_editor_posts_back_valid(node, "classifications-node-form")

    def test_clinvar_node_editor_posts_back_valid(self):
        node = ClinVarNode.objects.create(analysis=self.analysis)
        self._test_editor_posts_back_valid(node, "clinvar-node-form")

    def test_tag_node_editor_to_do_list_is_every_classify_queue_tag(self):
        """ The to-do list is the classify queue vocabulary, so a somatic lab's own queue tag is offered
            for classification the same way the seeded one is @see Tag.classify_queue_qs """
        create_fake_variants(self.grch37)
        queue_tag = create_classify_queue_tag("EditorToDo")
        label_tag = Tag.objects.create(pk="EditorArtefact")
        variant, other_variant = list(Variant.objects.order_by("pk")[:2])
        to_do = VariantTag.objects.create(analysis=self.analysis, genome_build=self.grch37,
                                          variant=variant, tag=queue_tag, user=self.user)
        VariantTag.objects.create(analysis=self.analysis, genome_build=self.grch37,
                                  variant=other_variant, tag=label_tag, user=self.user)

        node = TagNode.objects.create(analysis=self.analysis)
        response = self._get_editor(node)
        self.assertEqual([vt.pk for vt in response.context["requires_classification_tags"]], [to_do.pk])
