"""
DamageNode raw-score sliders are driven by TOOLS filtered through the annotation version's
visible_columns (#1808), so a tool only gets a slider where its column is actually populated.
The VEP 116 plugin tools are the first to need this: ProtVar is columns_version 5 on any build,
while EVE / popEVE / PromoterAI are GRCh38 + VEP >= 116 as well.

They also cover the three ways a raw score can run - popEVE is a log-likelihood ratio (lower is
worse) and PromoterAI is signed (magnitude is what counts), unlike every dbNSFP score.
"""
from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.test.client import Client
from django.urls.base import reverse

from analysis.forms.forms_nodes import DamageNodeForm
from analysis.models import Analysis
from analysis.models.nodes.filters.damage_node import DamageNode
from annotation.fake_annotation import (
    get_fake_annotation_settings_dict,
    get_fake_annotation_version,
)
from annotation.models import VariantAnnotationVersion
from annotation.pathogenicity_predictions import (
    TOOLS,
    RawScoreDirection,
    raw_score_pathogenic_funcs,
)
from library.guardian_utils import assign_permission_to_user_and_groups
from snpdb.models.models_genome import GenomeBuild

PROTVAR = {"protvar_stability"}  # columns_version 5, every build
GRCH38_VEP116 = {"eve_score", "popeve_score", "promoter_ai_score"}
DBNSFP_RAW_FIELDS = {t.raw_field for t in TOOLS if t.raw_field} - PROTVAR - GRCH38_VEP116


class _RawScoreSliderMixin:
    """ columns_version / vep are read off the VAV, so forcing it is the only way to move the node
        between the versions that populate the plugin columns and those that don't. """
    BUILD_NAME = "GRCh38"
    COLUMNS_VERSION = 5
    VEP_VERSION = 116

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        genome_build = GenomeBuild.get_name_or_alias(cls.BUILD_NAME)
        get_fake_annotation_version(genome_build)
        vav = VariantAnnotationVersion.objects.filter(genome_build=genome_build).first()
        vav.columns_version = cls.COLUMNS_VERSION
        vav.vep = cls.VEP_VERSION
        vav.save()
        cls.user = User.objects.get_or_create(username=f"test_{cls.__name__}")[0]
        cls.analysis = Analysis(genome_build=genome_build)
        cls.analysis.set_defaults_and_save(cls.user)
        assign_permission_to_user_and_groups(cls.user, cls.analysis)

    def _raw_fields(self) -> set[str]:
        return {t.raw_field for t in DamageNode(analysis=self.analysis).get_raw_score_tools()}

    def _pred_fields(self) -> set[str]:
        return {t.pred_field for t in DamageNode(analysis=self.analysis).get_pred_tools()}

    def _editor_html(self) -> str:
        node = DamageNode.objects.create(analysis=self.analysis)
        client = Client()
        client.force_login(self.user)
        url = reverse("node_view", kwargs={"analysis_id": self.analysis.pk,
                                           "analysis_version": self.analysis.version,
                                           "node_id": node.pk,
                                           "node_version": node.version,
                                           "extra_filters": "default"})
        response = client.get(url)
        self.assertEqual(200, response.status_code)
        return response.content.decode()


@override_settings(**get_fake_annotation_settings_dict(columns_version=4))
class TestPluginSlidersOffered(_RawScoreSliderMixin, TestCase):
    """ GRCh38 / columns_version 5 / VEP 116 - the only combination populating all four """

    def test_all_tools_offered(self):
        self.assertEqual(self._raw_fields(), DBNSFP_RAW_FIELDS | PROTVAR | GRCH38_VEP116)

    def test_eve_class_offered_as_prediction_dropdown(self):
        self.assertIn("eve_class", self._pred_fields())

    def test_eve_min_filters(self):
        node = DamageNode(analysis=self.analysis, eve_score_min=0.7)
        self.assertTrue(node.modifies_parents())
        self.assertTrue(node.has_individual_pathogenic_predictions())
        self.assertIn("variantannotation__eve_score__gte", str(node._get_node_q()))

    def test_eve_required_allows_null(self):
        node = DamageNode(analysis=self.analysis, eve_score_min=0.7, eve_score_required=True,
                          eve_score_allow_null=True)
        self.assertTrue(node.has_required())
        self.assertIn("variantannotation__eve_score__isnull", str(node._get_node_q()))

    def test_eve_class_filters_on_the_word(self):
        node = DamageNode(analysis=self.analysis, eve_class="Pathogenic")
        self.assertTrue(node.modifies_parents())
        q_str = str(node._get_node_q())
        self.assertIn("variantannotation__eve_class", q_str)
        self.assertIn("Pathogenic", q_str)

    def test_popeve_keeps_scores_at_or_below(self):
        """ popEVE runs the other way to every dbNSFP score - more negative is more damaging """
        node = DamageNode(analysis=self.analysis, popeve_score_max=-3.0)
        q_str = str(node._get_node_q())
        self.assertIn("variantannotation__popeve_score__lte", q_str)
        self.assertNotIn("variantannotation__popeve_score__gte", q_str)

    def test_promoter_ai_matches_either_sign(self):
        """ Signed score - a strong drop in expression is as interesting as a strong rise """
        node = DamageNode(analysis=self.analysis, promoter_ai_score_min=0.5)
        q_str = str(node._get_node_q())
        self.assertIn("variantannotation__promoter_ai_score__gte", q_str)
        self.assertIn("variantannotation__promoter_ai_score__lte", q_str)
        self.assertIn("-0.5", q_str)

    def test_protvar_stability_keeps_destabilising(self):
        node = DamageNode(analysis=self.analysis, protvar_stability_min=2.0)
        self.assertIn("variantannotation__protvar_stability__gte", str(node._get_node_q()))

    def test_plugin_tools_stay_out_of_the_damage_count(self):
        """ None are ClinGen-calibrated, so predictions_num_pathogenic is unchanged - no re-annotation """
        counted = set(raw_score_pathogenic_funcs())
        self.assertEqual(set(), counted & (PROTVAR | GRCH38_VEP116))

    def test_form_row_carries_slider_bounds(self):
        """ Uncalibrated, so bounds but no colour-band data attrs """
        form = DamageNodeForm(instance=DamageNode(analysis=self.analysis))
        row = next(r for r in form.get_raw_score_rows() if r["field"] == "eve_score")
        attrs = row["min_field"].field.widget.attrs
        self.assertEqual(attrs["min"], 0.0)
        self.assertEqual(attrs["max"], 1.0)
        self.assertNotIn("data-pathogenic-min", attrs)
        self.assertNotIn("data-benign-max", attrs)

    def test_inverted_and_signed_tools_are_labelled(self):
        """ 'Popeve score max' / 'Promoter ai score min' would not say which side is kept """
        form = DamageNodeForm(instance=DamageNode(analysis=self.analysis))
        labels = {r["field"]: str(r["min_field"].label) for r in form.get_raw_score_rows()}
        self.assertEqual("popEVE score ≤", labels["popeve_score"])
        self.assertEqual("PromoterAI |score| ≥", labels["promoter_ai_score"])

    def test_editor_renders_plugin_sliders(self):
        html = self._editor_html()
        for raw_field in GRCH38_VEP116 | PROTVAR:
            self.assertIn(f'id="{raw_field}-slider"', html)
        self.assertIn('name="popeve_score_max"', html)
        self.assertIn('name="eve_class"', html)


@override_settings(**get_fake_annotation_settings_dict(columns_version=4))
class TestGRCh38OnlyToolsNotOfferedOnGRCh37(_RawScoreSliderMixin, TestCase):
    BUILD_NAME = "GRCh37"

    def test_protvar_offered_but_not_the_grch38_tools(self):
        self.assertEqual(self._raw_fields(), DBNSFP_RAW_FIELDS | PROTVAR)

    def test_eve_class_not_offered(self):
        self.assertNotIn("eve_class", self._pred_fields())

    def test_eve_filter_warns_when_column_not_populated(self):
        """ A node cloned in from a build that has EVE keeps the filter - warn it has no data """
        node = DamageNode(analysis=self.analysis, eve_score_min=0.7)
        self.assertTrue(any("EVE filter has no data" in w for w in node.get_warnings()))

    def test_popeve_filter_warns_when_column_not_populated(self):
        node = DamageNode(analysis=self.analysis, popeve_score_max=-3.0)
        self.assertTrue(any("popEVE filter has no data" in w for w in node.get_warnings()))

    def test_dbnsfp_filter_does_not_warn(self):
        node = DamageNode(analysis=self.analysis, alphamissense_score_min=0.5)
        self.assertEqual([], [w for w in node.get_warnings() if "has no data" in w])

    def test_editor_renders_dbnsfp_sliders_without_the_grch38_tools(self):
        html = self._editor_html()
        for raw_field in GRCH38_VEP116:
            self.assertNotIn(f'id="{raw_field}-slider"', html)
        for raw_field in DBNSFP_RAW_FIELDS | PROTVAR:
            self.assertIn(f'id="{raw_field}-slider"', html)


@override_settings(**get_fake_annotation_settings_dict(columns_version=4))
class TestGRCh38ToolsNotOfferedBeforeVep116(_RawScoreSliderMixin, TestCase):
    VEP_VERSION = 112

    def test_protvar_offered_but_not_the_vep116_tools(self):
        self.assertEqual(self._raw_fields(), DBNSFP_RAW_FIELDS | PROTVAR)


@override_settings(**get_fake_annotation_settings_dict(columns_version=4))
class TestPluginToolsNotOfferedAtColumnsVersion4(_RawScoreSliderMixin, TestCase):
    COLUMNS_VERSION = 4

    def test_dbnsfp_tools_only(self):
        self.assertEqual(self._raw_fields(), DBNSFP_RAW_FIELDS)


class TestRawScoreDirection(TestCase):
    """ The model field is named for the side it keeps, so a reader of the node can't mistake it """

    def test_threshold_field_named_for_direction(self):
        by_name = {t.name: t for t in TOOLS}
        self.assertEqual("alphamissense_score_min", by_name["AlphaMissense"].node_threshold_field)
        self.assertEqual("popeve_score_max", by_name["popEVE"].node_threshold_field)
        self.assertEqual("promoter_ai_score_min", by_name["PromoterAI"].node_threshold_field)

    def test_every_tool_has_damage_node_fields(self):
        """ A tool added to TOOLS without model fields would silently drop out of the editor """
        for tool in TOOLS:
            if tool.raw_field:
                for field in (tool.node_threshold_field, f"{tool.raw_field}_required",
                              f"{tool.raw_field}_allow_null"):
                    DamageNode._meta.get_field(field)
            if tool.pred_field:
                for field in (tool.pred_field, f"{tool.pred_field}_required",
                              f"{tool.pred_field}_allow_null"):
                    DamageNode._meta.get_field(field)

    def test_only_higher_is_worse_tools_can_be_counted(self):
        """ predictions_num_pathogenic applies score >= threshold, so an inverted or signed
            score must not carry a PP3 cutoff """
        counted = set(raw_score_pathogenic_funcs())
        for tool in TOOLS:
            if tool.raw_direction != RawScoreDirection.HIGHER:
                self.assertNotIn(tool.raw_field, counted,
                                 f"{tool.name} does not run higher=worse, so it can't use a PP3 cutoff")
