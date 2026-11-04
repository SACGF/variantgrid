"""
ACMG criterion presets (#2605) fill a node editor's form client side, so what the server owns is the preset data:
values that are legal for the form's fields and complete (every field reset, then the criterion's own cutoffs), which
only the deployment's data (gnomAD gene constraint, the Population preset rows in admin) and the annotation version
(which predictors) vary.
"""
from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.test.client import Client
from django.urls.base import reverse

from analysis.enums import GroupOperation, NodeMatchInput
from analysis.forms.forms_nodes import ClinVarNodeForm, DamageNodeForm, PopulationNodeForm
from analysis.models import Analysis
from analysis.models.nodes.filters.clinvar_node import ClinVarNode
from analysis.models.nodes.filters.damage_node import DamageNode
from analysis.models.nodes.filters.population_node import PopulationNode, PopulationNodePreset
from annotation.enums import PathogenicityImpact
from annotation.fake_data import get_fake_annotation_settings_dict, get_fake_annotation_version
from annotation.models import VariantAnnotationVersion
from annotation.pathogenicity_predictions import (
    GENE_LOEUF_CUTOFF,
    GENE_MIS_Z_CUTOFF,
    TOOLS_BY_RAW_FIELD,
)
from genes.models import GeneSymbol, GnomADGeneConstraint
from library.guardian_utils import assign_permission_to_user_and_groups
from snpdb.enums import AlleleOriginFilterDefault
from snpdb.models import CachedWebResource
from snpdb.models.models_genome import GenomeBuild


def _post_data(values: dict) -> dict:
    """ A preset's values as the browser would post the filled form """
    data = {}
    for name, value in values.items():
        if value is True:
            data[name] = "on"
        elif value is False or value is None:
            continue
        elif isinstance(value, list):
            data[name] = value
        else:
            data[name] = str(value)
    return data


@override_settings(**get_fake_annotation_settings_dict(columns_version=4))
class NodePresetsTest(TestCase):
    """ GRCh38 / columns_version 5 / VEP 116 scores every predictor, and the deployment has gene constraint """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        genome_build = GenomeBuild.get_name_or_alias("GRCh38")
        get_fake_annotation_version(genome_build)
        vav = VariantAnnotationVersion.objects.filter(genome_build=genome_build).first()
        vav.columns_version = 5
        vav.vep = 116
        vav.save()
        cls.user = User.objects.get_or_create(username=f"test_{cls.__name__}")[0]
        cls.analysis = Analysis(genome_build=genome_build)
        cls.analysis.set_defaults_and_save(cls.user)
        assign_permission_to_user_and_groups(cls.user, cls.analysis)
        cached_web_resource = CachedWebResource.objects.create(name="test_gnomad_gene_constraint")
        GnomADGeneConstraint.objects.create(gene_symbol=GeneSymbol.objects.get_or_create(symbol="GENE1")[0],
                                            cached_web_resource=cached_web_resource,
                                            mane_select=True, lof_oe_ci_upper=0.2, mis_z_score=4.0)

    def _presets(self, form_class, node) -> dict[str, dict]:
        return {p["key"]: p for p in form_class(instance=node).get_preset_data()}

    def _node_url(self, node) -> str:
        return reverse("node_view", kwargs={"analysis_id": self.analysis.pk,
                                            "analysis_version": self.analysis.version,
                                            "node_id": node.pk, "node_version": node.version,
                                            "extra_filters": "default"})

    def _save_preset(self, form_class, node, preset: dict):
        """ The browser posts the whole form, kept fields included """
        kept = {name: getattr(node, name) for name in form_class.PRESET_KEEP_FIELDS}
        form = form_class(data=_post_data(preset["values"] | kept), instance=node)
        self.assertTrue(form.is_valid(), f"{preset['key']}: {form.errors}")
        return form.save(commit=False)

    def test_effect_presets(self):
        presets = self._presets(DamageNodeForm, DamageNode(analysis=self.analysis))
        self.assertEqual(["pvs1", "pp3", "pp2"], list(presets))

        pvs1 = self._save_preset(DamageNodeForm, DamageNode(analysis=self.analysis), presets["pvs1"])
        self.assertEqual((PathogenicityImpact.HIGH, True), (pvs1.impact_min, pvs1.impact_required))
        self.assertEqual((GENE_LOEUF_CUTOFF, True, True),
                         (pvs1.gene_loeuf_max, pvs1.gene_loeuf_required, pvs1.gene_loeuf_allow_null))
        self.assertIsNone(pvs1.ptc_nmd)

        pp3 = self._save_preset(DamageNodeForm, DamageNode(analysis=self.analysis), presets["pp3"])
        self.assertEqual((PathogenicityImpact.MODERATE, True), (pp3.impact_min, pp3.impact_required))
        self.assertEqual(TOOLS_BY_RAW_FIELD["revel_score"].raw_pathogenic_threshold, pp3.revel_score_min)
        self.assertFalse(pp3.revel_score_required)
        # Author cutoffs and uncalibrated scores aren't PP3 evidence
        self.assertIsNone(pp3.clinpred_score_min)
        self.assertIsNone(pp3.popeve_score_max)
        self.assertIsNone(pp3.gene_loeuf_max)

        pp2 = self._save_preset(DamageNodeForm, DamageNode(analysis=self.analysis), presets["pp2"])
        self.assertEqual((GENE_MIS_Z_CUTOFF, True, False),
                         (pp2.gene_mis_z_min, pp2.gene_mis_z_required, pp2.gene_mis_z_allow_null))

    def test_effect_preset_resets_what_the_node_had(self):
        node = DamageNode(analysis=self.analysis, revel_score_min=0.9, splice_required=True)
        values = self._presets(DamageNodeForm, node)["pvs1"]["values"]
        self.assertIsNone(values["revel_score_min"])
        self.assertFalse(values["splice_required"])
        self.assertEqual([], values["variant_class"])

    def test_effect_presets_without_gene_constraint(self):
        GnomADGeneConstraint.objects.all().delete()
        presets = self._presets(DamageNodeForm, DamageNode(analysis=self.analysis))
        self.assertEqual(["pvs1", "pp3"], list(presets))
        self.assertNotIn("gene_loeuf_max", {k for k, v in presets["pvs1"]["values"].items() if v is not None})

    def test_population_presets_are_admin_rows(self):
        PopulationNodePreset.objects.all().delete()
        PopulationNodePreset.objects.create(name="disabled", enabled=False)
        row = PopulationNodePreset.objects.create(name="PM2 - X-linked", percent=0.05, gnomad_popmax_af=True,
                                                  gnomad_fafmax_faf95_max=True, gnomad_hemi_count_max=1)
        key = f"population_{row.pk}"
        presets = self._presets(PopulationNodeForm, PopulationNode(analysis=self.analysis))
        self.assertEqual([key], list(presets))
        self.assertEqual("PM2 - X-linked (≤ 0.05%, hemi ≤ 1)", presets[key]["label"])
        self.assertIn("gnomAD hemizygous ≤ 1", presets[key]["description"])
        self.assertEqual([], presets[key]["values"]["gnomad_population"])

        node = PopulationNode.objects.create(analysis=self.analysis, af_1kg=True, use_internal_counts=True,
                                             gnomad_hom_alt_max=3)
        pm2 = self._save_preset(PopulationNodeForm, node, presets[key])
        self.assertEqual(0.05, pm2.percent)
        self.assertEqual(GroupOperation.ANY, pm2.group_operation)
        self.assertEqual((False, True), (pm2.gnomad_af, pm2.gnomad_popmax_af))
        # Only drawn (and so only applied) where the annotation version has gnomAD filtering AF
        self.assertEqual(node.has_filtering_allele_frequency, pm2.gnomad_fafmax_faf95_max)
        self.assertEqual((1, None), (pm2.gnomad_hemi_count_max, pm2.gnomad_hom_alt_max))
        self.assertTrue(pm2.keep_internally_classified_pathogenic)
        self.assertEqual((False, False), (pm2.af_1kg, pm2.use_internal_counts))

    def test_population_editor_links_admins_to_presets(self):
        admin_url = reverse("admin:analysis_populationnodepreset_changelist")
        superuser = User.objects.create_superuser(username="test_node_presets_admin")
        assign_permission_to_user_and_groups(superuser, self.analysis)
        node = PopulationNode.objects.create(analysis=self.analysis)
        for user, shows_link in [(self.user, False), (superuser, True)]:
            client = Client()
            client.force_login(user)
            html = client.get(self._node_url(node)).content.decode()
            self.assertEqual(shows_link, admin_url in html)

    def test_clinvar_pp5(self):
        node = ClinVarNode(analysis=self.analysis, node_input=NodeMatchInput.MATCHING_VARIANTS, stars_min=0)
        presets = self._presets(ClinVarNodeForm, node)
        self.assertNotIn("node_input", presets["pp5"]["values"])
        pp5 = self._save_preset(ClinVarNodeForm, node, presets["pp5"])
        self.assertEqual(AlleleOriginFilterDefault.GERMLINE, pp5.allele_origin)
        self.assertEqual(2, pp5.stars_min)
        self.assertEqual([True, True, False, False, False, False],
                         [pp5.germline_pathogenic, pp5.germline_likely_pathogenic, pp5.germline_uncertain,
                          pp5.germline_likely_benign, pp5.germline_benign, pp5.germline_other])
        self.assertFalse(pp5.conflicting)
        self.assertEqual(NodeMatchInput.MATCHING_VARIANTS, pp5.node_input)

    def test_editors_offer_the_dropdown(self):
        client = Client()
        client.force_login(self.user)
        for node_class, label in [(DamageNode, "PVS1"), (PopulationNode, "PM2 - recessive"), (ClinVarNode, "PP5")]:
            node = node_class.objects.create(analysis=self.analysis)
            html = client.get(self._node_url(node)).content.decode()
            self.assertIn(f'id="node-preset-{node.pk}"', html)
            self.assertIn(label, html)
