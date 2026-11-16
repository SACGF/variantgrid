""" The splicing filter's "keep if not scored" (#2836): a variant passes as unscored only with none of the splicing
    scores the filter reads. Nodes saved under the earlier rule (any one score missing passes - and dbscSNV only
    scores splice-consensus SNVs, so nearly every variant did) keep it through splice_allow_null_any_missing """
import json
import tempfile
from importlib import import_module

from django.apps import apps
from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.analysis_import_export import analysis_export_to_dict, analysis_import
from analysis.models.nodes.filters.damage_node import DamageNode
from analysis.tests.utils import AnalysisSetupMixin, create_annotated_variants, node_keeps

REQUIRED_ALLOW_NULL = {"splice_min": 0.5, "splice_required": True, "splice_allow_null": True}


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class SpliceAllowNullTest(AnalysisSetupMixin, TestCase):
    ROWS = {
        "spliceai_high": {"spliceai_max_ds": 0.6, "dbscsnv_ada_score": 0.1, "dbscsnv_rf_score": 0.1},
        "spliceai_low_no_dbscsnv": {"spliceai_max_ds": 0.1},
        "no_splicing_score": {},
    }

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        vav = cls.analysis.annotation_version.variant_annotation_version
        vav.backfilled_spliceai_max_ds = True
        vav.save()
        cls.variants = create_annotated_variants(vav, cls.ROWS)

    def _keeps(self, **kwargs) -> set[str]:
        return node_keeps(DamageNode(analysis=self.analysis, **kwargs), self.variants)

    def test_corrected_rule_keeps_only_variants_with_no_score(self):
        self.assertEqual({"spliceai_high", "no_splicing_score"}, self._keeps(**REQUIRED_ALLOW_NULL))

    def test_earlier_rule_keeps_a_variant_missing_any_score(self):
        self.assertEqual(set(self.ROWS), self._keeps(splice_allow_null_any_missing=True, **REQUIRED_ALLOW_NULL))

    def test_rules_agree_without_allow_null(self):
        kwargs = {"splice_min": 0.5, "splice_required": True, "splice_allow_null": False}
        corrected = DamageNode(analysis=self.analysis, **kwargs)
        earlier = DamageNode(analysis=self.analysis, splice_allow_null_any_missing=True, **kwargs)
        self.assertEqual(str(earlier._get_node_q()), str(corrected._get_node_q()))
        self.assertEqual({"spliceai_high"}, node_keeps(corrected, self.variants))

    def test_description_words_each_rule(self):
        corrected = DamageNode(analysis=self.analysis, **REQUIRED_ALLOW_NULL)
        earlier = DamageNode(analysis=self.analysis, splice_allow_null_any_missing=True, **REQUIRED_ALLOW_NULL)
        self.assertEqual("Keeps variants with splicing ≥ 0.5 (or no splicing score).", corrected.get_node_description())
        self.assertEqual("Keeps variants with splicing ≥ 0.5 (or any splicing score missing).",
                         earlier.get_node_description())

    def test_migration_marks_nodes_where_the_earlier_rule_applied(self):
        nodes = {
            "earlier_rule": DamageNode.objects.create(analysis=self.analysis, **REQUIRED_ALLOW_NULL),
            "no_allow_null": DamageNode.objects.create(analysis=self.analysis, splice_min=0.5, splice_required=True,
                                                       splice_allow_null=False),
            "not_required": DamageNode.objects.create(analysis=self.analysis, splice_min=0.5),
            "no_splice_min": DamageNode.objects.create(analysis=self.analysis, alphagenome_splicing_min=1.0,
                                                       splice_required=True, splice_allow_null=True),
        }
        migration = import_module("analysis.migrations.0021_damagenode_splice_allow_null_any_missing")
        migration._keep_earlier_splice_allow_null(apps, None)
        marked = {name for name, node in nodes.items()
                  if DamageNode.objects.get(pk=node.pk).splice_allow_null_any_missing}
        self.assertEqual({"earlier_rule"}, marked)

    def test_export_from_before_the_field_imports_on_the_earlier_rule(self):
        DamageNode.objects.create(analysis=self.analysis, **REQUIRED_ALLOW_NULL)
        analysis_json = analysis_export_to_dict(self.analysis)
        for node_record in analysis_json["nodes"]:
            node_record["fields"].pop("splice_allow_null_any_missing", None)
        with tempfile.NamedTemporaryFile("w", suffix=".json") as f:
            json.dump(analysis_json, f)
            f.flush()
            imported = analysis_import(User.objects.get(username=f"test_{self.__class__.__name__}"), None, f.name)
        node = DamageNode.objects.get(analysis=imported)
        self.assertTrue(node.splice_allow_null_any_missing)
