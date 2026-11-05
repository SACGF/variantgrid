"""
A count badge is a node configuration (#2662): a BuiltInFilters label counts the node
AnalysisNode.count_label_node configures for it, clicking the badge drops that node, and the dropped node
reuses its parent's stored count rather than counting again.
"""
from django.test import TestCase, override_settings
from django.urls import reverse

from analysis.enums import NodeStatus
from analysis.models import Analysis
from analysis.models.nodes.analysis_node import AnalysisNode
from analysis.models.nodes.filters.classifications_node import ClassificationsNode
from analysis.models.nodes.filters.clinvar_node import ClinVarNode
from analysis.models.nodes.filters.damage_node import DamageNode
from analysis.models.nodes.node_counts import get_node_counts_and_labels_dict
from analysis.models.nodes.node_version import NodeVersion
from analysis.models.nodes.sources.sample_node import SampleNode
from analysis.tests.utils import annotate_variants
from annotation.enums import ClinVarOncogenicity, Pathogenicity, PathogenicityImpact
from annotation.fake_data import get_fake_annotation_version
from annotation.models.models_clinvar import ClinVar
from annotation.models.models_gene_annotation import GeneAnnotation
from annotation.tasks.calculate_sample_stats import calculate_cohort_stats
from classification.enums import SpecialEKeys, SubmissionSource
from classification.models.classification import Classification
from classification.tests.models.test_utils import ClassificationTestUtils
from genes.enums import AnnotationConsortium
from genes.models.models_gene import Gene
from library.django_utils.django_partition import temporary_db_table
from snpdb.enums import BuiltInFilters, SomaticClinicalSignificance
from snpdb.fake_data import create_fake_cohort, make_cohort_genotype
from snpdb.models.models_genome import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import create_mock_allele, slowly_create_test_variant

LABELS = [label for label, _ in BuiltInFilters.FILTER_CHOICES]

# name: (annotation columns, ClinVar columns or None, classifications as (allele origin, evidence))
ROWS = {
    "high_omim_clinvar_p": ({"impact": PathogenicityImpact.HIGH, "omim": True},
                            {"highest_pathogenicity": Pathogenicity.PATHOGENIC},
                            [("germline", {SpecialEKeys.CLINICAL_SIGNIFICANCE: "P"})]),
    "moderate_cosmic": ({"impact": PathogenicityImpact.MODERATE, "cosmic_id": "COSV1"}, None, []),
    "low_cosmic_omim_somatic_tier_2": ({"impact": PathogenicityImpact.LOW, "cosmic_id": "COSV2", "omim": True},
                                       {"somatic_tier": SomaticClinicalSignificance.TIER_2},
                                       [("somatic", {SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE:
                                                     SomaticClinicalSignificance.TIER_1})]),
    "modifier_likely_oncogenic": ({"impact": PathogenicityImpact.MODIFIER},
                                  {"highest_oncogenicity": ClinVarOncogenicity.LIKELY_ONCOGENIC},
                                  [("germline", {SpecialEKeys.CLINICAL_SIGNIFICANCE: "VUS"})]),
    "low_clinvar_vus": ({"impact": PathogenicityImpact.LOW}, {"highest_pathogenicity": Pathogenicity.UNCERTAIN}, []),
    "plain": ({}, None, []),
}

EXPECTED_COUNTS = {
    BuiltInFilters.IMPACT_HIGH_OR_MODERATE: 2,
    BuiltInFilters.COSMIC: 2,
    BuiltInFilters.OMIM: 2,
    BuiltInFilters.CLINVAR: 3,  # P, somatic Tier II, likely oncogenic - not the VUS
    BuiltInFilters.CLASSIFIED: 3,
    BuiltInFilters.CLASSIFIED_PATHOGENIC: 1,
    BuiltInFilters.CLASSIFIED_TIER_1_2: 1,
}


@override_settings(ANALYSIS_NODE_CACHE_Q=False, CLASSIFICATION_ALLOW_UNKNOWN_KEYS=True,
                   CLINGEN_ALLELE_REGISTRY_LOGIN=None)
class CountLabelNodeTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        ClassificationTestUtils.setUp()
        lab, cls.user = ClassificationTestUtils.lab_and_user()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version = get_fake_annotation_version(cls.grch37)
        cls.cohort = create_fake_cohort(cls.user, cls.grch37, name="count_label_nodes")
        cls.sample = cls.cohort.cohortsample_set.get(sample__name='proband').sample
        cgc = cls.cohort.cohort_genotype_collection

        omim_gene = Gene.objects.create(identifier="ENSG_COUNT_LABEL_OMIM",
                                        annotation_consortium=AnnotationConsortium.ENSEMBL)
        gav = cls.annotation_version.gene_annotation_version
        with temporary_db_table(GeneAnnotation, gav.get_partition_table()):
            GeneAnnotation.objects.create(version=gav, gene=omim_gene, omim_terms="OMIM:123456")

        variants_and_columns = {}
        clinvar_version = cls.annotation_version.clinvar_version
        for i, (name, (columns, clinvar_columns, classifications)) in enumerate(ROWS.items()):
            variant = slowly_create_test_variant("1", 1000 + i * 100, "A", "T", cls.grch37)
            make_cohort_genotype(cgc, variant, "E..")
            columns = dict(columns)
            if columns.pop("omim", False):
                columns["gene"] = omim_gene
            variants_and_columns[name] = (variant, columns)
            if clinvar_columns:
                with temporary_db_table(ClinVar, clinvar_version.get_partition_table()):
                    ClinVar.objects.create(version=clinvar_version, variant=variant, clinvar_variation_id=i,
                                           clinvar_allele_id=i, **clinvar_columns)
            if classifications:
                allele = create_mock_allele(variant, cls.grch37)
                for allele_origin, evidence in classifications:
                    data = {SpecialEKeys.C_HGVS: {'value': 'c.301A>C'},
                            SpecialEKeys.ALLELE_ORIGIN: {'value': allele_origin}}
                    data.update({key: {'value': value} for key, value in evidence.items()})
                    classification = Classification.create(user=cls.user, lab=lab, data=data, save=True,
                                                           source=SubmissionSource.API)
                    classification.publish_latest(user=cls.user)
                    classification.allele = allele
                    classification.save()
        annotate_variants(cls.annotation_version.variant_annotation_version, variants_and_columns)

        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(cls.user)

    def _sample_node(self) -> SampleNode:
        node = SampleNode.objects.create(analysis=self.analysis, sample=self.sample)
        node.count = node.get_queryset().count()
        node.status = NodeStatus.READY
        node.save()
        return node

    def _drop(self, parent, label) -> AnalysisNode:
        node = AnalysisNode.count_label_node(self.analysis, label)
        node.save()
        node.add_parent(parent)
        node._cached_parents = None  # Clear stale cache from the first save()
        node.save()
        return node

    def test_badge_count_is_the_dropped_nodes_count(self):
        parent = self._sample_node()
        badge_counts = get_node_counts_and_labels_dict(parent, LABELS)
        self.assertEqual(EXPECTED_COUNTS, badge_counts)
        for label in LABELS:
            with self.subTest(label=label):
                self.assertEqual(badge_counts[label], self._drop(parent, label).get_queryset().count())

    def test_dropped_node_reuses_the_parents_count(self):
        parent = self._sample_node()
        stored_counts = {BuiltInFilters.TOTAL: 100} | {label: 10 + i for i, label in enumerate(LABELS)}
        NodeVersion.objects.filter(pk=parent.node_version.pk).update(load_data={"counts": stored_counts})
        for label in LABELS:
            with self.subTest(label=label):
                node = self._drop(parent, label)
                self.assertEqual(label, node.matches_count_label())
                self.assertEqual(stored_counts[label], node._get_cached_label_count(BuiltInFilters.TOTAL))
                self.assertEqual(stored_counts[label], node._get_cached_label_count(label))

    def test_editing_a_field_stops_the_match(self):
        parent = self._sample_node()
        edits = {
            BuiltInFilters.IMPACT_HIGH_OR_MODERATE: {"impact_required": True},
            BuiltInFilters.OMIM: {"published": True},
            BuiltInFilters.CLINVAR: {"stars_min": 1},
            BuiltInFilters.CLASSIFIED_PATHOGENIC: {"vus": True},
        }
        for label, fields in edits.items():
            with self.subTest(label=label):
                node = self._drop(parent, label)
                for field, value in fields.items():
                    setattr(node, field, value)
                self.assertIsNone(node.matches_count_label())

    def test_classifications_node_restricted_to_a_lab_stops_the_match(self):
        node = self._drop(self._sample_node(), BuiltInFilters.CLASSIFIED)
        lab, _ = ClassificationTestUtils.lab_and_user()
        node.classificationsnodelab_set.create(lab=lab)
        self.assertIsNone(node.matches_count_label())

    def test_badge_click_drops_the_labels_node(self):
        expected_classes = {
            BuiltInFilters.IMPACT_HIGH_OR_MODERATE: DamageNode,
            BuiltInFilters.COSMIC: DamageNode,
            BuiltInFilters.OMIM: DamageNode,
            BuiltInFilters.CLINVAR: ClinVarNode,
            BuiltInFilters.CLASSIFIED: ClassificationsNode,
            BuiltInFilters.CLASSIFIED_PATHOGENIC: ClassificationsNode,
            BuiltInFilters.CLASSIFIED_TIER_1_2: ClassificationsNode,
        }
        self.client.force_login(self.user)
        for label, node_class in expected_classes.items():
            with self.subTest(label=label):
                parent = self._sample_node()
                url = reverse("create_extra_filter_child", kwargs={"analysis_id": self.analysis.pk,
                                                                   "node_id": parent.pk,
                                                                   "extra_filters": label})
                response = self.client.post(url)
                self.assertEqual(200, response.status_code)
                child = AnalysisNode.objects.get_subclass(analysisnode_parent__parent=parent)
                self.assertIsInstance(child, node_class)
                self.assertEqual(label, child.matches_count_label())

    def test_sample_stats_counts_agree_with_the_labels_node(self):
        """ A sample node reads I, O and C from the cohort stats rather than counting them """
        calculate_cohort_stats(self.cohort, self.analysis.annotation_version)
        parent = self._sample_node()
        for label in (BuiltInFilters.IMPACT_HIGH_OR_MODERATE, BuiltInFilters.OMIM, BuiltInFilters.CLINVAR):
            with self.subTest(label=label):
                stats_count = parent._get_cached_label_count(label)
                self.assertEqual(EXPECTED_COUNTS[label], stats_count)
                self.assertEqual(stats_count, self._drop(parent, label).get_queryset().count())
