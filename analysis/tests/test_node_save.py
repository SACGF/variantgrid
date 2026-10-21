"""
AnalysisNode.save stores `valid` and cascades version bumps to every descendant - what each save has to
query (#2200).
"""
from django.contrib.auth.models import User
from django.db import connection
from django.test import TestCase, override_settings
from django.test.utils import CaptureQueriesContext

from analysis.models import Analysis, AnalysisNode, SampleNode, ZygosityNode
from annotation.fake_data import get_fake_annotation_version
from snpdb.fake_data import create_fake_trio
from snpdb.models import GenomeBuild


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestNodeSave(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        user = User.objects.get_or_create(username="test_node_save")[0]
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(grch37)
        cls.sample = create_fake_trio(user, grch37).proband.sample
        cls.analysis = Analysis(genome_build=grch37)
        cls.analysis.set_defaults_and_save(user)

    def _create_chain(self, length: int) -> list[AnalysisNode]:
        """ SampleNode with a line of ZygosityNodes under it, each loaded fresh as a request would """
        nodes = [SampleNode.objects.create(analysis=self.analysis, sample=self.sample)]
        for _ in range(length):
            node = ZygosityNode.objects.create(analysis=self.analysis)
            node = ZygosityNode.objects.get(pk=node.pk)
            node.add_parent(nodes[-1])
            node.save()
            nodes.append(node)
        return nodes

    def _count_cascade_queries(self, length: int) -> int:
        root = self._create_chain(length)[0]
        root = SampleNode.objects.get(pk=root.pk)
        root.queryset_dirty = True
        with CaptureQueriesContext(connection) as ctx:
            root.save()
        return len(ctx.captured_queries)

    def test_move_does_not_revalidate(self):
        leaf = ZygosityNode.objects.get(pk=self._create_chain(2)[-1].pk)
        self.assertTrue(leaf.valid)
        leaf.x = 100
        # SAVEPOINT/RELEASE, analysis lock, sample (auto node name), auditlog's copy of the old row and its
        # entry, UPDATE node + zygositynode, NodeVersion insert, analysis modified, VCF filter carry over (3).
        # Revalidating would add the parent, analysis and field lookups of get_errors
        with self.assertNumQueries(13):
            leaf.save()
        leaf.refresh_from_db()
        self.assertTrue(leaf.valid)

    def test_internal_error_invalidates_without_revalidating(self):
        leaf = ZygosityNode.objects.get(pk=self._create_chain(1)[-1].pk)
        leaf.errors = "Traceback..."
        leaf.save()
        leaf.refresh_from_db()
        self.assertFalse(leaf.valid)

    def test_cascade_cost_is_linear_in_depth(self):
        """ Each descendant is handed its parent rather than re-walking every ancestor for its samples """
        q2, q4, q6 = (self._count_cascade_queries(n) for n in (2, 4, 6))
        self.assertEqual(q6 - q4, q4 - q2)

    def test_cascade_sets_descendant_samples_and_validity(self):
        nodes = self._create_chain(3)
        root = SampleNode.objects.get(pk=nodes[0].pk)
        root.queryset_dirty = True
        root.save()
        for node in ZygosityNode.objects.filter(pk__in=[n.pk for n in nodes[1:]]):
            self.assertEqual(node.sample, self.sample)
            self.assertTrue(node.valid)

    def test_clone_validates_joined_up_nodes(self):
        self._create_chain(2)
        copy = Analysis.objects.get(pk=self.analysis.pk).clone()
        self.assertEqual([True, True, True], list(copy.analysisnode_set.values_list("valid", flat=True)))
