""" Regression tests for SACGF/variantgrid_com#22.

When a Trio (or other source-node input) is deleted, source nodes that referenced it
must:
  - report a configuration error rather than silently producing stale querysets,
  - have their cached q-dicts invalidated (via version bump from a pre_delete signal),
    so downstream nodes do not raise FieldError for now-missing CohortGenotype
    annotation aliases.
"""
from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.models import Analysis, FilterNode, TrioNode
from analysis.models.enums import NodeStatus, TrioInheritance
from annotation.fake_annotation import get_fake_annotation_version
from snpdb.models import GenomeBuild
from snpdb.tests.utils.fake_cohort_data import create_fake_trio


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestSourceNodeTolerance(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser_source_tolerance')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)

        cls.trio = create_fake_trio(cls.user, cls.grch37)
        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(cls.user)

    def test_trio_deletion_bumps_node_version_and_reports_error(self):
        node = TrioNode.objects.create(
            analysis=self.analysis, trio=self.trio,
            inheritance=TrioInheritance.RECESSIVE,
        )
        node.save()
        original_version = node.version

        self.trio.delete()

        node.refresh_from_db()
        self.assertGreater(node.version, original_version,
                           "Trio deletion should bump TrioNode version to invalidate caches")
        self.assertIsNone(node.trio, "Trio FK should be SET_NULL after deletion")

        errors = node.get_errors(flat=True)
        self.assertTrue(errors, "TrioNode with deleted trio should report a configuration error")

        status = node.get_status_from_errors(node.get_errors())
        self.assertEqual(status, NodeStatus.ERROR_CONFIGURATION)

    def test_trio_deletion_invalidates_downstream_node(self):
        """Downstream node version is bumped via the cascade in AnalysisNode.save()."""
        trio_node = TrioNode.objects.create(
            analysis=self.analysis, trio=self.trio,
            inheritance=TrioInheritance.RECESSIVE,
        )
        trio_node.save()

        child = FilterNode.objects.create(analysis=self.analysis)
        child.add_parent(trio_node)
        child.save()
        child_original_version = child.version

        self.trio.delete()

        child.refresh_from_db()
        self.assertGreater(child.version, child_original_version,
                           "Child node version should also be bumped via save() cascade")

    def test_trio_deletion_schedules_analysis_update(self):
        """Affected analyses must have create_and_launch_analysis_tasks queued so the
        DIRTY nodes actually re-evaluate to ERROR_CONFIGURATION (otherwise the grid
        spinner hangs forever waiting for a task that was never launched)."""
        trio_node = TrioNode.objects.create(
            analysis=self.analysis, trio=self.trio,
            inheritance=TrioInheritance.RECESSIVE,
        )
        trio_node.save()

        target = "analysis.tasks.analysis_update_tasks.create_and_launch_analysis_tasks.si"
        with patch(target) as mock_si:
            with self.captureOnCommitCallbacks(execute=True):
                self.trio.delete()

        called_analysis_ids = [call.args[0] for call in mock_si.call_args_list]
        self.assertIn(self.analysis.pk, called_analysis_ids,
                      "Affected analysis should have an update task scheduled on commit")
