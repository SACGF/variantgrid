import uuid
from datetime import timedelta
from unittest import mock

from django.contrib.auth.models import User
from django.db.models import F
from django.test import TestCase
from django.urls import reverse
from django.utils import timezone

from analysis.enums import NodeStatus
from analysis.models import AllVariantsNode, Analysis, AnalysisNode, NodeTask
from analysis.models.nodes.node_version import NodeVersion
from analysis.scheduler_status import get_analyses_with_work, get_scheduler_status
from analysis.tests.utils import AnalysisSetupMixin

QUEUE_DEPTHS = {"analysis_workers": {"messages": 3, "consumers": 2},
                "scheduling_single_worker": {"error": "OperationalError"}}


@mock.patch("analysis.scheduler_status.broker_queue_depths", return_value=QUEUE_DEPTHS)
class TestAnalysisStatus(AnalysisSetupMixin, TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.now = timezone.now()
        live = cls.now + timedelta(minutes=5)
        expired = cls.now - timedelta(minutes=5)
        future = cls.now + timedelta(minutes=5)
        past = cls.now - timedelta(minutes=5)

        cls._node(cls.analysis, NodeStatus.DIRTY)  # waiting
        cls._node(cls.analysis, NodeStatus.DIRTY, run_after=past)  # waiting - backoff over
        cls._node(cls.analysis, NodeStatus.DIRTY, run_after=future)  # backing off
        cls._node(cls.analysis, NodeStatus.QUEUED, lease_expires=live)  # dispatched
        cls._node(cls.analysis, NodeStatus.LOADING, lease_expires=live)  # running
        cls._node(cls.analysis, NodeStatus.LOADING_CACHE, lease_expires=live)  # running
        cls._node(cls.analysis, NodeStatus.QUEUED, lease_expires=expired)  # lapsed
        cls._node(cls.analysis, NodeStatus.LOADING)  # lapsed - no lease at all
        cls._node(cls.analysis, NodeStatus.READY, lease_expires=live)  # not loading - ignored

        # A live lease on a superseded version says nothing about the current one
        superseded = cls._node(cls.analysis, NodeStatus.LOADING, lease_expires=live)
        AnalysisNode.objects.filter(pk=superseded.pk).update(version=F("version") + 1)

        user = User.objects.get_or_create(username="test_analysis_status_other")[0]
        cls.other_analysis = Analysis(genome_build=cls.grch37, name="Other analysis")
        cls.other_analysis.set_defaults_and_save(user)
        old_node = cls._node(cls.other_analysis, NodeStatus.DIRTY)
        cls.oldest_modified = cls.now - timedelta(days=1)
        AnalysisNode.objects.filter(pk=old_node.pk).update(modified=cls.oldest_modified)

    @classmethod
    def _node(cls, analysis, status, lease_expires=None, run_after=None):
        node = AllVariantsNode.objects.create(analysis=analysis)
        if lease_expires or run_after:
            node_version, _ = NodeVersion.objects.get_or_create(node_id=node.pk, version=node.version)
            NodeTask.objects.create(node_version=node_version, analysis_update_uuid=uuid.uuid4(),
                                    lease_expires=lease_expires, run_after=run_after)
        AnalysisNode.objects.filter(pk=node.pk).update(status=status)
        return node

    def test_scheduler_status_buckets(self, _broker_queue_depths):
        status = get_scheduler_status(self.now)
        self.assertEqual(status.waiting, 3)  # 2 here + the other analysis's node
        self.assertEqual(status.backing_off, 1)
        self.assertEqual(status.dispatched, 1)
        self.assertEqual(status.running, 2)
        self.assertEqual(status.lapsed, 3)
        self.assertEqual(status.analyses_with_work, 2)
        self.assertEqual(status.queues, QUEUE_DEPTHS)

    def test_analyses_with_work(self, _broker_queue_depths):
        works = get_analyses_with_work(self.now)
        self.assertEqual([w.analysis_id for w in works], [self.other_analysis.pk, self.analysis.pk])
        other, this = works
        self.assertEqual(other.oldest_modified, self.oldest_modified)
        self.assertEqual((other.waiting, other.dispatched, other.running, other.lapsed), (1, 0, 0, 0))
        # per analysis, waiting includes the node backing off
        self.assertEqual((this.waiting, this.dispatched, this.running, this.lapsed), (3, 1, 2, 3))

    def test_page_and_old_url_redirect(self, _broker_queue_depths):
        superuser = User.objects.create_superuser("test_analysis_status_admin")
        self.client.force_login(superuser)
        url = reverse("analysis_status")
        response = self.client.get(url)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Other analysis")
        self.assertContains(response, "OperationalError")

        response = self.client.get("/analysis/analysis_issues?error_contains=x")
        self.assertRedirects(response, f"{url}?error_contains=x", status_code=301)
