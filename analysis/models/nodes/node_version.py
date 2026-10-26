""" NodeVersion records a node's load products (counts, timings, stored pks, live data sources) for one
    version; editing a node bumps its version, and everything scoped to a version cascades from this row:
    NodeTask (the scheduler's lease record), NodeCache (stored variants, only IntersectionNode) and the
    column summary cache behind the grid's column counts.

    AnalysisNode is named by string so that analysis_node.py can import these - it reads and writes all of them. """
from typing import TYPE_CHECKING

from django.contrib.postgres.fields import ArrayField
from django.db import models
from django.db.models.aggregates import Count
from django.db.models.deletion import CASCADE
from django.db.models.signals import post_delete
from django.dispatch import receiver
from model_utils.models import TimeStampedModel

from analysis.models.nodes.node_counts import get_node_extra_filters_q
from snpdb.models import VariantCollection

if TYPE_CHECKING:
    from analysis.models.nodes.analysis_node import AnalysisNode


class NodeVersion(TimeStampedModel):
    """ This will be deleted once a node updates, so make all version specific caches cascade delete from this """
    node = models.ForeignKey("AnalysisNode", on_delete=CASCADE)
    version = models.IntegerField(null=False)
    # {source_key: data_version} of the mutable tables this node read at load. Empty = deterministic
    live_data_sources = models.JSONField(default=dict)
    # The exact PK set stored at load, only for nodes <= ANALYSIS_NODE_STORE_ID_SIZE_MAX.
    # When present, load_data["counts"][TOTAL] == len(variant_ids)
    variant_ids = ArrayField(models.IntegerField(), null=True)
    # Products of the node's load:
    #   "counts":  {node count label: count} - the DAG badge counts, eg {"T": 1234, "C": 4, "tag_artefact": 3}.
    #              Labels only ever live under this key, so they can never collide with the keys beside it
    #   "timings": {phase: seconds} of the load - where load_seconds went @see AnalysisNode.load
    load_data = models.JSONField(default=dict)

    class Meta:
        unique_together = ("node", "version")

    @staticmethod
    def get(node: "AnalysisNode"):
        try:
            return NodeVersion.objects.get(node=node, version=node.version)
        except NodeVersion.DoesNotExist:
            node.check_still_valid()
            raise

    @property
    def counts(self) -> dict[str, int]:
        return self.load_data.get("counts", {})

    def __str__(self):
        return f"{self.node.pk} (v{self.version})"


class NodeTask(TimeStampedModel):
    """ Tracks/locks the celery update task for a single node-version, and acts as the lease
        record for the state-driven scheduler (issue #346, mocha-style): the row is claimed by a
        dispatcher (lease_ready_nodes), stamped with a lease_expires, and reclaimed if the owning
        worker dies.

        OneToOne against NodeVersion so it shares that row's lifecycle - the (node, version)
        identity comes from the NodeVersion, and stale tasks are cleaned up automatically when old
        versions are deleted (delete_analysis_old_node_versions cascades). """

    node_version = models.OneToOneField("NodeVersion", on_delete=CASCADE)
    analysis_update_uuid = models.UUIDField()
    celery_task = models.CharField(max_length=36, null=True)
    db_pid = models.IntegerField(null=True)
    # Lease + backoff fields (mocha lease record):
    leased_by = models.CharField(max_length=64, null=True)
    lease_expires = models.DateTimeField(null=True)
    attempt_count = models.IntegerField(default=0)
    last_attempt = models.DateTimeField(null=True)
    run_after = models.DateTimeField(null=True)  # backoff gate: not leasable before this

    def __str__(self):
        return f"NodeTask: {self.analysis_update_uuid} - {self.node_version}"


class NodeCache(models.Model):
    """ This is only used by Intersection node now """
    node_version = models.OneToOneField(NodeVersion, on_delete=CASCADE)
    variant_collection = models.OneToOneField(VariantCollection, on_delete=CASCADE)

    @staticmethod
    def get_or_create_for_node(node: "AnalysisNode") -> tuple['NodeCache', bool]:
        variant_collection = VariantCollection.objects.create(name=f"NodeCache {node.node_version}")
        defaults = {"variant_collection": variant_collection}
        node_cache, created = NodeCache.objects.get_or_create(node_version=node.node_version,
                                                              defaults=defaults)
        if not created:
            variant_collection.delete()
        return node_cache, created

    def __str__(self):
        return f"NodeCache {self.node_version}: {self.variant_collection.get_status_display()}"


@receiver(post_delete, sender=NodeCache)
def post_delete_node_cache(sender, instance, **kwargs):  # pylint: disable=unused-argument
    """ This can sometimes be called multiple times - if node updated again before previous delete is done """
    try:
        if instance.variant_collection:
            instance.variant_collection.delete_related_objects()
            instance.variant_collection.delete()
    except VariantCollection.DoesNotExist:
        # Deleted already
        pass


class NodeColumnSummaryCacheCollection(models.Model):
    node_version = models.ForeignKey(NodeVersion, on_delete=CASCADE)
    variant_column = models.TextField(null=False)
    extra_filters = models.TextField(null=False)

    @staticmethod
    def get_counts_for_node(node, variant_column, extra_filters):
        ncscc, created = NodeColumnSummaryCacheCollection.objects.get_or_create(node_version=node.node_version,
                                                                                variant_column=variant_column,
                                                                                extra_filters=extra_filters)
        if created:
            extra_filters_q = get_node_extra_filters_q(node, extra_filters)
            queryset = node.get_queryset(extra_filters_q)
            count_qs = queryset.values_list(variant_column).distinct().annotate(Count('id'))
            data_list = []
            for value, count in count_qs:
                data = NodeColumnSummaryData(collection=ncscc,
                                             value=value,
                                             count=count)
                data_list.append(data)

            if data_list:
                NodeColumnSummaryData.objects.bulk_create(data_list)
        else:
            data_list = ncscc.nodecolumnsummarydata_set.all()

        counts = {}
        for ncsd in data_list:
            counts[ncsd.value] = ncsd.count

        return counts


class NodeColumnSummaryData(models.Model):
    collection = models.ForeignKey(NodeColumnSummaryCacheCollection, on_delete=CASCADE)
    value = models.TextField(null=True)
    count = models.IntegerField(null=False)
