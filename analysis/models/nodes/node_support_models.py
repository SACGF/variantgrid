""" The models hanging off an AnalysisNode or its NodeVersion: NodeVersion (a node's load products), NodeTask
    (the scheduler's lease record), NodeCache and the column summary cache, NodeVCFFilter and the allele frequency
    filters several nodes share, plus NodeWiki, AnalysisNodeAlleleSource and AnalysisClassification.
    NodeAuditLogMixin puts the node's analysis/node ids on django-auditlog entries for node-owned rows.

    These name AnalysisNode by string so that analysis_node.py can import them - AnalysisNode reads and writes
    every one of them. AnalysisEdge and AnalysisDag stay with AnalysisNode, as the edge's cycle check walks nodes. """
import abc
import operator
from collections import defaultdict
from functools import reduce
from typing import TYPE_CHECKING, Optional

from auditlog.registry import auditlog
from django.contrib.postgres.fields import ArrayField
from django.db import models
from django.db.models.aggregates import Count
from django.db.models.deletion import CASCADE, SET_NULL
from django.db.models.query_utils import Q
from django.db.models.signals import post_delete
from django.dispatch import receiver
from django_extensions.db.models import TimeStampedModel

from analysis.models.enums import GroupOperation
from analysis.models.models_analysis import Analysis
from analysis.models.nodes.node_counts import get_node_extra_filters_q
from classification.models import Classification
from library.django_utils import thread_safe_unique_together_get_or_create
from library.utils import format_percent
from snpdb.models import (
    AlleleSource,
    Sample,
    Variant,
    VariantCollection,
    VCFFilter,
    Wiki,
)

if TYPE_CHECKING:
    from analysis.models.nodes.analysis_node import AnalysisNode


class NodeAuditLogMixin:
    @abc.abstractmethod
    def _get_node(self):
        pass

    def get_additional_data(self):
        """ For django-audit-log """
        node = self._get_node()
        return {
            "analysis_id": node.analysis_id,
            "node_id": node.pk,
        }


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


class NodeWiki(Wiki):
    node = models.OneToOneField("AnalysisNode", on_delete=CASCADE)

    def _get_restricted_object(self):
        return self.node.analysis


class AnalysisNodeAlleleSource(AlleleSource):
    """ Used to link a nodes variants to alleleles and liftover to other builds """
    node = models.ForeignKey("AnalysisNode", null=True, on_delete=SET_NULL)

    def get_genome_build(self):
        if self.node:
            genome_build = self.node.analysis.genome_build
        else:
            genome_build = None
        return genome_build

    def get_variant_qs(self):
        if self.node:
            qs = self.node.get_subclass().get_queryset()
        else:
            qs = Variant.objects.none()
        return qs


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


class NodeCache(models.Model):
    """ This is only used by Intersection node now """
    node_version = models.OneToOneField(NodeVersion, on_delete=CASCADE)
    variant_collection = models.OneToOneField(VariantCollection, on_delete=CASCADE)

    @staticmethod
    def get_or_create_for_node(node: "AnalysisNode") -> tuple['NodeCache', bool]:
        variant_collection = VariantCollection.objects.create(name=f"NodeCache {node.node_version}")
        defaults = {"variant_collection": variant_collection}
        node_cache, created = thread_safe_unique_together_get_or_create(NodeCache, node_version=node.node_version,
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


class NodeVCFFilter(NodeAuditLogMixin, models.Model):
    """ If these exist, they mean use that filter """
    node = models.ForeignKey("AnalysisNode", on_delete=CASCADE)
    vcf_filter = models.ForeignKey(VCFFilter, on_delete=CASCADE, null=True)  # null = 'PASS'

    def _get_node(self):
        return self.node

    @staticmethod
    def get_filter_ids(node):
        """ Returns PASS as None """
        all_nvf_qs = NodeVCFFilter.objects.filter(node_id=node.pk)
        return set(all_nvf_qs.values_list("vcf_filter__filter_id", flat=True))

    @staticmethod
    def has_pass(node) -> bool:
        """ PASS is stored as the node level row with no vcf_filter - it's the one FILTER value
            that means the same thing in every VCF """
        return NodeVCFFilter.objects.filter(node_id=node.pk, vcf_filter__isnull=True).exists()

    @staticmethod
    def get_vcf_filter_ids(node, vcf=None) -> list[tuple]:
        """ (vcf_id, filter_id) for the node's non-PASS rows, optionally for one VCF """
        nvf_qs = NodeVCFFilter.objects.filter(node_id=node.pk, vcf_filter__isnull=False)
        if vcf is not None:
            nvf_qs = nvf_qs.filter(vcf_filter__vcf=vcf)
        return list(nvf_qs.values_list("vcf_filter__vcf_id", "vcf_filter__filter_id"))

    @staticmethod
    def get_filter_codes(node, vcf, pass_only: Optional[bool] = None):
        """ What this VCF lets through: its own ticked codes, plus PASS if the node's PASS row is set.

            PASS is the one FILTER value that means the same thing in every VCF, so a code is never
            translated across them - 'LowDepth' in a DRAGEN small variant VCF is not 'LowDepth' in
            its CNV VCF. pass_only stands in for the node's PASS row, which is how SampleNode lets
            one of its samples decide for itself. """
        nvf_qs = NodeVCFFilter.objects.filter(node_id=node.pk, vcf_filter__vcf=vcf)
        filter_codes = set(nvf_qs.values_list("vcf_filter__filter_code", flat=True))
        if pass_only is None:
            pass_only = NodeVCFFilter.has_pass(node)
        if pass_only:
            filter_codes.add(None)
        return filter_codes

    @staticmethod
    def get_filter_ids_for_other_vcfs(node, vcfs) -> dict:
        """ {vcf: filter_ids} for rows on VCFs the node no longer reads - ticked before its sample or
            cohort changed, and not applied to anything """
        filter_ids_by_vcf = defaultdict(set)
        nvf_qs = NodeVCFFilter.objects.filter(node_id=node.pk, vcf_filter__isnull=False)
        for nvf in nvf_qs.exclude(vcf_filter__vcf__in=vcfs).select_related("vcf_filter__vcf"):
            filter_ids_by_vcf[nvf.vcf_filter.vcf].add(nvf.vcf_filter.filter_id)
        return filter_ids_by_vcf

    @staticmethod
    def carry_over_to_vcf(node, vcf):
        """ A node moved onto another VCF (sample swap, template run) keeps the selection made on its
            old one if the new VCF declares every filter by the same name - a newer pipeline version
            that added filters, or common ones like LowQual. Otherwise the rows stay on the old VCF,
            unapplied, and the node warns (@see CohortMixin.get_warnings) """
        other_vcf_filter_ids = NodeVCFFilter.get_filter_ids_for_other_vcfs(node, [vcf])
        if not other_vcf_filter_ids:
            return

        other_vcf_qs = NodeVCFFilter.objects.filter(node_id=node.pk, vcf_filter__isnull=False) \
            .exclude(vcf_filter__vcf=vcf)
        if NodeVCFFilter.objects.filter(node_id=node.pk, vcf_filter__vcf=vcf).exists():
            other_vcf_qs.delete()  # Chosen on this VCF, which replaces anything left from the old one
            return

        if len(other_vcf_filter_ids) != 1:
            return
        filter_ids = next(iter(other_vcf_filter_ids.values()))
        vcf_filters = list(VCFFilter.objects.filter(vcf=vcf, filter_id__in=filter_ids))
        if len(vcf_filters) != len(filter_ids):
            return

        other_vcf_qs.delete()
        for vcf_filter in vcf_filters:
            NodeVCFFilter.objects.create(node_id=node.pk, vcf_filter=vcf_filter)


class NodeAlleleFrequencyFilter(NodeAuditLogMixin, models.Model):
    """ Used for various nodes """
    node = models.OneToOneField("AnalysisNode", on_delete=CASCADE)
    group_operation = models.CharField(max_length=1, choices=GroupOperation.choices, default=GroupOperation.ANY)

    def _get_node(self):
        return self.node

    def get_q(self, allele_frequency_path: str, allele_frequency_percent: bool) -> Optional[Q]:
        af_q = None
        try:
            filters = []
            for af_range in self.nodeallelefrequencyrange_set.all():
                # Only apply filter if restricted range.
                # Missing value (historical data) == -1 so those will come through
                and_filters = []
                if af_range.min > 0:
                    min_value = af_range.min
                    if allele_frequency_percent:
                        min_value *= 100.0
                    and_filters.append(Q(**{allele_frequency_path + "__gte": min_value}))
                if af_range.max < 1:
                    max_value = af_range.max
                    if allele_frequency_percent:
                        max_value *= 100.0
                    and_filters.append(Q(**{allele_frequency_path + "__lte": max_value}))

                if and_filters:
                    and_q = reduce(operator.and_, and_filters)
                    filters.append(and_q)
            if filters:
                group_op = GroupOperation.get_operation(self.group_operation)
                af_q = reduce(group_op, filters)
        except NodeAlleleFrequencyFilter.DoesNotExist:
            pass

        return af_q

    @staticmethod
    def get_sample_arg_q_dict(node: "AnalysisNode", sample: Sample) -> dict[Optional[str], dict[str, Q]]:
        arg_q_dict = {}
        if sample:
            try:
                alias, allele_frequency_path = sample.get_cohort_genotype_alias_and_field("allele_frequency")
                allele_frequency_percent = sample.vcf.allele_frequency_percent
                if af_q := node.nodeallelefrequencyfilter.get_q(allele_frequency_path, allele_frequency_percent):
                    arg_q_dict[alias] = {str(af_q): af_q}
            except NodeAlleleFrequencyFilter.DoesNotExist:
                pass
        return arg_q_dict

    def get_description(self):
        # TODO: do this properly with group operators etc
        af_ranges = list(self.nodeallelefrequencyrange_set.all())
        if len(af_ranges) == 1:
            description = str(af_ranges[0])
        else:
            description = f"{self.get_group_operation_display()} of {len(af_ranges)} filters"
        return description

    def __str__(self):
        return self.get_description()


class NodeAlleleFrequencyRange(NodeAuditLogMixin, models.Model):
    MIN_VALUE = 0
    MAX_VALUE = 1

    filter = models.ForeignKey(NodeAlleleFrequencyFilter, on_delete=CASCADE)
    min = models.FloatField(null=False)
    max = models.FloatField(null=False)

    def _get_node(self):
        return self.filter.node

    def __str__(self):
        has_min = self.min is not None and self.min > self.MIN_VALUE
        has_max = self.max is not None and self.max < self.MAX_VALUE

        min_perc = format_percent(self.min, is_unit=True)
        max_perc = format_percent(self.max, is_unit=True)

        if has_min and has_max:
            return f"{min_perc} - {max_perc}"
        if has_min:
            return f">={min_perc}"
        if has_max:
            return f"<={max_perc}"
        return ""


class AnalysisClassification(models.Model):
    analysis = models.ForeignKey(Analysis, on_delete=CASCADE)
    classification = models.ForeignKey(Classification, on_delete=CASCADE)



auditlog.register(NodeVCFFilter)
auditlog.register(NodeAlleleFrequencyFilter)
auditlog.register(NodeAlleleFrequencyRange)
