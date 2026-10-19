""" NodeVCFFilter - which VCF FILTER values a node lets through (CohortMixin and SampleNode read it).
    PASS is a node level row, every other code belongs to one VCF, and carry_over_to_vcf moves a
    selection when the node's sample or cohort changes. """
from collections import defaultdict
from typing import Optional

from auditlog.registry import auditlog
from django.db import models
from django.db.models.deletion import CASCADE

from analysis.models.nodes.node_audit_log_mixin import NodeAuditLogMixin
from snpdb.models import VCFFilter


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


auditlog.register(NodeVCFFilter)
