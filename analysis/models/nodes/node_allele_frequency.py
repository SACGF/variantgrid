""" NodeAlleleFrequencyFilter and its NodeAlleleFrequencyRange rows - the sample allele frequency ranges
    several nodes share (AlleleFrequencyNode, SampleNode, cohort nodes), combined with the filter's
    GroupOperation into a Q on the cohort genotype's allele frequency. """
import operator
from functools import reduce
from typing import TYPE_CHECKING, Optional

from auditlog.registry import auditlog
from django.db import models
from django.db.models.deletion import CASCADE
from django.db.models.query_utils import Q

from analysis.enums import GroupOperation
from analysis.models.nodes.node_audit_log_mixin import NodeAuditLogMixin
from library.utils import format_percent
from snpdb.models import Sample

if TYPE_CHECKING:
    from analysis.models.nodes.analysis_node import AnalysisNode


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


auditlog.register(NodeAlleleFrequencyFilter)
auditlog.register(NodeAlleleFrequencyRange)
