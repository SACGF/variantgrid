import json
from typing import Optional

from auditlog.registry import auditlog
from django.db import models
from django.db.models import Q
from django.db.models.deletion import CASCADE

from analysis.models.nodes.analysis_node import AnalysisNode, NodeAuditLogMixin
from analysis.models.nodes.node_display import NodeIcon
from library.django_utils.filter_rules import format_operation, rules_to_q
from snpdb.models import VariantGridColumn


# TODO: This node has quite a few redundant operations - e.g. it will filter the queryset
# By the filter, and then have it re-applied by the grid engine. Maybe we could override the filter
# method if already filtered, to save work.
# But: Get it right, then get it fast.....
class FilterNode(AnalysisNode):
    group_operation = models.CharField(max_length=3)  # "OR" or "AND"

    def modifies_parents(self):
        return self.filternodeitem_set.exists()

    def _get_node_q(self) -> Optional[Q]:
        # The same rule -> Q conversion the grid's column filters use, so a FilterNode and a grid
        # column filter on the same rules produce the same queryset
        return rules_to_q(self.get_filters())

    def get_grid_post_data(self) -> dict:
        post_data = super().get_grid_post_data()
        if self.filternodeitem_set.exists():
            # The grid filters on these rules, and the editor's filter builder loads them back
            post_data['filters'] = self.get_filters_json()
        return post_data

    def get_rules(self):
        rules = []
        for rule in self.filternodeitem_set.order_by('sort_order'):
            rules.append({'field': rule.field, 'op': rule.operation, 'data': rule.data})
        return rules

    def get_filters_json(self):
        """ The grid sends the rules up as a JSON string, so they go into postData already encoded """
        return json.dumps(self.get_filters())

    def get_filters(self):
        return {'groupOp': self.group_operation, 'rules': self.get_rules()}

    def _get_method_summary(self):
        rules_summary = []
        for rule in self.get_rules():
            rule["op"] = format_operation(rule["op"])
            rules_summary.append(f"{rule['field']} {rule['op']} {rule['data']}")

        joiner = f" {self.group_operation} "
        return joiner.join(rules_summary)

    def get_node_name(self):
        node_name = ''
        if self.pk:
            fn_items = list(self.filternodeitem_set.all())  # Evaluate once
            num_filters = len(fn_items)
            if num_filters == 0:
                node_name = ''
            elif num_filters == 1:
                node_name = str(fn_items[0])
            else:
                node_name = f"{num_filters} filters"

        return node_name

    @staticmethod
    def get_help_text() -> str:
        return "Filter based on column values"

    def _get_inherited_columns(self):
        # Don't offer a filter on inherited columns, as this causes an extra join
        extra_columns = super()._get_inherited_columns()
        for rich_column in extra_columns:
            rich_column.column_filter = None
        return extra_columns

    @staticmethod
    def get_node_class_label():
        return "Filter"

    @classmethod
    def get_node_class_icon(cls) -> NodeIcon:
        return NodeIcon(fa="fa-solid fa-filter")

    def save_clone(self):
        filter_items = list(self.filternodeitem_set.all())
        copy = super().save_clone()
        for fi in filter_items:
            copy.filternodeitem_set.create(sort_order=fi.sort_order,
                                           operation=fi.operation,
                                           field=fi.field,
                                           data=fi.data)
        return copy


class FilterNodeItem(NodeAuditLogMixin, models.Model):
    filter_node = models.ForeignKey(FilterNode, on_delete=CASCADE)
    sort_order = models.IntegerField()
    operation = models.CharField(max_length=2)
    field = models.TextField()
    data = models.TextField()

    @property
    def column(self):
        return VariantGridColumn.objects.get(variant_column=self.field)

    def _get_node(self):
        return self.filter_node

    def __str__(self):
        op = format_operation(self.operation)
        description = f"{self.column.grid_column_name} {op}"
        if self.data:
            description += " " + self.data
        return description


auditlog.register(FilterNode)
auditlog.register(FilterNodeItem)
