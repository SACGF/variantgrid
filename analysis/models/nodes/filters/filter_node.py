from typing import Optional

from django.db import models
from django.db.models import Q
from django.db.models.deletion import CASCADE
import json

from analysis.models.nodes.analysis_node import AnalysisNode
from library import jqgrid
from snpdb.models import VariantGridColumn


# TODO: This node has quite a few redundant operations - eg it will filter the queryset
# By the filter, and then have it re-applied by jqgrid. Maybe we could override the filter
# method if already filtered, to save work.
# But: Get it right, then get it fast.....
class FilterNode(AnalysisNode):
    group_operation = models.CharField(max_length=3)  # "OR" or "AND"

    def modifies_parents(self):
        return self.filternodeitem_set.count()

    def _get_node_q(self) -> Optional[Q]:
        class FakeFilterGrid(jqgrid.JqGrid):
            fields = ["id"]  # Needs to be defined

        # This filter uses JQGrid's built in query filter.
        # Load stored params from the DB, convert to JSON and send to a fake request
        fake_filter_grid = FakeFilterGrid(case_sensitive_search=False)
        return fake_filter_grid.get_q(self.get_filters())

    def get_extra_grid_config(self):
        existing_extra_config = super().get_extra_grid_config()

        if self.filternodeitem_set.count():
            existing_extra_config['search'] = True
            post_data = existing_extra_config.get('postData', {})
            post_data.update({'filters': self.get_filters_json()})
            existing_extra_config['postData'] = post_data
        return existing_extra_config

    def get_rules(self):
        rules = []
        for rule in self.filternodeitem_set.order_by('sort_order'):
            rules.append({'field': rule.field, 'op': rule.operation, 'data': rule.data})
        return rules

    def get_filters_json(self):
        """ Grid expects quoted json string (so we'll double json this after it is put in dict) """
        return json.dumps(self.get_filters())

    def get_filters(self):
        return {'groupOp': self.group_operation, 'rules': self.get_rules()}

    def _get_method_summary(self):
        rules_summary = []
        for rule in self.get_rules():
            rule["op"] = jqgrid.format_operation(rule["op"])
            rules_summary.append("%(field)s %(op)s %(data)s" % rule)

        joiner = f" {self.group_operation} "
        return joiner.join(rules_summary)

    def get_node_name(self):
        fn_items = list(self.filternodeitem_set.all())  # Evaluate once
        num_filters = len(fn_items)
        if num_filters == 0:
            node_name = ''
        elif num_filters == 1:
            node_name = str(fn_items[0])
        else:
            node_name = f"{num_filters} filters"

        return node_name

    def _get_inherited_colmodel_overrides(self):
        # Don't allow searching on inherited columns, as this causes an extra join
        extra_overrides = super()._get_inherited_colmodel_overrides()
        extra_columns = self._get_inherited_columns()

        for col in extra_columns:
            data = extra_overrides.get(col) or {}
            data['search'] = False
            extra_overrides[col] = data

        return extra_overrides

    @staticmethod
    def get_node_class_label():
        return "Filter"

    def save_clone(self):
        filter_items = list(self.filternodeitem_set.all())
        copy = super().save_clone()
        for fi in filter_items:
            copy.filternodeitem_set.create(sort_order=fi.sort_order,
                                           operation=fi.operation,
                                           field=fi.field,
                                           data=fi.data)
        return copy


class FilterNodeItem(models.Model):
    filter_node = models.ForeignKey(FilterNode, on_delete=CASCADE)
    sort_order = models.IntegerField()
    operation = models.CharField(max_length=2)
    field = models.TextField()
    data = models.TextField()

    @property
    def column(self):
        return VariantGridColumn.objects.get(variant_column=self.field)

    def __str__(self):
        op = jqgrid.format_operation(self.operation)
        description = f"{self.column.grid_column_name} {op}"
        if self.data:
            description += " " + self.data
        return description
