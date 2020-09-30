from typing import Optional

from django.db import models
from django.db.models import Q

from analysis.models.nodes.analysis_node import AnalysisNode, NodeCount
from analysis.models.nodes.node_counts import get_extra_filters_q
from annotation.models import ClinVarReviewStatus
from snpdb.models.models_enums import BuiltInFilters


class BuiltInFilterNode(AnalysisNode):
    built_in_filter = models.CharField(max_length=1, choices=BuiltInFilters.FILTER_CHOICES, null=True)
    min_clinvar_stars = models.IntegerField(default=0, null=True, blank=True)

    def modifies_parents(self):
        return self.built_in_filter is not None

    def get_extra_filters(self):
        return dict(BuiltInFilters.FILTER_CHOICES)[self.built_in_filter]

    def get_clinvar_stars_q(self):
        review_statuses = []
        for rs, stars in ClinVarReviewStatus.STARS.items():
            if stars >= self.min_clinvar_stars:
                review_statuses.append(rs)
        return Q(clinvar__clinvar_review_status__in=review_statuses)

    def _get_node_q(self) -> Optional[Q]:
        q = get_extra_filters_q(self.analysis.user, self.analysis.genome_build, self.built_in_filter)
        if self.built_in_filter == BuiltInFilters.CLINVAR and self.min_clinvar_stars:
            q &= self.get_clinvar_stars_q()
        return q

    def _get_method_summary(self):
        if self.modifies_parents():
            extra_filters = self.get_extra_filters()
            method_summary = f"Filtering to '{extra_filters}'"
        else:
            method_summary = 'No filters applied as no built in filter selected.'

        return method_summary

    def get_node_name(self):
        name = ''
        if self.modifies_parents():
            extra_filters = self.get_extra_filters()
            name = extra_filters.replace("_", " ")
        return name

    def get_css_classes(self):
        css_classes = super().get_css_classes()
        if self.modifies_parents():
            css_classes.append(f"node-count-{self.built_in_filter}")
        return css_classes

    def get_cached_label_count(self, label):
        count = super().get_cached_label_count(label)
        if count is None:
            if label in [BuiltInFilters.TOTAL, self.built_in_filter]:
                try:
                    parent = self.get_single_parent()
                    parent_node_count = NodeCount.load_for_node(parent, self.built_in_filter)
                    count = parent_node_count.count
                except:
                    pass
        return count

    @staticmethod
    def get_node_class_label():
        return "Built In Filter Node"
