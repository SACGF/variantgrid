from typing import Optional

from django.db import models
from django.db.models import Q

from analysis.models.nodes.analysis_node import AnalysisNode, NodeCount
from analysis.models.nodes.node_counts import get_extra_filters_q
from annotation.models import ClinVarReviewStatus
from snpdb.models.models_enums import BuiltInFilters


class BuiltInFilterNode(AnalysisNode):
    built_in_filter = models.CharField(max_length=1, choices=BuiltInFilters.FILTER_CHOICES, null=True)
    clinvar_stars_min = models.IntegerField(default=0)
    cosmic_count_min = models.IntegerField(default=0)

    def modifies_parents(self):
        return self.built_in_filter is not None

    def get_extra_filters(self):
        return dict(BuiltInFilters.FILTER_CHOICES)[self.built_in_filter]

    def get_clinvar_stars_q(self):
        review_statuses = ClinVarReviewStatus.statuses_gte_stars(self.clinvar_stars_min)
        return Q(clinvar__clinvar_review_status__in=review_statuses)

    def _get_node_q(self) -> Optional[Q]:
        q = get_extra_filters_q(self.analysis.user, self.analysis.genome_build, self.built_in_filter)
        if self.built_in_filter == BuiltInFilters.CLINVAR and self.clinvar_stars_min:
            q &= self.get_clinvar_stars_q()
        elif self.built_in_filter == BuiltInFilters.COSMIC and self.cosmic_count_min:
            q &= Q(variantannotation__cosmic_count__gte=self.cosmic_count_min)
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
            if self.built_in_filter == BuiltInFilters.CLINVAR and self.clinvar_stars_min:
                name += f"\n★ >= {self.clinvar_stars_min}"
            elif self.built_in_filter == BuiltInFilters.COSMIC and self.cosmic_count_min:
                name += f"\n>= {self.cosmic_count_min}"
        return name

    def get_css_classes(self):
        css_classes = super().get_css_classes()
        if self.modifies_parents():
            css_classes.append(f"node-count-{self.built_in_filter}")
        return css_classes

    def _get_cached_label_count(self, label):
        count = super()._get_cached_label_count(label)
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
        return "Built In Filter"
