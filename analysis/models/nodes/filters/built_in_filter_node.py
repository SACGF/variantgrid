from typing import Optional

from auditlog.registry import auditlog
from django.db import models
from django.db.models import Q

from analysis.models.nodes.analysis_node import AnalysisNode
from analysis.models.nodes.node_display import NodeIcon
from annotation.models import ClinVarReviewStatus
from snpdb.enums import BuiltInFilters


class BuiltInFilterNode(AnalysisNode):
    """ Kept for existing analyses and templates - a count badge now drops the node its label is
        configured as (@see AnalysisNode.count_label_node), so no new ones are made """
    disabled = True
    built_in_filter = models.CharField(max_length=1, choices=BuiltInFilters.FILTER_CHOICES, null=True)
    clinvar_stars_min = models.IntegerField(default=0)
    cosmic_count_min = models.IntegerField(default=0)

    def modifies_parents(self):
        return self.built_in_filter is not None

    def get_extra_filters(self):
        return dict(BuiltInFilters.FILTER_CHOICES)[self.built_in_filter]

    def get_clinvar_stars_q(self):
        review_statuses = ClinVarReviewStatus.statuses_gte_stars(self.clinvar_stars_min)
        return Q(clinvar__review_status__in=review_statuses)

    def _get_extra_q(self) -> Optional[Q]:
        if self.built_in_filter == BuiltInFilters.CLINVAR and self.clinvar_stars_min:
            return self.get_clinvar_stars_q()
        if self.built_in_filter == BuiltInFilters.COSMIC and self.cosmic_count_min:
            return Q(variantannotation__cosmic_count__gte=self.cosmic_count_min)
        return None

    def _get_node_q(self) -> Optional[Q]:
        if self.built_in_filter is None:
            return None
        q = self.count_label_q(self.analysis, self.built_in_filter)
        if (extra_q := self._get_extra_q()) is not None:
            q &= extra_q
        return q

    def matches_count_label(self) -> Optional[str]:
        if self._get_extra_q() is not None:
            return None
        return self.built_in_filter

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
                name += "★" * self.clinvar_stars_min
            elif self.built_in_filter == BuiltInFilters.COSMIC and self.cosmic_count_min:
                name += f"\n>= {self.cosmic_count_min}"
        return name

    @staticmethod
    def get_help_text() -> str:
        return "Kept for existing analyses - the Effect, ClinVar and Classifications nodes replace it"

    def get_css_classes(self):
        css_classes = super().get_css_classes()
        if self.modifies_parents():
            css_classes.append(f"node-count-{self.built_in_filter}")
        return css_classes

    @staticmethod
    def get_node_class_label():
        return "Built In Filter"

    @classmethod
    def get_node_class_icon(cls) -> NodeIcon:
        return NodeIcon(fa="fa-solid fa-sliders")

    @classmethod
    def get_node_class_label_short(cls) -> str:
        return "Built-in"


auditlog.register(BuiltInFilterNode)
