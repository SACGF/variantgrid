from typing import Optional

from django.db import models
from django.db.models import Q
from django.db.models.deletion import SET_NULL

from analysis.models.nodes.analysis_node import AnalysisNode
from snpdb.models import Variant


class AllVariantsNode(AnalysisNode):
    max_variant = models.ForeignKey(Variant, null=True, on_delete=SET_NULL)  # Store this, so we know when we are out of date
    min_inputs = 0
    max_inputs = 0

    def _get_node_q(self) -> Optional[Q]:
        """ Restrict to analysis genome build """
        return Variant.get_contigs_q(self.analysis.genome_build)

    def get_node_name(self):
        return ''

    @staticmethod
    def get_node_class_label():
        return "All Variants"

    def _get_method_summary(self):
        class_name = AllVariantsNode.get_node_class_label()
        max_id = self.max_variant.pk
        method_summary = f"{class_name}, date={self.modified}, max_variant={max_id}"
        return method_summary

    def _get_configuration_errors(self):
        errors = []
        if not self.max_variant:
            errors.append("Not Saved")
        return errors
