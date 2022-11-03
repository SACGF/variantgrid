from typing import Optional, List, Dict, Set

from django.db import models
from django.db.models import Q
from django.db.models.deletion import SET_NULL

from analysis.models.nodes.analysis_node import AnalysisNode, NodeAlleleFrequencyFilter
from analysis.models.nodes.cohort_mixin import AncestorSampleMixin
from snpdb.models import Sample


class AlleleFrequencyNode(AncestorSampleMixin, AnalysisNode):
    sample = models.ForeignKey(Sample, null=True, on_delete=SET_NULL)

    def _get_configuration_errors(self) -> List:
        errors = super()._get_configuration_errors()
        if not self.sample:
            errors.append("No sample selected.")
        return errors

    def modifies_parents(self):
        return NodeAlleleFrequencyFilter.get_sample_arg_q_dict(self, self.sample)

    def _get_node_arg_q_dict(self) -> Dict[Optional[str], Dict[str, Q]]:
        return NodeAlleleFrequencyFilter.get_sample_arg_q_dict(self, self.sample)

    def _get_method_summary(self):
        if self.modifies_parents():
            af_name = self.nodeallelefrequencyfilter.get_description()
            method_summary = f"Filtering to '{af_name}'"
        else:
            method_summary = 'No filters applied'

        return method_summary

    def get_node_name(self):
        name = ''
        if self.modifies_parents():
            name = self.nodeallelefrequencyfilter.get_description()
        return name

    @staticmethod
    def get_help_text() -> str:
        return "Variant Allele Frequency filter"

    def save(self, **kwargs):
        inital_save = not self.pk

        super().save(**kwargs)
        if inital_save:
            # Create a NodeAlleleFrequencyFilter so the filter shows
            naff, created = NodeAlleleFrequencyFilter.objects.get_or_create(node=self)
            if created:
                naff.nodeallelefrequencyrange_set.create(min=0, max=100)

    @staticmethod
    def get_node_class_label():
        return "Allele Frequency"
