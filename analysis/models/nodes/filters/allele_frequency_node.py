from typing import Optional

from auditlog.registry import auditlog
from django.db import models
from django.db.models import Q
from django.db.models.deletion import SET_NULL

from analysis.models.nodes.analysis_node import AnalysisNode, NodeAlleleFrequencyFilter
from analysis.models.nodes.cohort_mixin import (
    AncestorSampleMixin,
    get_sample_any_zygosity_arg_q_dict,
)
from analysis.models.nodes.node_display import NodeIcon
from patients.models import Patient
from snpdb.models import Sample


class AlleleFrequencyNode(AncestorSampleMixin, AnalysisNode):
    """ AF is stored on NodeAlleleFrequencyFilter so can reuse code for all nodes that have AF filters """
    sample = models.ForeignKey(Sample, null=True, on_delete=SET_NULL)
    patient = models.ForeignKey(Patient, null=True, blank=True, on_delete=SET_NULL)

    def _get_configuration_errors(self) -> list:
        errors = super()._get_configuration_errors()
        if not (self.sample or self.patient):
            errors.append("No sample or patient selected.")
        return errors

    def modifies_parents(self):
        if not (self.sample or self.patient):
            return False
        return self.has_restricted_range()

    def has_restricted_range(self) -> bool:
        """ A full 0-1 range is no filter - @see NodeAlleleFrequencyFilter.get_q """
        try:
            naff = self.nodeallelefrequencyfilter
        except NodeAlleleFrequencyFilter.DoesNotExist:
            return False
        return any(af_range.min > 0 or af_range.max < 1
                   for af_range in naff.nodeallelefrequencyrange_set.all())

    def _get_sample_arg_q_dict(self, sample: Sample) -> dict[Optional[str], dict[str, Q]]:
        """ A VCF with no allele frequency column has nothing to filter on, so its rows come
            through as they are rather than the node emptying itself """
        if not sample.has_allele_frequency:
            return get_sample_any_zygosity_arg_q_dict(sample)
        return NodeAlleleFrequencyFilter.get_sample_arg_q_dict(self, sample)

    def _get_node_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        return self._get_filter_samples_arg_q_dict(self._get_sample_arg_q_dict)

    def _get_method_summary(self):
        if not self.modifies_parents():
            return 'No filters applied'

        af_name = self.nodeallelefrequencyfilter.get_description()
        method_summary = f"Filtering to '{af_name}'"
        if no_af := [s.name for s in self.get_filter_samples() if not s.has_allele_frequency]:
            method_summary += f". No allele frequency to filter on for {', '.join(no_af)}"
        return method_summary

    def get_node_name(self):
        name = ''
        if self.modifies_parents():
            name = self.nodeallelefrequencyfilter.get_description()
            if self.patient:
                name += f"\n{self.patient} ({len(self.get_filter_samples())} samples)"
        return name

    @staticmethod
    def get_help_text() -> str:
        return "Variant Allele Frequency filter"

    def save(self, *args, **kwargs):
        inital_save = not self.pk

        super().save(*args, **kwargs)
        if inital_save:
            # Create a NodeAlleleFrequencyFilter so the filter shows
            naff, created = NodeAlleleFrequencyFilter.objects.get_or_create(node=self)
            if created:
                naff.nodeallelefrequencyrange_set.create(min=0, max=100)

    @staticmethod
    def get_node_class_label():
        return "Allele Frequency"

    @classmethod
    def get_node_class_icon(cls) -> NodeIcon:
        return NodeIcon(fa="fa-solid fa-chart-simple")

    @classmethod
    def get_node_class_label_short(cls) -> str:
        return "AF%"


auditlog.register(AlleleFrequencyNode)
