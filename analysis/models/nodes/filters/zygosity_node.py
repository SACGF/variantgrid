from typing import Optional

from auditlog.registry import auditlog
from django.db import models
from django.db.models.deletion import SET_NULL
from django.db.models.query_utils import Q

from analysis.models.enums import ZygosityNodeZygosity
from analysis.models.gene_counts import NodeGenesCountCollection
from analysis.models.nodes.analysis_node import AnalysisNode
from analysis.models.nodes.cohort_mixin import (
    AncestorSampleMixin,
    get_sample_any_zygosity_arg_q_dict,
)
from analysis.models.nodes.node_display import NodeIcon
from annotation.models.models import VariantAnnotation
from patients.models import Patient
from snpdb.models import Sample


class ZygosityNode(AncestorSampleMixin, AnalysisNode):
    sample = models.ForeignKey(Sample, null=True, on_delete=SET_NULL)
    patient = models.ForeignKey(Patient, null=True, blank=True, on_delete=SET_NULL)
    zygosity = models.CharField(max_length=1, choices=ZygosityNodeZygosity.CHOICES, null=True)
    exclude = models.BooleanField(default=False)

    def _get_configuration_errors(self) -> list:
        errors = super()._get_configuration_errors()
        if not (self.sample or self.patient):
            errors.append("No sample or patient selected.")
        return errors

    def modifies_parents(self):
        return bool((self.sample or self.patient) and self.zygosity)

    def get_zygosity_name(self):
        return dict(ZygosityNodeZygosity.CHOICES)[self.zygosity]

    def _get_sample_arg_q_dict(self, sample: Sample) -> dict[Optional[str], dict[str, Q]]:
        """ One sample's zygosity filter, exclude applied inside it - so in patient mode
            "exclude HET" is each sample's rows that aren't HET, unioned """
        if not sample.has_genotype:
            return get_sample_any_zygosity_arg_q_dict(sample)

        alias, field = sample.get_cohort_genotype_alias_and_field("zygosity")
        q = Q(**{f"{field}": self.zygosity})
        if self.exclude:
            q = ~q
        return {alias: {str(q): q}}

    def _get_node_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        if self.zygosity == ZygosityNodeZygosity.MULTIPLE_HIT:
            parent = self.get_single_parent()
            # Need to pass in kwargs in case we have parent (eg VennNode) that doesn't have needed annotation kwargs
            annotation_kwargs = self.get_annotation_kwargs()
            parent_qs = parent.get_queryset(extra_annotation_kwargs=annotation_kwargs)
            gene_counts_qs = NodeGenesCountCollection.get_or_create_gene_counts_qs_for_node(self, parent_qs)
            genes_with_compound_het_qs = gene_counts_qs.filter(count__gte=2).values('gene_id')
            q = Q(**{VariantAnnotation.GENE_COLUMN + "__in": genes_with_compound_het_qs})
            if self.exclude:
                q = ~q
            return {None: {str(q): q}}

        return self._get_filter_samples_arg_q_dict(self._get_sample_arg_q_dict)

    def _get_method_summary(self):
        if not self.modifies_parents():
            return 'No filters applied as no zygosity selected.'

        zygosity_name = self.get_zygosity_name()
        method_summary = f"Filtering to '{zygosity_name}'"
        if self.zygosity == ZygosityNodeZygosity.MULTIPLE_HIT:
            return method_summary

        filtered = []
        passed_through = []
        for sample in self.get_filter_samples():
            if sample.has_genotype:
                filtered.append(sample.name)
            else:
                passed_through.append(sample.name)
        if filtered:
            method_summary += f" for {', '.join(filtered)}"
        if passed_through:
            # A caller that reports no GT has nothing to filter on, so its rows come through as they are
            method_summary += f". No genotype to filter on for {', '.join(passed_through)}"
        return method_summary

    def get_node_name(self):
        name = ''
        if self.modifies_parents():
            name = self.get_zygosity_name()
            if self.patient:
                name += f"\n{self.patient} ({len(self.get_filter_samples())} samples)"
        return name

    @staticmethod
    def get_help_text() -> str:
        return "Filter by sample zygosity"

    @staticmethod
    def get_node_class_label():
        return "Zygosity"

    @classmethod
    def get_node_class_icon(cls) -> NodeIcon:
        return NodeIcon(symbol="node-icon-zygosity")


auditlog.register(ZygosityNode)
