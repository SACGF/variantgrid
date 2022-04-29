from typing import Optional, List

from django.db import models
from django.db.models.deletion import SET_NULL
from django.db.models.query_utils import Q

from analysis.models import GroupOperation
from analysis.models.nodes.analysis_node import AnalysisNode
from annotation.models.models import HumanProteinAtlasTissueSample, \
    HumanProteinAtlasAnnotation, VariantTranscriptAnnotation
from annotation.models.models_enums import DetectedHumanProteinAtlasAbundance
from genes.models import Gene, GeneSymbol
from genes.models_enums import AnnotationConsortium


class TissueNode(AnalysisNode):
    HUMAN_PROTEIN_ATLAS = 0
    UNIPROTKB = 1

    tissue_sample = models.ForeignKey(HumanProteinAtlasTissueSample, null=True, blank=True, on_delete=SET_NULL)
    min_abundance = models.CharField(max_length=1, choices=DetectedHumanProteinAtlasAbundance.choices,
                                     default=DetectedHumanProteinAtlasAbundance.LOW)
    text_tissue = models.TextField(null=True, blank=True)
    accordion_panel = models.IntegerField(default=0)
    group_operation = models.CharField(max_length=1, choices=GroupOperation.choices, default=GroupOperation.ALL)

    def _get_configuration_errors(self) -> List:
        errors = super()._get_configuration_errors()
        if self.analysis.annotation_version.human_protein_atlas_version is None:
            msg = "To use the TissueNode you must use an AnnotationVersion with human_protein_atlas_version not equal to None. <a href='javascript:analysisSettings()'>Open Analysis Settings</a>"
            errors.append(msg)
        return errors

    def modifies_parents(self):
        if self.uses_hpa():
            return self.tissue_sample is not None
        return self.text_tissue is not None

    def uses_hpa(self):
        return self.accordion_panel == self.HUMAN_PROTEIN_ATLAS

    def get_hpa_qs(self):
        hpa_version = self.analysis.annotation_version.human_protein_atlas_version
        min_value = hpa_version.get_minimum_for_abundance_level(self.min_abundance)
        return HumanProteinAtlasAnnotation.objects.filter(version=hpa_version, tissue_sample=self.tissue_sample,
                                                          value__gte=min_value)

    def _get_hpa_q(self) -> Q:
        va_version = self.analysis.annotation_version.variant_annotation_version
        hpa_qs = self.get_hpa_qs()
        if va_version.annotation_consortium == AnnotationConsortium.ENSEMBL:
            # Use genes directly
            gene_qs = Gene.objects.filter(humanproteinatlasannotation__in=hpa_qs)
        else:
            gene_symbols_qs = GeneSymbol.objects.filter(humanproteinatlasannotation__in=hpa_qs)
            gene_qs = self.analysis.gene_annotation_release.genes_for_symbols(gene_symbols_qs)
        return VariantTranscriptAnnotation.get_overlapping_genes_q(va_version, gene_qs)

    def _get_node_q(self) -> Optional[Q]:
        if self.uses_hpa():
            q = self._get_hpa_q()
        else:
            filters = []
            for word in self.text_tissue.split():
                f = Q(variantannotation__transcript_version__gene_version__hgnc__uniprot__tissue_specificity__icontains=word)
                filters.append(f)
            q = GroupOperation.reduce(filters, self.group_operation)
        return q

    def _get_method_summary(self):
        if self.modifies_parents():
            if self.uses_hpa():
                min_abundance = self.get_min_abundance_display()
                hpa_version = self.analysis.annotation_version.human_protein_atlas_version
                min_value = hpa_version.get_minimum_for_abundance_level(self.min_abundance)
                method_summary = f"Filtering to '{self.tissue_sample}' >= {min_abundance}"
                method_summary += f" ({hpa_version.unit} >= {min_value})"
            else:
                words = ", ".join(self.text_tissue.split())
                method_summary = f"Filtering to tissue_specificity_from_uniprotkb contains '{words}'"
        else:
            method_summary = 'No filters applied.'

        return method_summary

    def get_node_name(self):
        name = ''
        if self.modifies_parents():
            if self.uses_hpa():
                min_abundance = self.get_min_abundance_display()
                name = f">= {min_abundance} in {self.tissue_sample}"
            else:
                words = ", ".join(self.text_tissue.split())
                name = f"contains: {words}"
        return name

    @staticmethod
    def get_help_text() -> str:
        return "Filter based on tissue specific expression (from Human Protein Atlas)"

    @staticmethod
    def get_node_class_label():
        return "Tissue Expression"
