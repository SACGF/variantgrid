from typing import Optional, List

from django.db import models
from django.db.models.deletion import SET_NULL
from django.db.models.query_utils import Q
from functools import reduce
import operator

from analysis.models.nodes.analysis_node import AnalysisNode
from annotation.models.models import HumanProteinAtlasTissueSample, \
    HumanProteinAtlasAnnotation
from annotation.models.models_enums import HumanProteinAtlasAbundance, DetectedHumanProteinAtlasAbundance
from genes.models import Gene


class TissueNode(AnalysisNode):
    HUMAN_PROTEIN_ATLAS = 0
    UNIPROTKB = 1

    tissue_sample = models.ForeignKey(HumanProteinAtlasTissueSample, null=True, blank=True, on_delete=SET_NULL)
    min_abundance = models.CharField(max_length=1, choices=DetectedHumanProteinAtlasAbundance.choices, default=DetectedHumanProteinAtlasAbundance.LOW)
    text_tissue = models.TextField(null=True, blank=True)
    accordion_panel = models.IntegerField(default=0)
    disabled = True  # Needs to be made per-genome build see Issue #9

    def _get_configuration_errors(self) -> List:
        errors = super()._get_configuration_errors()
        if self.analysis.annotation_version.human_protein_atlas_version is None:
            msg = "To use the TissueNode you must use an AnnotationVersion with human_protein_atlas_version not equal to None. <a href='javascript:analysisSettings()'>Open Analysis Settings</a>"
            errors.append(msg)
        return errors

    def modifies_parents(self):
        if self.accordion_panel == self.HUMAN_PROTEIN_ATLAS:
            return self.tissue_sample is not None
        return self.text_tissue is not None

    def get_gene_qs(self):
        abundances = HumanProteinAtlasAbundance.get_abundance_or_greater_levels(self.min_abundance)
        hpa_qs = HumanProteinAtlasAnnotation.objects.filter(tissue_sample=self.tissue_sample, abundance__in=abundances)
        Gene.objects.filter(humanproteinatlasannotation__in=hpa_qs)

    def _get_node_q(self) -> Optional[Q]:
        if self.accordion_panel == self.HUMAN_PROTEIN_ATLAS:
            q = Q(variantannotation__gene__in=self.get_gene_qs())
        else:
            filters = []
            for word in self.text_tissue.split():
                f = Q(variantannotation__gene__ensemblgeneannotation__tissue_specificity_from_uniprotkb__icontains=word)
                filters.append(f)

            q = reduce(operator.and_, filters)
        return q

    def _get_method_summary(self):
        if self.modifies_parents():
            if self.accordion_panel == self.HUMAN_PROTEIN_ATLAS:
                min_abundance = self.get_min_abundance_display()
                method_summary = f"Filtering to '{self.tissue_sample}' >= {min_abundance}"
            else:
                words = ", ".join(self.text_tissue.split())
                method_summary = f"Filtering to tissue_specificity_from_uniprotkb contains '{words}'"
        else:
            method_summary = 'No filters applied.'

        return method_summary

    def get_node_name(self):
        name = ''
        if self.modifies_parents():
            if self.accordion_panel == self.HUMAN_PROTEIN_ATLAS:
                min_abundance = self.get_min_abundance_display()
                name = f">= {min_abundance} in {self.tissue_sample}"
            else:
                words = ", ".join(self.text_tissue.split())
                name = f"contains: {words}"
        return name

    @staticmethod
    def get_node_class_label():
        return "Tissue Expression"
