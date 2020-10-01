from typing import Optional

from django.db import models
from django.db.models import Q

from analysis.models.nodes.analysis_node import AnalysisNode
from classification.enums import ClinicalSignificance
from classification.enums.variant_classification_enums import ClinicalSignificanceComparison
from classification.models.variant_classification import VariantClassification


class ClassificationsNode(AnalysisNode):
    clinical_significance = models.CharField(max_length=1, choices=ClinicalSignificance.CHOICES, null=True, blank=True)
    comparison = models.CharField(max_length=1, choices=ClinicalSignificanceComparison.CHOICES, default=ClinicalSignificanceComparison.EQUAL)
    min_inputs = 0
    max_inputs = 0

    def _get_node_q(self) -> Optional[Q]:
        cs_list = ClinicalSignificanceComparison.get_clinical_significance_list(self.clinical_significance, self.comparison)
        return VariantClassification.get_variant_q(self.analysis.user, self.analysis.genome_build, cs_list)

    def get_node_name(self):
        return self.get_node_class_label()

    @staticmethod
    def get_node_class_label():
        return "Classifications"

    def _get_method_summary(self):
        class_name = ClassificationsNode.get_node_class_label()
        method_summary = f"{class_name}, date={self.modified}"
        if self.clinical_significance:
            method_summary += f"clinical significance={self.clinical_significance}"
        else:
            method_summary += "Any classification"
        return method_summary
