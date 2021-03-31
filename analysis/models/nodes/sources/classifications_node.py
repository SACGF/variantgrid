from typing import Optional

from django.db import models
from django.db.models import Q

from analysis.models.nodes.analysis_node import AnalysisNode
from classification.enums import ClinicalSignificance
from classification.models.classification import Classification


class ClassificationsNode(AnalysisNode):
    # Default is to show all
    other = models.BooleanField(default=True, blank=True)
    benign = models.BooleanField(default=True, blank=True)
    likely_benign = models.BooleanField(default=True, blank=True)
    vus = models.BooleanField(default=True, blank=True)
    likely_pathogenic = models.BooleanField(default=True, blank=True)
    pathogenic = models.BooleanField(default=True, blank=True)

    FIELD_CLINICAL_SIGNIFICANCE = {
        'other': ClinicalSignificance.OTHER,
        'benign': ClinicalSignificance.BENIGN,
        'likely_benign': ClinicalSignificance.LIKELY_BENIGN,
        'vus': ClinicalSignificance.VUS,
        'likely_pathogenic': ClinicalSignificance.LIKELY_PATHOGENIC,
        'pathogenic': ClinicalSignificance.PATHOGENIC,
    }
    min_inputs = 0
    max_inputs = 0

    def _get_node_q(self) -> Optional[Q]:
        cs_list = []
        for field, cs in self.FIELD_CLINICAL_SIGNIFICANCE.items():
            if getattr(self, field):
                cs_list.append(cs)
        return Classification.get_variant_q(self.analysis.user, self.analysis.genome_build, cs_list)

    def get_node_name(self):
        return self.get_node_class_label()

    @staticmethod
    def get_node_class_label():
        return "Classifications"

    def _get_method_summary(self):
        class_name = ClassificationsNode.get_node_class_label()

        fields = []
        for field, cs in self.FIELD_CLINICAL_SIGNIFICANCE.items():
            if getattr(self, field):
                fields.append(field)

        method_summary = f"{class_name}, date={self.modified}"
        if len(fields) == len(self.FIELD_CLINICAL_SIGNIFICANCE):
            method_summary += "Any classification"
        else:
            method_summary += ", ".join(fields)
        return method_summary
