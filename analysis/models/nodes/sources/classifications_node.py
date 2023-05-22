from typing import Optional, List

from django.db import models
from django.db.models import Q, CASCADE

from analysis.models.nodes.analysis_node import AnalysisNode
from classification.enums import ClinicalSignificance
from classification.models.classification import Classification
from snpdb.models import Lab


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
        lab_list = list(self.get_labs())
        return Classification.get_variant_q(self.analysis.user, self.analysis.genome_build, cs_list, lab_list)

    def get_labs(self):
        return Lab.objects.filter(pk__in=self.classificationsnodelab_set.all().values_list("lab", flat=True))

    def get_node_name(self):
        return self.get_node_class_label()

    @staticmethod
    def get_help_text() -> str:
        return "Variants that have been classified. Can filter by pathogenicity."

    @staticmethod
    def get_node_class_label():
        return "Classifications"

    def _get_method_summary(self):
        class_name = ClassificationsNode.get_node_class_label()

        fields = []
        for field in self.FIELD_CLINICAL_SIGNIFICANCE:
            if getattr(self, field):
                fields.append(field)

        method_summary = f"{class_name}, date={self.modified}"
        if len(fields) == len(self.FIELD_CLINICAL_SIGNIFICANCE):
            method_summary += "Any classification"
        else:
            method_summary += ", ".join(fields)

        if labs := self.get_labs():
            method_summary += f". Restricted to labs: {','.join([l.name for l in labs])}"

        return method_summary


class ClassificationsNodeLab(models.Model):
    node = models.ForeignKey(ClassificationsNode, on_delete=CASCADE)
    lab = models.ForeignKey(Lab, on_delete=CASCADE)

