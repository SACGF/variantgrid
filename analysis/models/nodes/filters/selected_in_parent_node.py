from typing import Optional

from django.db import models
from django.db.models.deletion import CASCADE
from django.db.models.query_utils import Q

from analysis.models.nodes.analysis_node import AnalysisNode
from snpdb.models import Variant


class SelectedInParentNode(AnalysisNode):
    """ This node is enabled, it just can't be created via dropdown """
    disabled = True

    def modifies_parents(self):
        return True  # Always filter - as we never return parent unmodified (except if all variants selected)

    def _get_node_q(self) -> Optional[Q]:
        parent = self.get_single_parent()
        selected_in_parent_qs = NodeVariant.objects.filter(node_id=parent.pk).values_list("variant_id", flat=True)
        return Q(id__in=selected_in_parent_qs)

    def get_node_name(self):
        description = "Selected in parent"
        return description

    @staticmethod
    def get_node_class_label():
        return "Selected In Parent Filter"

    def _get_method_summary(self):
        return "Manually Selected"


class NodeVariant(models.Model):
    variant = models.ForeignKey(Variant, on_delete=CASCADE)
    node = models.ForeignKey(AnalysisNode, on_delete=CASCADE)
