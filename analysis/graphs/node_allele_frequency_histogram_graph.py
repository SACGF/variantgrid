from typing import Optional

from django.db.models import Q

from analysis.models.nodes.analysis_node import AnalysisNode
from library.utils import sha1_str
from snpdb.graphs.allele_frequency_graph import AlleleFrequencyHistogramGraph


class NodeAlleleFrequencyHistogramGraph(AlleleFrequencyHistogramGraph):
    MIN_READ_DEPTH = 10

    def __init__(self, _cmap, node_id):
        self.node = AnalysisNode.objects.get_subclass(pk=node_id)
        sample_ids = self.node.get_sample_ids()
        if len(sample_ids) != 1:
            msg = f"NodeAlleleFrequencyHistogramGraph expected node {self.node} to have exactly 1 sample (had {sample_ids})!"
            raise ValueError(msg)

        super().__init__(sample_ids[0], self.MIN_READ_DEPTH)

    def get_params_hash(self):
        description = f"{self.node.node_version}"
        return sha1_str(description)

    def _get_q(self) -> Optional[Q]:
        return Q(pk__in=self.node.get_queryset())
