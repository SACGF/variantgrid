from analysis.models import AnalysisNode
from library.utils import sha1_str
from snpdb.graphs.homozygosity_percent_graph import HomozygosityPercentGraph


class NodeHomozygosityPercentGraph(HomozygosityPercentGraph):

    def __init__(self, cmap, node_id):
        self.node = AnalysisNode.objects.get_subclass(pk=node_id)
        sample_ids = self.node.get_sample_ids()
        if len(sample_ids) != 1:
            msg = f"NodeHomozygosityPercentGraph expected node {self.node} to have exactly 1 sample (had {sample_ids})!"
            raise ValueError(msg)

        super().__init__(cmap, sample_ids[0])

    def get_params_hash(self):
        description = f"{self.node.node_version}{self.cmap}"
        return sha1_str(description)

    def get_queryset(self):
        qs = super().get_queryset()  # uses self.sample_id
        return qs.filter(pk__in=self.node.get_queryset())
