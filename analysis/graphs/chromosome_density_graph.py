from analysis.models import AnalysisNode
from library.utils import sha1_str
from snpdb.graphs.chromosome_density_graph import AbstractChromosomeDensityGraph


class ChromosomeDensityGraph(AbstractChromosomeDensityGraph):

    def __init__(self, cmap, node_id):
        super().__init__()
        self.node = AnalysisNode.objects.get_subclass(pk=node_id)
        self.cmap = cmap

    def get_params_hash(self):
        description = f"{self.node.node_version}.{self.cmap}"
        return sha1_str(description)

    def get_queryset(self):
        return self.node.get_queryset()

    def get_genome_build(self):
        return self.node.analysis.genome_build
