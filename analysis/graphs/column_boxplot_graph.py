import numpy as np
import pandas as pd

from analysis.models import AnalysisNode
from library.utils import sha1_str
from snpdb.graphs.graphcache import CacheableGraph


class ColumnBoxplotGraph(CacheableGraph):

    def __init__(self, node_id, label, variant_column):
        super().__init__()
        self.node = AnalysisNode.objects.get_subclass(pk=node_id)
        self.label = label
        self.variant_column = variant_column

    def get_params_hash(self):
        description = f"{self.node.node_version}{self.variant_column}"
        return sha1_str(description)

    def plot(self, ax):
        qs = self.node.get_queryset()
        values = qs.values(*[self.variant_column])
        df = pd.DataFrame.from_records(values)
        # I get: ImportError: No module named _backend_gdk
        # df.boxplot(ax=ax)

        data = []
        for c in df.columns:
            series = df[c]
            ok = np.isfinite(series * 1.01)  # Handle almost infinity
            series = series[ok]
            data.append(series)

        ax.boxplot(data)
        ax.set_xticklabels([self.label])
