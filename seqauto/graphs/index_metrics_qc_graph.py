import hashlib

import numpy as np
from django.shortcuts import get_object_or_404

from library.graphs.graph_utils import MandKIntFormatter
from seqauto.models import IlluminaFlowcellQC
from snpdb.graphs.graphcache import CacheableGraph


class IndexMetricsQCGraph(CacheableGraph):
    VERSION = 1
    def __init__(self, illumina_qc_id):
        super().__init__()

        self.illuminaflowcellqc = get_object_or_404(IlluminaFlowcellQC, pk=illumina_qc_id)

    def get_params_hash(self):
        """ This uses get_values_list rather than just hashing params as the underlying DB may have changed """
        sha1 = hashlib.sha1()
        sha1.update(str(IndexMetricsQCGraph.VERSION).encode())
        for values in self.get_values_list():
            for value in values:
                sha1.update(str(value).encode())

        return sha1.hexdigest()

    def figure(self, figure):
        figure.tight_layout()

    def plot(self, ax):
        labels = []
        val = []

        for l, v in self.get_values_list():
            labels.append(l)
            val.append(v)

        pos = np.arange(len(val), 0, -1) - 0.5

        ax.set_title("Index Metrics")

        ax.barh(pos, val, align='center')
        ax.yaxis.set_ticks(pos)
        ax.yaxis.set_ticklabels(labels)
        ax.xaxis.set_major_formatter(MandKIntFormatter())
        ax.set_xlabel("Reads")

    def get_values_list(self):
        qs = self.illuminaflowcellqc.illuminaindexqc_set.all()
        data = []
        for name, reads in qs.values_list("name", "reads"):
            data.append((name, reads))
        return data
