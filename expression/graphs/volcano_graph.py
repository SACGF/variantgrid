from numpy import log10

from expression.models import CuffDiffRecord, CuffDiffFile
from library.utils import sha1_str
import pandas as pd
from snpdb.graphs.graphcache import CacheableGraph


class VolcanoGraph(CacheableGraph):

    def __init__(self, expression_id):
        super().__init__()
        self.expression_id = expression_id

    def get_params_hash(self):
        return sha1_str(self.expression_id)

    def plot(self, ax):
        cuff_diff_file = CuffDiffFile.objects.get(pk=self.expression_id)
        qs = CuffDiffRecord.objects.filter(cuff_diff_file=cuff_diff_file)
        fold_change = 'log2_fold_change'
        p_value = "p_value"
        df = pd.DataFrame.from_records(qs.values(fold_change, p_value))
        x = df[fold_change]
        y = -log10(df[p_value])
        ax.scatter(x, y, color='red', alpha=0.5, s=4)

        significant = -log10(0.05)
        ax.axhline(y=significant, linestyle='--', color='purple')

        title = f"Expression {cuff_diff_file.name} ({cuff_diff_file.annotation_level})"
        ax.set_title(title)
        ax.set_xlabel("log2 fold change")
        ax.set_ylabel("-log10(p-value)")
        ax.set_ylim(0, 16)
        ax.set_xlim(-20, 20)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.yaxis.grid(True)
