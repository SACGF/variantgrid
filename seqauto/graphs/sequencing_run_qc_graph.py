from django.shortcuts import get_object_or_404
from lazy import lazy
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import hashlib

from library.graphs.graph_utils import ForceMandKIntFormatter
from seqauto.models import SequencingRun
from seqauto.models_enums import QCCompareType
from seqauto.qc.sequencing_run_utils import ILLUMINA_FLOWCELL_QC_COLUMNS, SEQUENCING_RUN_QC_COLUMNS, PAIRED_END_READS, get_q30_col_name, get_sequencing_run_data
from snpdb.graphs.graphcache import CacheableGraph


class SequencingRunQCGraph(CacheableGraph):
    VERSION = 11  # Change to force cache update
    POINT_SIZE = 40
    POINT_COLOR = 'red'
    DEFAULT_BOXPLOT_COLOR = 'blue'
    BOXPLOT_COLORS = {QCCompareType.GOLD_ENRICHMENT_KIT_RUNS: 'gold'}
    EDGE_COLOR = 'black'

    def __init__(self, sequencing_run_id, qc_compare_type):
        super().__init__()

        self.sequencing_run_id = sequencing_run_id
        self.qc_compare_type = qc_compare_type
        self.bp_color = SequencingRunQCGraph.BOXPLOT_COLORS.get(self.qc_compare_type, SequencingRunQCGraph.DEFAULT_BOXPLOT_COLOR)
        qc_types_dict = dict(QCCompareType.CHOICES)
        self.type_name = qc_types_dict[qc_compare_type]

    @lazy
    def sequencing_run(self):
        return get_object_or_404(SequencingRun, pk=self.sequencing_run_id)

    def get_params_hash(self):
        """ This uses get_values_list rather than just hashing params as the underlying DB may have changed """
        sha1 = hashlib.sha1()
        sha1.update(str(SequencingRunQCGraph.VERSION).encode())
        sha1.update(str(self.sequencing_run_id).encode())
        sha1.update(str(self.qc_compare_type).encode())
        for values in self.get_values_list():
            for value in values:
                sha1.update(str(value).encode())

        return sha1.hexdigest()

    def plot(self, ax):
        pass

    def get_sequencing_run_qc(self):
        qc = self.sequencing_run.get_current_sample_sheet().illuminaflowcellqc
        qc_data = {}
        for col in ILLUMINA_FLOWCELL_QC_COLUMNS:
            qc_data[col] = getattr(qc, col)

        for read_q30 in qc.readq30_set.filter(read__in=PAIRED_END_READS):
            col = get_q30_col_name(read_q30.read)
            qc_data[col] = read_q30.percent

        return qc_data

    def get_values_list(self):
        run_data = get_sequencing_run_data(self.sequencing_run, self.qc_compare_type)

        data = []
        for col in SEQUENCING_RUN_QC_COLUMNS:
            data.append((col, run_data[col]))
        return data

    def get_comparison_description(self, num_runs):
        title = self.type_name
        if self.qc_compare_type in QCCompareType.ENRICHMENT_KIT_TYPES:
            title += ": " + (self.sequencing_run.get_enrichment_kit_name())
        elif self.qc_compare_type == QCCompareType.SEQUENCER:
            title += f": {self.sequencing_run.sequencer} "

        title += f"({num_runs} other runs)"
        return title

    # All of the metrics are on different scales - so need to plot as different graphs
    def plot_figure(self, figure):
        values = self.get_values_list()
        sequencing_run_qc = self.get_sequencing_run_qc()

        try:
            num_runs = len(values[0][1])
        except:
            num_runs = 0

        comparison_description = self.get_comparison_description(num_runs)
        color = SequencingRunQCGraph.POINT_COLOR
        patches = [Line2D([0], [0], color='none', marker='o', markerfacecolor=color),
                   Rectangle((0, 0), 1, 1, fc=self.bp_color)]
        labels = [self.sequencing_run, comparison_description]
        figure.legend(patches, labels, loc='upper center', bbox_to_anchor=(0.5, 1), numpoints=1)

        num_rows = len(values)
        num_cols = 1

        for i, (column_name, column_values) in enumerate(values):
            ax = figure.add_subplot(num_rows, num_cols, i + 1)
            run_value = sequencing_run_qc[column_name]
            self.plot_enrichment_kit(ax, column_name, column_values, run_value)

        figure.tight_layout(rect=[0, 0, 1, 0.85])

    def plot_enrichment_kit(self, ax, column_name, column_values, run_value):
        NICER_COLUMN_NAMES = {
            "mean_cluster_density": "mean\ncluster density",
            "mean_pf_cluster_density": "mean P.F.\ncluster density",
            "percentage_of_clusters_pf": "% clusters P.F.",
        }

        bplot = ax.boxplot([column_values], vert=False, widths=0.9, patch_artist=True)
        for patch in bplot['boxes']:
            patch.set_facecolor(self.bp_color)
            patch.set_edgecolor(SequencingRunQCGraph.EDGE_COLOR)
            patch.set_zorder(1)

        for whisker in bplot['whiskers']:
            whisker.set(color=SequencingRunQCGraph.EDGE_COLOR, linewidth=1.2, linestyle='-')
            whisker.set_zorder(1)

        ax.scatter([run_value], [1], c=SequencingRunQCGraph.POINT_COLOR, s=[SequencingRunQCGraph.POINT_SIZE], zorder=2)
        if column_name in NICER_COLUMN_NAMES:
            column_name = NICER_COLUMN_NAMES[column_name]
        else:
            column_name = column_name.replace("_", " ")

        ax.set_yticklabels([column_name], weight='bold')
        ax.xaxis.set_major_formatter(ForceMandKIntFormatter())

        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        #ax.get_yaxis().set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)
