import hashlib
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

from library.graphs.graph_utils import ForceMandKIntFormatter
from seqauto.models import QCExecSummary
from seqauto.models.models_enums import QCCompareType
from seqauto.qc.sequencing_run_utils import get_qc_exec_summary_data
from snpdb.graphs.graphcache import CacheableGraph


MANDATORY_COLUMNS = ("mean_coverage_across_genes", "mean_coverage_across_kit", "uniformity_of_coverage",
                     "percent_read_enrichment", "duplicated_alignable_reads", "median_insert", "ts_to_tv_ratio",
                     "number_snps", "snp_dbsnp_percent", "number_indels", "indels_dbsnp_percent")

class QCExecSummaryGraph(CacheableGraph):
    VERSION = 2  # Change to force cache update
    POINT_SIZE = 20
    POINT_COLOR = 'red'
    DEFAULT_BOXPLOT_COLOR = 'blue'
    BOXPLOT_COLORS = {QCCompareType.GOLD_ENRICHMENT_KIT_RUNS: 'gold'}
    EDGE_COLOR = 'black'

    def __init__(self, qc_exec_summary_id, qc_compare_type):
        super().__init__()

        self.qc_exec_summary = QCExecSummary.objects.get(pk=qc_exec_summary_id)
        self.sequencing_run = self.qc_exec_summary.sequencing_run
        self.qc_compare_type = qc_compare_type
        self.bp_color = QCExecSummaryGraph.BOXPLOT_COLORS.get(self.qc_compare_type, QCExecSummaryGraph.DEFAULT_BOXPLOT_COLOR)
        self.type_name = QCCompareType(qc_compare_type).label

    def get_params_hash(self):
        """ This uses get_values_list rather than just hashing params as the underlying DB may have changed """
        sha1 = hashlib.sha1()
        sha1.update(str(QCExecSummaryGraph.VERSION).encode())
        sha1.update(str(self.qc_exec_summary.pk).encode())
        sha1.update(str(self.qc_compare_type).encode())
        for values in self.get_values_list():
            for value in values:
                sha1.update(str(value).encode())

        return sha1.hexdigest()

    def plot(self, ax):
        pass

    def get_columns(self):
        return self.qc_exec_summary.get_coverage_columns() + MANDATORY_COLUMNS

    def get_exec_summary_data(self):
        qc = self.qc_exec_summary
        qc_data = {}
        for col in self.get_columns():
            qc_data[col] = getattr(qc, col)
        return qc_data

    def get_values_list(self):
        run_data = get_qc_exec_summary_data(self.sequencing_run, self.qc_compare_type, self.qc_exec_summary)

        data = []
        for col in self.get_columns():
            data.append((col, run_data[col]))
        return data

    def get_comparison_description(self, num_runs):
        comparison_description = self.type_name
        if self.qc_compare_type in QCCompareType.ENRICHMENT_KIT_TYPES:
            comparison_description += ": " + (self.sequencing_run.get_enrichment_kit_name())
        elif self.qc_compare_type == QCCompareType.SEQUENCER:
            comparison_description += f": {self.sequencing_run.sequencer} "

        comparison_description += f"({num_runs} other QC samples)"
        return comparison_description

    # All of the metrics are on different scales - so need to plot as different graphs
    def plot_figure(self, figure):
        values = self.get_values_list()
        exec_summary_qc = self.get_exec_summary_data()

        try:
            num_runs = len(values[0][1])
        except:
            num_runs = 0

        comparison_description = self.get_comparison_description(num_runs)
        color = QCExecSummaryGraph.POINT_COLOR
        patches = [Line2D([0], [0], color='none', marker='o', markerfacecolor=color),
                   Rectangle((0, 0), 1, 1, fc=self.bp_color)]

        sample_name = self.qc_exec_summary.qc.bam_file.unaligned_reads.sequencing_sample.sample_name
        run_label = f"{sample_name} from {self.sequencing_run}"
        labels = [run_label, comparison_description]
        figure.legend(patches, labels, loc='upper center', bbox_to_anchor=(0.5, 1), numpoints=1)

        num_rows = len(self.get_columns())
        num_cols = 1

        for i, (column_name, column_values) in enumerate(values):
            ax = figure.add_subplot(num_rows, num_cols, i+1)
            run_value = exec_summary_qc[column_name]
            self.plot_enrichment_kit(ax, column_name, column_values, run_value)

        figure.subplots_adjust(left=0.2, bottom=0.05, right=0.95, top=0.85, wspace=0, hspace=1.5)

    def plot_enrichment_kit(self, ax, column_name, column_values, run_value):
        NICER_COLUMN_NAMES = {
            "mean_coverage_across_genes": "mean coverage\nacross genes",
            "mean_coverage_across_kit": "mean coverage\nacross kit",
            "percent_read_enrichment": "percentage\nread enrichment",
            "percent_duplication": "percent\nduplication",
        }

        bplot = ax.boxplot([column_values], vert=False, patch_artist=True, flierprops={"markersize": 1})
        for patch in bplot['boxes']:
            patch.set_facecolor(self.bp_color)
            patch.set_edgecolor(QCExecSummaryGraph.EDGE_COLOR)
            patch.set_zorder(1)

        for whisker in bplot['whiskers']:
            whisker.set(color=QCExecSummaryGraph.EDGE_COLOR, linestyle='-')
            whisker.set_zorder(1)

        ax.scatter([run_value], [1], c=QCExecSummaryGraph.POINT_COLOR, s=[QCExecSummaryGraph.POINT_SIZE], zorder=2)
        ax.xaxis.set_major_formatter(ForceMandKIntFormatter())
        ax.xaxis.set_ticks_position('bottom')
        ax.tick_params(labelsize=6)

        if column_name in NICER_COLUMN_NAMES:
            column_name = NICER_COLUMN_NAMES[column_name]
        else:
            column_name = column_name.replace("_", " ")
        ax.set_yticklabels([column_name], weight='bold')
        ax.yaxis.set_ticks_position('left')

        ax.spines["top"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)
