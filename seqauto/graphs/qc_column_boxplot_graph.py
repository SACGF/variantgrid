import hashlib
import logging
from collections import defaultdict

from seqauto.models import QC, QCColumn, EnrichmentKit
from seqauto.models.models_enums import QCGraphEnrichmentKitSeparationChoices
from snpdb.graphs.graphcache import CacheableGraph
from snpdb.models import DataState

SEQUENCING_SAMPLE_PATH = 'bam_file__unaligned_reads__sequencing_sample'
SAMPLE_SHEET_PATH = SEQUENCING_SAMPLE_PATH + '__sample_sheet'
PADDING = 0.5


class QCColumnBoxPlotGraph(CacheableGraph):

    def __init__(self, qc_column_id, enrichment_kit_separation, enrichment_kit_id, percent):
        super().__init__()
        self.qc_column = QCColumn.objects.get(pk=qc_column_id)
        self.enrichment_kit_separation = enrichment_kit_separation
        if enrichment_kit_id:
            self.enrichment_kit = EnrichmentKit.objects.get(pk=enrichment_kit_id)
        else:
            self.enrichment_kit = None
        self.percent = bool(percent)

        logging.info("QC_column: %s", self.qc_column)
        logging.info("EnrichmentKit Sep: %s", self.enrichment_kit_separation)
        logging.info("EnrichmentKit: %s", self.enrichment_kit)
        logging.info("Percent: %s", self.percent)

    def get_params_hash(self):
        """ This uses get_values_list rather than just hashing params as the underlying DB may have changed """
        sha1 = hashlib.sha1()
        for values in self.get_values_list():
            for value in values:
                sha1.update(str(value).encode())

        return sha1.hexdigest()

    def get_values_list(self):
        data_state = self.qc_column.qc_type.qc_object_path + "__data_state"
        base_qs = QC.objects.filter(**{data_state: DataState.COMPLETE})

        if self.enrichment_kit_separation == QCGraphEnrichmentKitSeparationChoices.ALL_ENRICHMENT_KITS:
            return self.get_values_for_qs(base_qs)
        elif self.enrichment_kit_separation == QCGraphEnrichmentKitSeparationChoices.SELECTED_ENRICHMENT_KIT:
            return self.get_values_for_enrichment_kit(base_qs, self.enrichment_kit)
        elif self.enrichment_kit_separation == QCGraphEnrichmentKitSeparationChoices.SEPARATED_ENRICHMENT_KITS:
            enrichment_kit_values = []  # (enrichment_kit name, data)
            for enrichment_kit in EnrichmentKit.objects.all().order_by("name"):
                enrichment_kit_data = self.get_values_for_enrichment_kit(base_qs, enrichment_kit)
                if enrichment_kit_data:
                    enrichment_kit_values.append((enrichment_kit.name, enrichment_kit_data))
            return enrichment_kit_values

    def get_values_for_enrichment_kit(self, qs, enrichment_kit):
        ENRICHMENT_KIT_PATH = SEQUENCING_SAMPLE_PATH + "__enrichment_kit"

        qs = qs.filter(**{ENRICHMENT_KIT_PATH: enrichment_kit})
        return self.get_values_for_qs(qs)

    def get_values_for_qs(self, qs):
        """ Returns array of tuples ('sequencing_run', array of values) """

        SEQUENCING_RUN = SAMPLE_SHEET_PATH + "__sequencing_run"

        def get_field(f):
            return self.qc_column.qc_type.qc_object_path + "__" + f

        path = get_field(self.qc_column.field)
        qs = qs.filter(**{path + "__isnull": False})

        args = [SEQUENCING_RUN, path]
        if self.percent:
            total_field = self.qc_column.qc_type.total_field
            if total_field is None:
                msg = f"Asked for percentage for {self.qc_column} ({self.qc_column.qc_type}) with no total field!"
                raise ValueError(msg)
            total_path = get_field(total_field)
            args.append(total_path)

        sequencing_run_values = defaultdict(list)
        for values in qs.values_list(*args):
            sequencing_run = values[0]
            if self.percent:
                (val, total) = values[1:]
                if total:
                    val = 100.0 * val / total
                elif val:
                    msg = f"Val was '{val}' with total field {total_field} of {total}"
                    raise ValueError(msg)
            else:
                val = values[1]

            # TODO: Right here we want to handle percent (2 values) and normal (1 value)
            # Being put into a scalar

            sequencing_run_values[sequencing_run].append(val)

        sr_data = []
        for sr in reversed(sorted(sequencing_run_values)):  # Want top to bottom:
            data = (sr, sequencing_run_values[sr])
            sr_data.append(data)

        return sr_data

    def plot_figure(self, figure):
        if self.enrichment_kit_separation == QCGraphEnrichmentKitSeparationChoices.SEPARATED_ENRICHMENT_KITS:
            self.plot_multiple(figure)
        else:
            super().plot_figure(figure)

    # This will only be called for single figure plots (ie not SEPARATED_ENRICHMENT_KITS)
    def plot(self, ax):
        enrichment_kit_values = self.get_values_list()
        enrichment_kit_name = self.enrichment_kit.name if self.enrichment_kit else None
        self.plot_enrichment_kit(ax, enrichment_kit_name, enrichment_kit_values)

    def plot_enrichment_kit(self, ax, enrichment_kit_name, enrichment_kit_values):
        if enrichment_kit_name:
            title = f"EnrichmentKit: {enrichment_kit_name}"
        else:
            title = "All Enrichment Kits"
        ax.set_title(title)

        y_values = []
        y_labels = []

        ax.set_xlabel(self.qc_column.name)

        for sequencing_run, values in enrichment_kit_values:
            y_labels.append(sequencing_run)
            y_values.append(values)

        if self.percent:
            title += " %% of %s" % self.qc_column.qc_type.total_field

        ax.boxplot(y_values, vert=False)
        ax.set_yticklabels(y_labels)

    def plot_multiple(self, figure):
        enrichment_kit_values = self.get_values_list()

        num_rows = len(enrichment_kit_values)
        num_cols = 1

        first_axis = None

        for i, (enrichment_kit_name, enrichment_kit_values) in enumerate(enrichment_kit_values):
            logging.info("enrichment_kit_name = %s", enrichment_kit_name)
            logging.info("enrichment_kit data: %s ", enrichment_kit_values)

            kwargs = {}
            if first_axis:
                kwargs['sharex'] = first_axis
            logging.info("figure.add_subplot(num_rows=%d, num_cols=%d, i=%d, **kwargs=%s)", num_rows, num_cols, i, kwargs)
            ax = figure.add_subplot(num_rows, num_cols, i + 1, **kwargs)

            self.plot_enrichment_kit(ax, enrichment_kit_name, enrichment_kit_values)
            if first_axis is None:
                first_axis = ax

    def figure(self, figure):
        figure.tight_layout()
