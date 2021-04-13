import hashlib
import logging

from snpdb.graphs.graphcache import CacheableGraph
from seqauto.models import QC, QCColumn
from snpdb.models import DataState


class QCColumnLineGraph(CacheableGraph):
    def __init__(self, qc_column_id, _enrichment_kit_separation, _enrichment_kit_name, percent):
        super().__init__()
        self.qc_column = QCColumn.objects.get(pk=qc_column_id)
        self.percent = bool(percent)

        logging.info("QC_column: %s ", self.qc_column)
        logging.info("Percent: %s", self.percent)

    def get_params_hash(self):
        """ This uses get_values_list rather than just hashing params as the underlying DB may have changed """
        sha1 = hashlib.sha1()
        for values in self.get_values_list():
            for value in values:
                sha1.update(str(value).encode())

        return sha1.hexdigest()

    def get_values_list(self):
        path = self.qc_column.qc_type.qc_object_path + "__" + self.qc_column.field
        data_state = self.qc_column.qc_type.qc_object_path + "__data_state"
        qs = QC.objects.filter(**{data_state: DataState.COMPLETE})
        args = ["sequencing_run", path]
        if self.percent:
            total_field = self.qc_column.qc_type.qc_object_path + "__" + self.qc_column.qc_type.total_field
            args.append(total_field)

        # Format date
        data = []
        for sequencing_run_id, path in qs.order_by("sequencing_run").values_list(*args):
            sr_date = sequencing_run_id.split("_")[0]
            data.append((sr_date, path))
        return data

    def plot(self, ax):
        PADDING = 0.5
        x = []
        y = []
        x_labels = []
        for i, values in enumerate(self.get_values_list()):
            date_string = values[0]
            value = values[1]
            x.append(i+1)
            x_labels.append(date_string)

            if self.percent:
                value *= 100.0 / values[2]
            y.append(value)

        title = self.qc_column.name
        if self.percent:
            title += " %% of %s" % self.qc_column.qc_type.total_field

        ax.set_title(title)
        if x:
            ax.plot(x, y)
            ax.set_xticks(x)
            ax.set_xticklabels(x_labels, rotation=-90)
            ax.set_xlim(1 - PADDING, x[-1] + PADDING)

    def figure(self, figure):
        figure.tight_layout()
