import hashlib
import logging

from snpdb.graphs.graphcache import CacheableGraph
from seqauto.models import QC, DataState, QCColumn


PADDING = 0.5

class QCColumnLineGraph(CacheableGraph):
    def __init__(self, qc_column_id, enrichment_kit_separation, enrichment_kit_name, percent):
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
        #values_list = [(datetime.datetime(2016, 2, 18), 41.2), (datetime.datetime(2016, 3, 18), 109.0)]

        SEQUENCING_DATE = 'bam_file__unaligned_reads__sequencing_sample__sample_sheet__date'

        path = self.qc_column.qc_type.qc_object_path + "__" + self.qc_column.field
        data_state = self.qc_column.qc_type.qc_object_path + "__data_state"
        qs = QC.objects.filter(**{data_state: DataState.COMPLETE})
        args = [SEQUENCING_DATE, path]
        if self.percent:
            total_field = self.qc_column.qc_type.qc_object_path + "__" + self.qc_column.qc_type.total_field
            args.append(total_field)

        values_list = qs.order_by(SEQUENCING_DATE).values_list(*args)

        # Format date
        data = []
        for values in values_list:
            l = list(values)
            date = l[0]
            l[0] = date.strftime("%Y-%m-%d")
            data.append(tuple(l))
        return data

    def plot(self, ax):
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
