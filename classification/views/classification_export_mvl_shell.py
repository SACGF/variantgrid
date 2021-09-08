import csv
import io

from classification.views.classification_export_mvl import ExportFormatterMVL
from classification.views.classification_export_utils import AlleleGroup


class ExportFormatterMVLShell(ExportFormatterMVL):

    def row(self, group: AlleleGroup) -> str:
        out = io.StringIO()
        writer = csv.writer(out, delimiter='\t')

        for c_parts, vcms_w_chgvs in group.iter_c_hgvs_versionless_transcripts():

            transcript = c_parts.transcript
            c_hgvs = c_parts.raw_c

            self.row_count += 1
            writer.writerow([transcript, c_hgvs, "VOUS", "This is a test", "This is a test"])

        return out.getvalue()

    def filename(self) -> str:
        return self.generate_filename(suffix='mvl_shell', extension='tsv')
