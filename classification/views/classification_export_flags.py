import csv
import io
from dataclasses import dataclass
from typing import Optional

from classification.enums.classification_enums import SpecialEKeys
from classification.models.classification import ClassificationModification
from classification.views.classification_export_utils import ExportFormatter
from flags.models.models import Flag, FlagComment
from library.django_utils import get_url_from_view_path
from library.utils import ExportRow, export_column


@dataclass(frozen=True)
class Problem:
    code: str
    message: str
    flag_comment: Optional[FlagComment] = None


@dataclass(frozen=True)
class ProblemRow(ExportRow):
    cm: ClassificationModification
    problem: Problem

    @property
    def classification(self):
        return self.cm.classification

    @export_column("URL")
    def url(self):
        return get_url_from_view_path(self.classification.get_absolute_url())

    @export_column("Lab")
    def lab(self):
        return str(self.classification.lab)

    @export_column("Lab Record ID")
    def lab_record_id(self):
        return self.classification.lab_record_id

    @export_column("Imported c.HGVS")
    def exported_c_hgvs(self):
        return self.cm.get(SpecialEKeys.C_HGVS)

    @export_column("Classification Last Updated", format={"tz": "default"})
    def classification_last_updated(self):
        return self.cm.created

    @export_column("Flag Last Updated", format={"tz": "default"})
    def flag_last_updated(self):
        if fc := self.problem.flag_comment:
            return fc.created

    @export_column("Problem Code")
    def problem_code(self):
        return self.problem.code

    @export_column("Problem")
    def problem(self):
        return self.problem.message


class ExportFormatterFlags(ExportFormatter):

    CODE_SUBSTITUTE = {'mandatory': 'Missing value'}

    @property
    def filter_out_not_lifted_over_to_desired(self) -> bool:
        return False

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def header(self) -> str:
        return ExportFormatter.write_single_row(ProblemRow.csv_header(), delimiter=',')

    def iterate_problems(self, vcm: ClassificationModification):
        evidence = vcm.evidence
        # validation messages
        for key, blob in evidence.items():
            validations = blob.get('validation')
            if validations:
                for validation in validations:
                    key_label = self.ekeys.get(key).pretty_label
                    code = validation.get('code', '')
                    code = ExportFormatterFlags.CODE_SUBSTITUTE.get(code, code)
                    code = (code[0].upper() + code[1:]).replace('_', ' ')
                    message = validation.get('message')
                    yield Problem(code=code, message=f'{key_label}: {message}')
        # flags
        open_flags = vcm.classification.flag_collection_safe.flags(only_open=True)\
            .select_related('flag_type', 'resolution')
        flag: Flag
        for flag in open_flags:
            flag_type = flag.flag_type
            parts = [flag_type.label]
            if flag.resolution_id != 'open':
                parts.append(flag.resolution.label)
            last_comment = FlagComment.last(flag)
            if last_comment and last_comment.text:
                parts.append(last_comment.text)

            yield Problem(code=flag_type.label, message='\n'.join(parts), flag_comment=last_comment)

    def row(self, vcm: ClassificationModification):
        out = io.StringIO()
        writer = csv.writer(out, delimiter=',')

        for problem in self.iterate_problems(vcm):
            self.row_count += 1
            writer.writerow(ProblemRow(vcm, problem).to_csv())

        return out.getvalue()

    def prepare_groups(self):
        pass

    def row_iterator(self):
        return self.qs

    def filename(self) -> str:
        return self.generate_filename(include_genome_build=False, suffix='issues')
