from dataclasses import dataclass
from functools import cached_property
from typing import Optional, List, Iterator

from django.http import HttpRequest

from classification.enums import SpecialEKeys
from classification.models import ClassificationModification, EvidenceKeyMap
from classification.views.exports.classification_export_decorator import register_classification_exporter
from classification.views.exports.classification_export_filter import ClassificationFilter, AlleleData, \
    ClassificationIssue
from classification.views.exports.classification_export_formatter import ClassificationExportFormatter2
from flags.models import FlagComment, Flag
from library.django_utils import get_url_from_view_path
from library.utils import export_column, ExportDataType, ExportRow, delimited_row


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

    @export_column("Classification Last Updated", data_type=ExportDataType.datetime)
    def classification_last_updated(self):
        return self.cm.created

    @export_column("Flag Last Updated", data_type=ExportDataType.datetime)
    def flag_last_updated(self):
        if fc := self.problem.flag_comment:
            return fc.created

    @export_column("Problem Code")
    def problem_code(self):
        return self.problem.code

    @export_column("Problem")
    def problem(self):
        return self.problem.message


@register_classification_exporter("flags")
class ClassificationExportFormatterFlags(ClassificationExportFormatter2):

    CODE_SUBSTITUTE = {'mandatory': 'Missing Value'}

    def __init__(self, classification_filter: ClassificationFilter):
        self.first_row = True
        self.e_keys = EvidenceKeyMap.cached()
        super().__init__(classification_filter=classification_filter)

    @classmethod
    def from_request(cls, request: HttpRequest) -> 'ClassificationExportFormatterFlags':
        return ClassificationExportFormatterFlags(
            classification_filter=ClassificationFilter.from_request(request)
        )

    def header(self) -> List[str]:
        return [delimited_row(ProblemRow.csv_header(), delimiter=',')]

    def iterate_problems(self, ci: ClassificationIssue) -> Iterator[Problem]:
        vcm = ci.classification
        evidence = vcm.evidence
        # validation messages
        for key, blob in evidence.items():
            validations = blob.get('validation')
            if validations:
                for validation in validations:
                    key_label = self.e_keys.get(key).pretty_label
                    code = validation.get('code', '')
                    code = ClassificationExportFormatterFlags.CODE_SUBSTITUTE.get(code, code)
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

        # allele validation
        if not ci.validation_include:
            if (allele_info := ci.classification.classification.allele_info) and (validation := allele_info.latest_validation):
                errors = [row for row in validation.validation_tags_list if row.severity == 'E']
                for error in errors:
                    yield Problem(code="Variant Matching", message=f"{error.category_pretty} {error.field_pretty}")

    def row(self, data: AlleleData) -> List[str]:
        output: List[str] = []
        for row in data.all_cms:
            if row.withdrawn:
                continue

            for problem in self.iterate_problems(row):
                output.append( delimited_row(ProblemRow(row.classification, problem).to_csv(), ",") )

        return output

    def extension(self) -> str:
        return "csv"

    def content_type(self) -> str:
        return "text/csv"
