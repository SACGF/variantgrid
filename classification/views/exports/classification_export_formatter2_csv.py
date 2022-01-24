from dataclasses import dataclass
from io import StringIO
from typing import List, Optional

from django.http import HttpRequest
from lazy import lazy

from classification.models import Classification, ClassificationModification, EvidenceKeyMap
from classification.views.classification_export_utils import UsedKeyTracker, KeyValueFormatter
from classification.views.exports.classification_export_decorator import register_classification_exporter
from classification.views.exports.classification_export_formatter2 import ClassificationExportFormatter2
from classification.views.exports.classification_export_filter import AlleleData, ClassificationFilter
from library.utils import delimited_row, export_column, ExportRow
from snpdb.models import GenomeBuild


@dataclass
class FormatDetailsCSV:
    pretty: bool = False

    @staticmethod
    def from_request(request: HttpRequest) -> 'FormatDetailsCSV':
        pretty = request.query_params.get('value_format') == 'labels'
        return FormatDetailsCSV(pretty=pretty)


class RowID(ExportRow):

    def __init__(self, cm: ClassificationModification, genome_build: GenomeBuild, message: Optional[str] = None):
        self.cm = cm
        self.vc = cm.classification
        self.message = message
        self.genome_build = genome_build

    @export_column()
    def id(self):
        return self.vc.id

    @export_column()
    def lab(self):
        return self.vc.lab.name

    @export_column()
    def lab_record_id(self):
        return self.vc.lab_record_id

    @export_column()
    def share_level_enum(self):
        return self.cm.share_level_enum.label

    @export_column()
    def version(self):
        return self.cm.created.timestamp()

    @export_column()
    def liftover_error(self):
        return self.message

    @export_column()
    def internal_allele_id(self):
        return self.vc.allele_id

    @export_column()
    def resolved_clingen_allele_id(self):
        if allele := self.vc.allele:
            return str(allele.clingen_allele)

    @export_column()
    def target_genome_build(self):
        return self.genome_build.name

    @export_column()
    def target_variant_coordinate(self):
        try:
            if allele := self.vc.allele:
                if variant := allele.variant_for_build(genome_build=self.genome_build, best_attempt=True):
                    return str(variant)
        except ValueError:
            pass

    @export_column
    def target_c_hgvs(self):
        return self.vc.get_c_hgvs(genome_build=self.genome_build, use_full=False)


class ClassificationMeta(ExportRow):

    def __init__(self, cm: ClassificationModification, discordant: bool, e_keys: EvidenceKeyMap):
        self.cm = cm
        self.vc = cm.classification
        self.discordant = discordant
        self.e_keys = e_keys

    @export_column()
    def resolved_condition(self):
        return (self.vc.condition_resolution_dict or {}).get('display_text')

    @export_column()
    def acmg_criteria(self):
        return self.cm.criteria_strength_summary(self.e_keys)

    @export_column()
    def evidence_weights(self):
        return Classification.summarize_evidence_weights(self.cm.evidence, self.e_keys)

    @export_column()
    def citations(self):
        return ', '.join([c.ref_id() for c in self.cm.citations])

    @export_column()
    def is_discordant(self):
        return 'TRUE' if self.discordant else 'FALSE'


@register_classification_exporter("csv")
class ClassificationExportFormatter2CSV(ClassificationExportFormatter2):

    def __init__(self, classification_filter: ClassificationFilter, format_details: FormatDetailsCSV):
        self.format_details = format_details
        self.errors_io: Optional[StringIO] = None
        self.e_keys = EvidenceKeyMap.cached()
        super().__init__(classification_filter=classification_filter)

    @staticmethod
    def from_request(request: HttpRequest) -> 'ClassificationExportFormatter2CSV':
        return ClassificationExportFormatter2CSV(
            classification_filter=ClassificationFilter.from_request(request),
            format_details=FormatDetailsCSV.from_request(request)
        )

    @lazy
    def used_keys(self) -> UsedKeyTracker:
        used_keys = UsedKeyTracker(self.classification_filter.user, self.e_keys, KeyValueFormatter(), pretty=self.format_details.pretty)
        for evidence in self.classification_filter.cms_qs().values_list('published_evidence', flat=True):
            used_keys.check_evidence(evidence)
        return used_keys

    def content_type(self) -> str:
        return "text/csv"

    def extension(self) -> str:
        return "csv"

    def header(self) -> List[str]:

        self.errors_io = StringIO()
        header = RowID.csv_header() + ClassificationMeta.csv_header() + self.used_keys.header()
        return [delimited_row(header, delimiter=',')]

    def row(self, allele_data: AlleleData) -> List[str]:
        # record error to report them in the footer
        if issues := allele_data.issues:
            for issue in issues:
                self.errors_io.writelines(self.to_row(issue.classification, message=issue.message))
        rows = []
        for vcm in allele_data.cms:
            rows.append(self.to_row(vcm))

        return rows

    def footer(self) -> List[str]:
        footer_content = self.errors_io.getvalue()
        self.errors_io.close()
        if footer_content:
            return [footer_content]
        else:
            return []

    def to_row(self, vcm: ClassificationModification, message=None) -> str:
        row_data = \
            RowID(cm=vcm, genome_build=self.classification_filter.genome_build, message=message).to_csv() + \
            ClassificationMeta(cm=vcm, discordant=self.classification_filter.is_discordant(vcm), e_keys=self.e_keys).to_csv() + \
            self.used_keys.row(classification_modification=vcm)

        return delimited_row(row_data, delimiter=',')
