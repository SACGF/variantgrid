import itertools
from dataclasses import dataclass
from functools import cached_property
from typing import Iterator, Any

from classification.models import EvidenceKeyMap, ClassificationGrouping
from classification.views.classification_export_utils import UsedKeyTracker, KeyValueFormatter
from classification.views.exports.classification_export_formatter_csv import CSVCellFormatting
from classification.views.exports_grouping.classification_grouping_export_filter import \
    ClassificationGroupingExportFormat, ClassificationGroupingExportFormatProperties, \
    ClassificationGroupingExportFilter, ClassificationGroupingExportFileSettings
from library.django_utils import get_url_from_view_path
from library.utils import delimited_row, ExportRow, export_column, ExportTweak
from snpdb.models import GenomeBuild


@dataclass(frozen=True)
class CSVFormatDetails:
    full_detail = False
    html_handling: CSVCellFormatting = CSVCellFormatting.PURE_TEXT


@dataclass(frozen=True)
class CSVNonEvidence(ExportRow):
    classification_grouping: ClassificationGrouping
    date_str: str

    @export_column(label="URL")
    def url(self):
        return get_url_from_view_path(self.classification_grouping.get_absolute_url()) + "?refer=csv&seen=" + self.date_str

    @export_column(label="Lab")
    def lab(self):
        return str(self.classification_grouping.lab)

    @export_column(label="Resolved GRCh37 c.HGVS", categories={"GRCh37": True})
    def grch37_hgvs(self):
        if allele_info := self.classification_grouping.latest_allele_info:
            if grch := allele_info.grch37:
                return grch.c_hgvs
        return None

    @export_column(label="Resolved GRCh38 c.HGVS", categories={"GRCh38": True})
    def grch38_hgvs(self):
        if allele_info := self.classification_grouping.latest_allele_info:
            if grch := allele_info.grch38:
                return grch.c_hgvs
        return None

    @export_column(label="Resolved ClinGen Allele")
    def clingen_allele(self):
        return self.classification_grouping.allele.clingen_allele

    @export_column(label="Classification Count")
    def record_count(self):
        return self.classification_grouping.classification_count


class ClassificationGroupingExportFormatterCSV(ClassificationGroupingExportFormat):

    @classmethod
    def format_properties(cls) -> ClassificationGroupingExportFormatProperties:
        return ClassificationGroupingExportFormatProperties(
            http_content_type="text/csv",
            extension="csv"
        )

    def __init__(self,
                 classification_grouping_filter: ClassificationGroupingExportFilter,
                 csv_format_details: CSVFormatDetails = CSVFormatDetails()):
        self.csv_format_details = csv_format_details
        super().__init__(classification_grouping_filter)

    @cached_property
    def export_tweak(self):
        return ExportTweak(
            categories={
                "GRCh37": GenomeBuild.grch37().enabled,
                "GRCh38": GenomeBuild.grch38().enabled
            }
        )

    @cached_property
    def used_keys(self) -> UsedKeyTracker:
        e_keys = EvidenceKeyMap.cached()
        consider_only = None
        if not self.csv_format_details.full_detail:
            consider_only = [e_key.key for e_key in e_keys.vital()]

        used_keys = UsedKeyTracker(
            user=self.classification_grouping_filter.user,
            ekeys=e_keys,
            key_value_formatter=KeyValueFormatter(),
            pretty_headers=True,
            pretty_values=False,
            cell_formatter=self.csv_format_details.html_handling.format,
            # ignore_evidence_keys=self.csv_format_details.ignore_evidence_keys,
            include_only_evidence_keys=consider_only,
            include_explains_and_notes=self.csv_format_details.full_detail
        )
        if self.csv_format_details.full_detail:
            # This is significantly quicker than the attempt to use an aggregate
            for evidence in self.queryset().values_list('latest_classification_modification__published_evidence', flat=True).iterator(chunk_size=1000):
                used_keys.check_evidence(evidence)
        else:
            used_keys.check_evidence_enable_all_considered()

        return used_keys

    def header(self) -> list[str]:
        return [delimited_row(CSVNonEvidence.csv_header(export_tweak=self.export_tweak) + self.used_keys.header(), include_new_line=False)]

    def single_row_generator(self) -> Iterator[str]:
        for row in self.queryset().iterator():
            cm = row.latest_classification_modification
            row_data = []
            row_data.extend(CSVNonEvidence(row, date_str=self.classification_grouping_filter.date_str).to_csv(export_tweak=self.export_tweak))
            row_data.extend(self.used_keys.row(cm, formatter=self.csv_format_details.html_handling.format))
            yield delimited_row(row_data, include_new_line=False)
