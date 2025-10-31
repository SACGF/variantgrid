from dataclasses import dataclass
from functools import cached_property
from typing import Iterator

from classification.models import EvidenceKeyMap
from classification.views.classification_export_utils import UsedKeyTracker, KeyValueFormatter
from classification.views.exports.classification_export_formatter_csv import CSVCellFormatting
from classification.views.exports_grouping.classification_grouping_export_filter import \
    ClassificationGroupingExportFormat, ClassificationGroupingExportFormatProperties, \
    ClassificationGroupingExportFilter, ClassificationGroupingExportFileSettings
from library.utils import delimited_row


@dataclass(frozen=True)
class CSVFormatDetails:
    pretty = False
    full_detail = False
    html_handling: CSVCellFormatting = CSVCellFormatting.PURE_TEXT


class ClassificationGroupingExportFormatterCSV(ClassificationGroupingExportFormat):

    @classmethod
    def format_properties(cls) -> ClassificationGroupingExportFormatProperties:
        return ClassificationGroupingExportFormatProperties(
            is_genome_build_relevant=True,
            http_content_type="text/csv",
            extension="csv"
        )

    def __init__(self,
                 classification_grouping_filter: ClassificationGroupingExportFilter,
                 csv_format_details: CSVFormatDetails = CSVFormatDetails()):
        self.csv_format_details = csv_format_details
        super().__init__(classification_grouping_filter)

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
            pretty=self.csv_format_details.pretty,
            cell_formatter=self.csv_format_details.html_handling.format,
            # ignore_evidence_keys=self.csv_format_details.ignore_evidence_keys,
            include_only_evidence_keys=consider_only,
            include_explains_and_notes=self.csv_format_details.full_detail
        )
        if self.csv_format_details.full_detail:
            # apparently this is significantly quicker than the attempt to use an aggregate
            for evidence in self.queryset().values_list('latest_classification_modification__published_evidence', flat=True).iterator(chunk_size=1000):
                used_keys.check_evidence(evidence)
        else:
            used_keys.check_evidence_enable_all_considered()

        # below took up to 3 minutes in Shariant test, vs 7 seconds of just iterating through the evidence twice
        # used_keys.check_evidence_qs(self.classification_filter.cms_qs)

        return used_keys

    def header(self) -> list[str]:
        return [delimited_row(self.used_keys.header(), include_new_line=False)]

    def single_row_generator(self) -> Iterator[str]:
        for row in self.queryset().iterator():
            cm = row.latest_classification_modification
            row_data = []
            row_data.extend(self.used_keys.row(cm, formatter=self.csv_format_details.html_handling.format))
            yield delimited_row(row_data, include_new_line=False)

            # def to_row(self, vcm: ClassificationModification, allele_data: AlleleData, message=None) -> str:
            #     row_data = \
            #         RowID(cm=vcm, allele_data=allele_data, date_str=self.classification_filter.date_str, message=message).to_csv(export_tweak=self._export_tweak) + \
            #         ClassificationMeta(
            #             cm=vcm,
            #             discordance_status=self.classification_filter.is_discordant(vcm),
            #             pending_clin_sig=self.grouping_utils.pending_changes_for(vcm),
            #             e_keys=self.e_keys
            #         ).to_csv(export_tweak=self._export_tweak) + \
            #         self.used_keys.row(classification_modification=vcm, formatter=self.format_details.html_handling.format)
            #
            #
            #
            #     return delimited_row(row_data, delimiter=',')
