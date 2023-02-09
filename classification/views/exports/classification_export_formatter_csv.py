from dataclasses import dataclass
from functools import cached_property
from typing import List, Optional, Dict, Set

from django.http import HttpRequest

from classification.models import Classification, ClassificationModification, EvidenceKeyMap
from classification.models.classification_groups import ClassificationGroupUtils
from classification.views.classification_export_utils import UsedKeyTracker, KeyValueFormatter
from classification.views.exports.classification_export_decorator import register_classification_exporter
from classification.views.exports.classification_export_filter import AlleleData, ClassificationFilter, \
    DiscordanceReportStatus
from classification.views.exports.classification_export_formatter import ClassificationExportFormatter
from classification.views.exports.classification_export_utils import CitationCounter
from library.log_utils import AdminNotificationBuilder
from library.utils import delimited_row, export_column, ExportRow, ExportDataType, DebugTimer
from snpdb.models import GenomeBuild


@dataclass(frozen=True)
class FormatDetailsCSV:
    # format evidence keys to nice human labels or leave as raw codes easier handled by code
    pretty: bool = False

    # include the explain keys (text that lets labs explain their process), notes (human entered text) are always included
    include_explains: bool = False

    # exclude fields that change between environments, makes it easier to compare changes if the same data is in both environments
    exclude_transient: bool = False

    @staticmethod
    def from_request(request: HttpRequest) -> 'FormatDetailsCSV':
        pretty = request.query_params.get('value_format') == 'labels'
        include_explains = request.query_params.get('include_explains') == 'true'
        exclude_transient = request.query_params.get('exclude_transient') == 'true'
        return FormatDetailsCSV(pretty=pretty, include_explains=include_explains, exclude_transient=exclude_transient)

    @property
    def ignore_evidence_keys(self) -> Set[str]:
        if self.exclude_transient:
            return {
                "owner",
                "source_id",
                "source_data",
            }
        else:
            return set()


class RowID(ExportRow):

    def __init__(self, cm: ClassificationModification, allele_data: AlleleData, message: Optional[str] = None):
        self.cm = cm
        self.vc = cm.classification
        self.message = message
        self.allele_data = allele_data

    @property
    def genome_build(self) -> GenomeBuild:
        return self.allele_data.genome_build

    @export_column(categories={"transient": True})
    def id(self):
        return self.vc.id

    @export_column()
    def lab(self):
        return str(self.vc.lab)

    @export_column()
    def lab_record_id(self):
        return self.vc.lab_record_id

    @export_column(categories={"transient": True}, data_type=ExportDataType.datetime)
    def server_created_date(self):
        return self.vc.created

    @export_column(categories={"transient": True})
    def share_level(self):
        return self.cm.share_level_enum.label

    @export_column(categories={"transient": True})
    def version(self):
        return self.cm.created.timestamp()

    @export_column()
    def liftover_error(self):
        return self.message

    @export_column(categories={"transient": True})
    def internal_allele_id(self):
        return self.allele_data.allele_id

    @export_column()
    def resolved_clingen_allele_id(self):
        if allele := self.allele_data.allele:
            return str(allele.clingen_allele)

    @export_column()
    def target_genome_build(self):
        return self.genome_build.name

    @export_column()
    def target_c_hgvs(self):
        if c_hgvs := self.vc.get_c_hgvs(self.genome_build):
            return c_hgvs

    @export_column()
    def target_variant_coordinate(self):
        if variant := self.allele_data.variant:
            return str(variant)


class ClassificationMeta(ExportRow):
    """
    Deals with the static columns - to be followed by all the evidence key columns
    """
    def __init__(
            self,
            cm: ClassificationModification,
            discordance_status: DiscordanceReportStatus,
            pending_clin_sig: Optional[str],
            e_keys: EvidenceKeyMap):
        self.cm = cm
        self.discordance_status = discordance_status
        self.pending_clin_sig = pending_clin_sig
        self.e_keys = e_keys

    @property
    def vc(self) -> Classification:
        return self.cm.classification

    @export_column(categories={"transient": True})
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
        cc = CitationCounter()
        cc.reference_citations(self.cm)
        return ', '.join(cc.citation_ids())
        # return ', '.join([c.ref_id() for c in sorted(set(cc.citations()), key=lambda c:c.sort_key)])

    @export_column(categories={"transient": True})
    def discordance_status(self):
        if not self.discordance_status:
            return ''
        elif self.discordance_status == DiscordanceReportStatus.ON_GOING:
            return 'active discordance'
        elif self.discordance_status == DiscordanceReportStatus.CONTINUED:
            return 'continued discordance'
        elif self.discordance_status == DiscordanceReportStatus.PENDING_CONCORDANCE:
            return 'pending concordance'
        else:
            return self.discordance_status

    @export_column(categories={"transient": True, "pending_changes": True})
    def pending_clinical_significance(self):
        return self.pending_clin_sig


@register_classification_exporter("csv")
class ClassificationExportFormatterCSV(ClassificationExportFormatter):

    def __init__(self, classification_filter: ClassificationFilter, format_details: FormatDetailsCSV):
        self.format_details = format_details
        self.error_rows: List[str] = []
        self.e_keys = EvidenceKeyMap.cached()
        self.grouping_utils = ClassificationGroupUtils()
        super().__init__(classification_filter=classification_filter)

    @classmethod
    def from_request(cls, request: HttpRequest) -> 'ClassificationExportFormatterCSV':
        return ClassificationExportFormatterCSV(
            classification_filter=ClassificationFilter.from_request(request),
            format_details=FormatDetailsCSV.from_request(request)
        )

    @cached_property
    def used_keys(self) -> UsedKeyTracker:

        used_keys = UsedKeyTracker(
            self.classification_filter.user,
            self.e_keys, KeyValueFormatter(),
            pretty=self.format_details.pretty,
            include_explains=self.format_details.include_explains,
            ignore_evidence_keys=self.format_details.ignore_evidence_keys
        )
        # apparently this is signficantly quicker than the attempt to use an aggregate
        for evidence in self.classification_filter.cms_qs.values_list('published_evidence', flat=True).iterator(chunk_size=1000):
           used_keys.check_evidence(evidence)
        # below took up to 3 minutes in Shariant test, vs 7 seconds of just iterating through the evidence twice
        # used_keys.check_evidence_qs(self.classification_filter.cms_qs)

        return used_keys

    def content_type(self) -> str:
        return "text/csv"

    def extension(self) -> str:
        return "csv"

    def header(self) -> List[str]:
        self.error_rows = []
        header = RowID.csv_header(self._categories) + ClassificationMeta.csv_header(self._categories) + self.used_keys.header()
        return [delimited_row(header, delimiter=',')]

    def row(self, allele_data: AlleleData) -> List[str]:
        # record error to report them in the footer
        if issues := allele_data.issues:
            for issue in issues:
                if not issue.withdrawn:
                    self.error_rows.append(self.to_row(issue.classification,  allele_data=allele_data, message=issue.message))
        rows: List[str] = []
        for vcm in allele_data.cms:
            rows.append(self.to_row(vcm, allele_data=allele_data))

        return rows

    def footer(self) -> List[str]:
        return self.error_rows

    @cached_property
    def _categories(self) -> Optional[Dict]:
        categories = {}
        if self.format_details.exclude_transient:
            categories["transient"] = None
        if not self.grouping_utils.any_pending_changes:
            categories["pending_changes"] = None
        return categories

    def to_row(self, vcm: ClassificationModification, allele_data: AlleleData, message=None) -> str:
        row_data = \
            RowID(cm=vcm, allele_data=allele_data, message=message).to_csv(categories=self._categories) + \
            ClassificationMeta(
                cm=vcm,
                discordance_status=self.classification_filter.is_discordant(vcm),
                pending_clin_sig=self.grouping_utils.pending_changes_for(vcm),
                e_keys=self.e_keys
            ).to_csv(categories=self._categories) + \
            self.used_keys.row(classification_modification=vcm)

        return delimited_row(row_data, delimiter=',')
