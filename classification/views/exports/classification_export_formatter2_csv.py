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
from library.utils import delimited_row


@dataclass
class FormatDetailsCSV:
    pretty: bool = False

    @staticmethod
    def from_request(request: HttpRequest) -> 'FormatDetailsCSV':
        pretty = request.query_params.get('value_format') == 'labels'
        return FormatDetailsCSV(pretty=pretty)


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
        header = [
                 'id',
                 'lab',
                 'lab_record_id',
                 'share_level',
                 'version',
                 'liftover_error',
                 'internal_allele_id',
                 'target_genome_build',
                 'target_c_hgvs',
                 'resolved_condition',
                 'acmg_criteria',
                 'evidence_weights',
                 'citations',
                 'in_discordance'
             ] + self.used_keys.header()
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
        vc = vcm.classification

        acmg_criteria = vcm.criteria_strength_summary(self.e_keys)
        evidence_weights = Classification.summarize_evidence_weights(vcm.evidence, self.e_keys)
        citations = ', '.join([c.ref_id() for c in vcm.citations])

        full_chgvs = vc.get_c_hgvs(genome_build=self.classification_filter.genome_build, use_full=False)

        row = [
            vc.id,
            vc.lab.name,  # row
            vc.lab_record_id,
            vcm.share_level_enum.label,
            vcm.created.timestamp(),
            message,
            vc.allele_id,
            self.classification_filter.genome_build.name,
            full_chgvs,
            (vc.condition_resolution_dict or {}).get('display_text'),
            acmg_criteria,
            evidence_weights,
            citations,
            'TRUE' if self.classification_filter.is_discordant(vcm) else 'FALSE',
        ] + self.used_keys.row(classification_modification=vcm)
        return delimited_row(row, delimiter=',')
