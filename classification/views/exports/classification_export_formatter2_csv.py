from dataclasses import dataclass
from typing import List

from django.http import HttpRequest
from lazy import lazy

from classification.models import Classification, ClassificationModification
from classification.views.classification_export_utils import UsedKeyTracker, KeyValueFormatter
from classification.views.exports.classification_export_formatter2 import ClassificationExportFormatter2
from classification.views.exports.classification_export_filter import AlleleData, ClassificationFilter
from classification.views.exports.classification_exporter import register_classification_exporter
from library.utils import delimited_row


"""
Doesn't have the ability to include errors without a bit of inefficient reworking of parent class.
Not sure I want to do that
"""


@dataclass
class FormatDetailsCSV:
    pretty: bool = False

    @staticmethod
    def from_request(request: HttpRequest) -> 'FormatDetailsCSV':
        pretty = request.query_params.get('value_format') == 'labels'
        return FormatDetailsCSV(pretty=pretty)


@register_classification_exporter("csv")
class ClassificationExportFormatter2CSV(ClassificationExportFormatter2):

    def __init__(self, filter: ClassificationFilter, format: FormatDetailsCSV):
        self.format = format
        super().__init__(filter=filter)

    @staticmethod
    def from_request(request: HttpRequest) -> 'ClassificationExportFormatter2CSV':
        return ClassificationExportFormatter2CSV(
            filter=ClassificationFilter.from_request(request),
            format=FormatDetailsCSV.from_request(request)
        )

    @lazy
    def used_keys(self) -> UsedKeyTracker:
        used_keys = UsedKeyTracker(self.filter.user, self.ekeys, KeyValueFormatter(), pretty=self.format.pretty)
        for evidence in self.filter.cms_qs().values_list('published_evidence', flat=True):
            used_keys.check_evidence(evidence)
        return used_keys

    def content_type(self) -> str:
        return "text/csv"

    def extension(self) -> str:
        return "csv"

    def header(self) -> List[str]:
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
        rows = []
        for vcm in allele_data.cms:
            rows += self.to_row(vcm)
        return rows

    def footer(self) -> List[str]:
        # FIXME need to print out all the errors here
        return []

    def to_row(self, vcm: ClassificationModification, message=None) -> str:
        vc = vcm.classification

        acmg_criteria = vcm.criteria_strength_summary(self.ekeys)
        evidence_weights = Classification.summarize_evidence_weights(vcm.evidence, self.ekeys)
        citations = ', '.join([c.ref_id() for c in vcm.citations])

        full_chgvs = vc.get_c_hgvs(genome_build=self.filter.genome_build, use_full=False)

        row = [
            vc.id,
            vc.lab.name,  # row
            vc.lab_record_id,
            vcm.share_level_enum.label,
            vcm.created.timestamp(),
            message,
            vc.allele_id,
            self.genome_build.name,
            full_chgvs,
            (vc.condition_resolution_dict or {}).get('display_text'),
            acmg_criteria,
            evidence_weights,
            citations,
            'TRUE' if self.is_discordant(vc) else 'FALSE',
        ] + self.used_keys.row(classification_modification=vcm)
        return delimited_row(row, delimiter=',')
