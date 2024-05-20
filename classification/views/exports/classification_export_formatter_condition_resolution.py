from typing import List

from django.http import HttpRequest
from django.urls.base import reverse

from classification.enums import SpecialEKeys
from classification.models import ClassificationModification, Classification
from classification.views.exports.classification_export_decorator import register_classification_exporter
from classification.views.exports.classification_export_filter import ClassificationFilter, AlleleData
from classification.views.exports.classification_export_formatter import ClassificationExportFormatter
from library.django_utils import get_url_from_view_path
from library.utils import ExportRow, export_column, delimited_row
from ontology.models import OntologyTerm, OntologyTermRelation, OntologyRelation, OntologyImportSource, \
    PanelAppClassification


class ClassificationConditionResolutionRow(ExportRow):

    def __init__(self, classificationid, condition, gene_symbol, strongest_classification, panel_app, other):
        self.cm = classificationid
        self.vc = self.cm.classification
        self.condition = condition
        self.gene_symbol = gene_symbol
        self.strongest_classification = strongest_classification
        self.panel_app = panel_app
        self.other = other

    @export_column('Classification ID')
    def classification_id(self):
        return self.vc.id

    @export_column('C.HGVS')
    def c_hgvs(self):
        return self.vc.c_parts.full_c_hgvs

    @export_column('Condition')
    def condition(self):
        return self.condition.get('display_text')

    @export_column('Gene Symbol')
    def gene_symbol(self):
        return self.vc.allele_info.gene_symbols[0].symbol if self.vc.allele_info.gene_symbols else None

    @export_column('Strongest Classification')
    def strongest_classification(self):
        return self.strongest_classification

    @export_column('Panel App Relation')
    def panel_app_resolution(self):
        return self.panel_app

    @export_column('Other')
    def other(self):
        return self.other


@register_classification_exporter("condition_resolution")
class ClassificationExportFormatterConditionResolution(ClassificationExportFormatter):

    @classmethod
    def from_request(cls, request: HttpRequest) -> 'ClassificationExportFormatterConditionResolution':
        return ClassificationExportFormatterConditionResolution(
            classification_filter=ClassificationFilter.from_request(request),
        )

    def content_type(self) -> str:
        return "text/csv"

    def extension(self) -> str:
        return "csv"

    def header(self) -> List[str]:
        return [delimited_row(ClassificationConditionResolutionRow.csv_header())]

    def footer(self) -> List[str]:
        return []

    def row(self, allele_data: AlleleData) -> list[str]:
        rows: list[str] = []
        for vcm in allele_data.cms:
            if condition := vcm.classification.condition_resolution:
                gene = vcm.classification.allele_info.gene_symbols[
                    0].symbol if vcm.classification.allele_info.gene_symbols else None
                panel_app_data = set()
                other_relations = set()
                strongest_classification = ''
                ontology_term = OntologyTerm.get_gene_symbol(gene)
                panel_app_relations = OntologyTermRelation.objects.filter(dest_term_id=ontology_term.id)
                for rel in panel_app_relations:
                    if rel.source_term.name:
                        if rel.from_import.import_source == OntologyImportSource.PANEL_APP_AU:
                            if rel.extra.get('strongest_classification', '') != 'Expert Review Green':
                                strongest_classification = rel.extra.get('strongest_classification', '')
                                panel_app_data.add(rel.source_term.name)
                        elif rel.from_import.import_source in [OntologyImportSource.MONDO, OntologyImportSource.GENCC]:
                            other_relations.add(f'({rel.from_import.import_source}): {rel.source_term.name} ')

                panel_app_data_str = ','.join(panel_app_data)
                others_str = ','.join(other_relations)

                row = ClassificationConditionResolutionRow(vcm, condition, gene, strongest_classification,
                                                           panel_app_data_str, others_str)
                rows.append(delimited_row(row.to_csv()))

        return rows
