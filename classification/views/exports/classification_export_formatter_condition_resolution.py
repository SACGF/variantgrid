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
    PanelAppClassification, OntologySnake, OntologyService, GeneDiseaseClassification


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

    @export_column('c.HGVS')
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

                # JAMES EDIT STARTS
                for gene_symbol in vcm.classification.allele_info.gene_symbols:
                    # on the off chance there are 2 gene symbols
                    all_relationships = list(OntologySnake.get_all_term_to_gene_relationships(condition, gene))
                    # we now have a list of potential PanelApp, MONDO, GenCC relationships all of different strengths to the exact condition term
                    # find the best strength for each kind of relationship and put it in one row
                    has_direct_panel_app_relationship = True  # TODO - determine if there actually is a panel app relationship

                    if not has_direct_panel_app_relationship:
                        if all_relationships_snakes := OntologySnake.terms_for_gene_symbol(gene_symbol=gene_symbol, desired_ontology=OntologyService.MONDO, quality_q=GeneDiseaseClassification.DISPUTED):
                            all_relationships = all_relationships_snakes.leaf_relations()
                            # we now have a list of potential PanelApp, MONDO, GenCC relationships all of different strengths to various condition terms
                            # get a list of all conditions and make a row for each
                            # mark the fact that it's not a direct match
                # JAMES EDIT ENDS

                # I have not edited the below
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
