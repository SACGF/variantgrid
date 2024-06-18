from dataclasses import dataclass
from typing import List, Optional

from django.http import HttpRequest
from classification.models import ClassificationModification
from classification.views.exports.classification_export_decorator import register_classification_exporter
from classification.views.exports.classification_export_filter import ClassificationFilter, AlleleData
from classification.views.exports.classification_export_formatter import ClassificationExportFormatter
from library.utils import ExportRow, export_column, delimited_row
from ontology.models import OntologyTerm, OntologyRelation, OntologyImportSource, \
    PanelAppClassification, OntologySnake, OntologyService, GeneDiseaseClassification, \
    ONTOLOGY_RELATIONSHIP_MINIMUM_QUALITY_FILTER


@dataclass(frozen=True)
class ClassificationConditionResolutionRow(ExportRow):
    classification: ClassificationModification
    condition: OntologyTerm
    gene_symbol: str
    panel_app_strength: Optional[set[PanelAppClassification]]
    gencc_strength: Optional[set[GeneDiseaseClassification]]
    mondo_strength: Optional[set[str]]
    all_relations: Optional[set[OntologyTerm]] = None

    @property
    def cm(self):
        return self.classification

    @property
    def vc(self):
        return self.classification.classification

    @export_column('Classification ID')
    def classification_id(self):
        return self.vc.id

    @export_column('Lab')
    def lab(self):
        return self.vc.lab.name

    @export_column('c.HGVS')
    def c_hgvs(self):
        return self.vc.c_parts.full_c_hgvs

    @export_column('Allele Origin')
    def allele_origin(self):
        return self.vc.allele_origin_bucket

    @export_column('Provided Condition')
    def condition(self):
        return self.condition

    @export_column('Date Curated')
    def date_curated(self):
        return self.cm.curated_date.date()

    @export_column('Gene Symbol')
    def gene_symbol(self):
        return self.gene_symbol

    @export_column('Panel App Strength')
    def panel_app_strength(self):
        if self.panel_app_strength:
            return ', '.join(relation.label for relation in self.panel_app_strength)

    @export_column('Gencc Strength')
    def gencc_strength(self):
        if self.gencc_strength:
            return ', '.join(relation.label for relation in self.gencc_strength)

    @export_column('MONDO Strength')
    def mondo_strength(self):
        if self.mondo_strength:
            return ', '.join(self.mondo_strength)

    @export_column('All Panel App Relationships For Gene Symbol')
    def matching_relations(self):
        if self.all_relations:
            return ', '.join(str(relation) for relation in self.all_relations)


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
            if condition_obj := vcm.classification.condition_resolution_obj:
                if not condition_obj.is_multi_condition and condition_obj.terms and condition_obj.terms[0].ontology_service in {OntologyService.MONDO, OntologyService.OMIM}:
                    # only work for single conditions in MONDO or OMIM
                    condition_term = condition_obj.terms[0]
                    for gene_symbol in vcm.classification.allele_info.gene_symbols:
                        # on the off chance there are 2 gene symbols, we can make results for each gene symbol individually
                        all_relationships = []
                        for snake in OntologySnake.get_all_term_to_gene_relationships(condition_term, gene_symbol):
                            if snake_relations := snake.get_import_relations:
                                all_relationships.append(snake_relations)

                        # we now have a list of potential PanelApp, MONDO, GenCC relationships all of different strengths to the exact condition term
                        panel_app_strength, gencc_strength, mondo_strength = set(), set(), set()
                        for rel in all_relationships:
                            source = rel.from_import.import_source
                            if source in OntologyImportSource.PANEL_APP_AU:
                                panel_app_strength.add(rel.relationship_quality)
                            elif source == OntologyImportSource.MONDO:
                                mondo_strength.add(rel.relation)
                            elif source == OntologyImportSource.GENCC:
                                gencc_strength.add(rel.relationship_quality)

                        # filter to only be panel app
                        all_gene_symbol_relationships = OntologySnake.terms_for_gene_symbol(gene_symbol=gene_symbol, desired_ontology=OntologyService.MONDO, quality_filter=ONTOLOGY_RELATIONSHIP_MINIMUM_QUALITY_FILTER)
                        panel_app_snakes = [os for os in all_gene_symbol_relationships if os.get_import_relations.from_import.import_source == OntologyImportSource.PANEL_APP_AU]
                        panel_app_relationships = [os.get_import_relations for os in panel_app_snakes]
                        panel_app_relationship_strs = [(rel.relationship_quality.label if rel.relationship_quality else "Green?") + f" {rel.source_term}" for rel in panel_app_relationships]

                        row = ClassificationConditionResolutionRow(vcm, condition_term, gene_symbol,
                                                                   panel_app_strength, gencc_strength,
                                                                   mondo_strength, panel_app_relationship_strs)

                        rows.append(delimited_row(row.to_csv()))

        return rows
