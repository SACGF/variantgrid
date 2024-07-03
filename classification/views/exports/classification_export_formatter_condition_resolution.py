from dataclasses import dataclass
from typing import List, Optional

from django.http import HttpRequest

from classification.enums import SpecialEKeys
from classification.models import ClassificationModification
from classification.views.exports.classification_export_decorator import register_classification_exporter
from classification.views.exports.classification_export_filter import ClassificationFilter, AlleleData
from classification.views.exports.classification_export_formatter import ClassificationExportFormatter
from genes.models import GeneSymbol
from library.utils import ExportRow, export_column, delimited_row
from ontology.models import OntologyTerm, OntologyRelation, OntologyImportSource, \
    PanelAppClassification, OntologySnake, OntologyService, GeneDiseaseClassification


@dataclass(frozen=True)
class GeneSymbolCondition:
    gene_symbol: GeneSymbol
    condition: OntologyTerm


@dataclass(frozen=True)
class ClassificationConditionResolutionRow(ExportRow):
    classification: ClassificationModification
    condition: OntologyTerm
    gene_symbol: GeneSymbol
    panel_app_strength: Optional[set[PanelAppClassification]] = None
    gencc_strength: Optional[set[GeneDiseaseClassification]] = None
    mondo_strength: Optional[set[str]] = None
    all_relations: Optional[list[OntologyRelation]] = None

    def duplicate_for(self, new_classification: ClassificationModification):
        return ClassificationConditionResolutionRow(
            classification=new_classification,
            condition=self.condition,
            gene_symbol=self.gene_symbol,
            panel_app_strength=self.panel_app_strength,
            gencc_strength=self.gencc_strength,
            mondo_strength=self.mondo_strength,
            all_relations=self.all_relations
        )

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

    @export_column("Classification")
    def classification_value(self):
        return self.classification.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)

    @export_column("Somatic Clinical Significance")
    def somatic_clinical_significance_value(self):
        return self.classification.get(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE)

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

    @export_column('All PanelApp Relationships For Gene Symbol')
    def panelapp_all(self):
        if filtered_relationships := [otr for otr in self.all_relations if otr.from_import.import_source == OntologyImportSource.PANEL_APP_AU]:
            relationship_strs = [
                (rel.relationship_quality.label if rel.relationship_quality else "Green?") + f" {rel.source_term}" for rel
                in filtered_relationships]
            return ', '.join(str(relation) for relation in relationship_strs)

    @export_column('All GenCC Relationships For Gene Symbol')
    def gencc_all(self):
        if filtered_relationships := [otr for otr in self.all_relations if otr.from_import.import_source == OntologyImportSource.GENCC]:
            relationship_strs = [
                (rel.relationship_quality.label if rel.relationship_quality else "Unknown?") + f" {rel.source_term}" for rel
                in filtered_relationships]
            return ', '.join(str(relation) for relation in relationship_strs)


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
        cache: dict[GeneSymbolCondition, ClassificationConditionResolutionRow] = {}
        for vcm in allele_data.cms_regardless_of_issues:
            if (allele_info := vcm.classification.allele_info) and (gene_symbols := allele_info.gene_symbols):
                if condition_obj := vcm.classification.condition_resolution_obj:
                    if not condition_obj.is_multi_condition and condition_obj.terms and condition_obj.terms[0].ontology_service in {OntologyService.MONDO, OntologyService.OMIM}:
                        # only work for single conditions in MONDO or OMIM
                        condition_term = condition_obj.terms[0]
                        for gene_symbol in gene_symbols:

                            row: ClassificationConditionResolutionRow
                            key = GeneSymbolCondition(gene_symbol=gene_symbol, condition=condition_term)

                            # can re-use the last row as cache
                            if cached_row := cache.get(key):
                                row = cached_row.duplicate_for(vcm)
                            else:
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

                                # grab all values, don't convert to MONDO
                                direct_relationships = OntologySnake.direct_relationships_for_gene_symbol(gene_symbol)
                                row = ClassificationConditionResolutionRow(
                                    classification=vcm,
                                    condition=condition_term,
                                    gene_symbol=gene_symbol,
                                    panel_app_strength=panel_app_strength,
                                    gencc_strength=gencc_strength,
                                    mondo_strength=mondo_strength,
                                    all_relations=direct_relationships
                                )
                                cache[key] = row

                            rows.append(delimited_row(row.to_csv()))

        return rows
