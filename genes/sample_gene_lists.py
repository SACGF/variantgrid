""" The cards on a sample's Genes Of Interest tab (genes/views/views.py:sample_gene_lists_tab): each SampleGeneList
    with where it came from (pipeline or the user who applied it, and the gene list's PanelApp provenance) and the
    sample's coverage of its genes (#2846). Entry point: get_sample_gene_list_cards """
from dataclasses import dataclass
from typing import Optional

from django.utils.timezone import localtime

from genes.models import (
    GeneCoverageCollection,
    GeneListCoverageSummary,
    PanelAppProvenance,
    SampleGeneList,
)
from snpdb.models import Sample


@dataclass(frozen=True)
class SampleGeneListCard:
    sample_gene_list: SampleGeneList
    origin: str
    provenance: Optional[PanelAppProvenance]
    gene_coverage_collection: Optional[GeneCoverageCollection]
    coverage: Optional[GeneListCoverageSummary]


def get_sample_gene_list_origin(sample_gene_list: SampleGeneList, provenance: Optional[PanelAppProvenance]) -> str:
    date = localtime(sample_gene_list.created).date()
    if sample_gene_list.qcgenelist_set.exists():
        return f"Pipeline (API), {date}"
    if sample_gene_list.user:
        gene_list = sample_gene_list.gene_list
        if provenance:
            source = provenance.label
        elif hasattr(gene_list, "customtextgenelist"):
            source = "custom text"
        else:
            source = f"gene list '{gene_list.name}'"
        return f"Applied by {sample_gene_list.user} {date}, {source}"
    return f"Created {date}"


def get_sample_gene_list_cards(sample: Sample) -> list[SampleGeneListCard]:
    gene_coverage_collection = GeneCoverageCollection.get_gene_coverage_for_sample(sample)
    if gene_coverage_collection and gene_coverage_collection.data_archived:
        gene_coverage_collection = None  # Only the QC tab says it has been archived
    min_coverage = sample.get_minimum_coverage_required()

    cards = []
    for sample_gene_list in sample.samplegenelist_set.select_related("gene_list", "user"):
        provenance = sample_gene_list.gene_list.get_panel_app_provenance()
        coverage = None
        if gene_coverage_collection:
            coverage = gene_coverage_collection.get_gene_lists_coverage_summary([sample_gene_list.gene_list],
                                                                               min_coverage)
        cards.append(SampleGeneListCard(sample_gene_list=sample_gene_list,
                                        origin=get_sample_gene_list_origin(sample_gene_list, provenance),
                                        provenance=provenance,
                                        gene_coverage_collection=gene_coverage_collection,
                                        coverage=coverage))
    return cards
