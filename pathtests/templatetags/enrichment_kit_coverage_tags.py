from django.conf import settings
from django.template.library import Library

from seqauto.models import EnrichmentKit

register = Library()


@register.inclusion_tag("pathtests/tags/enrichment_kit_coverage.html")
def enrichment_kit_coverage(pathology_test_versions,
                            sorted_enrichment_kits=None):
    if sorted_enrichment_kits is None:
        sorted_enrichment_kits = EnrichmentKit.get_enrichment_kits(settings.PATHOLOGY_TEST_SORTED_ENRICHMENT_KITS)

    enrichment_kits_names = map(str, sorted_enrichment_kits)
    # kit_coverage_rows = name of path test, then enrichment_kit, missing genes and gene list for each kit
    # eg: (pathology test version, [(enrichment_kit, 0, []), (3, ['GATA2', 'MSN', 'MAP2K4'])] """
    kit_coverage_rows = []

    kit_and_gene_symbol_sets = []
    for enrichment_kit in sorted_enrichment_kits:
        if enrichment_kit.gene_list:
            gene_names_set = set(enrichment_kit.gene_list.get_gene_names())
            kit_and_gene_symbol_sets.append((enrichment_kit, gene_names_set))

    missing_data_id = 0
    for ptv in pathology_test_versions:
        ptv_gene_names = set(ptv.gene_list.get_gene_names())
        kit_data = []
        for enrichment_kit, kit_gene_names_set in kit_and_gene_symbol_sets:
            missing = ptv_gene_names - kit_gene_names_set
            missing_data_id += 1
            missing_genes = ", ".join(sorted(missing))
            kit_data.append((enrichment_kit, len(missing), missing_genes))

        row_data = [ptv, kit_data]
        kit_coverage_rows.append(row_data)

    return {"enrichment_kits_names": enrichment_kits_names, "kit_coverage_rows": kit_coverage_rows}
