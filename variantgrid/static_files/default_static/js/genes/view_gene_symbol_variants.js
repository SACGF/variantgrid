// @ts-check
// genes/templates/genes/view_gene_symbol.html
/* global geneVariantsGridParams:writable */ // genes/templates/genes/view_gene_symbol.html
/* global setupTagCountsSummary, clearTagCountsSummary */ // tag_counts_summary.js
// Called by DataTables before each ajax call - the grid class reads them as extra_filters
function geneVariantsDatatableFilter(data) {
    data["extra_filters"] = JSON.stringify(geneVariantsGridParams);
    data["show_clinvar"] = $("#gene-variants-show-clinvar").is(":checked");
}

function reloadGeneVariantsTable() {
    const table = $("#gene-variants-datatable");
    if ($.fn.DataTable.isDataTable(table)) {
        table.DataTable().ajax.reload();
    }
}

function clearGeneVariantsGrid() {
    geneVariantsGridParams = {};
    clearTagCountsSummary(GENE_VARIANTS_TAG_COUNTS);
    reloadGeneVariantsTable();
    $("#gene-variants-grid-filtering-message").empty().parent().hide();
}

function geneVariantsTagsChanged(tagIds) {
    geneVariantsGridParams = tagIds.length ? {"tags": tagIds} : {};
    reloadGeneVariantsTable();
}

const GENE_VARIANTS_TAG_COUNTS = "[data-tag-counts-summary]";

$(document).ready(() => {
    setupTagCountsSummary(GENE_VARIANTS_TAG_COUNTS, geneVariantsTagsChanged);
    $("#gene-variants-show-clinvar").change(reloadGeneVariantsTable);
});
