// @ts-check
// genes/templates/genes/genes.html
let gridExtraFilters = {};

// The page level filters (release, and the OMIM terms shortcuts) - the CSV download picks these
// up too, as it goes out with the table's own ajax params
function genesDatatableFilter(data) {
    Object.assign(data, gridExtraFilters);
    data['gene_annotation_release_id'] = $("#id_release").val();
}

function filterGrid(extra_filters) {
    gridExtraFilters = extra_filters || {};
    $("#gene-annotation-versions-grid").DataTable().ajax.reload();
}

$(document).ready(function() {
    $("#id_release").change(function() {
        filterGrid();
    });

    $('#id_gene_symbol').change(function() {
        const geneSymbol = $("#id_gene_symbol").val();
        window.location = Urls.view_gene_symbol(geneSymbol);
    });
});
