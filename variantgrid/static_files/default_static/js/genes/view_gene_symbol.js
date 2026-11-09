// @ts-check
// genes/templates/genes/view_gene_symbol.html
/* global classificationGroupingRedraw */ // classification/classification_groupings.js
/* global setupTagCountsSummary, clearTagCountsSummary */ // tag_counts_summary.js

let geneSymbolPage = null;  // options passed to initViewGeneSymbol
let geneClassificationsGridParams = {};
let geneVariantsGridParams = {};

const GENE_VARIANTS_TAG_COUNTS = "[data-tag-counts-summary]";

/**
 * @param {{geneSymbol: string, genomeBuildName: string, quickLinksAllowList: Array, quickLinksBlockList: Array,
 *          showClassificationsHotspotGraph: boolean, showHotspotGraph: boolean, geneCoverageUrl: ?string,
 *          ideogram: ?{contig: Object, contigUrl: string, cytobands: Array, marker: Object}}} options
 */
function initViewGeneSymbol(options) {
    geneSymbolPage = options;
    $(document).ready(setupViewGeneSymbol);
}

function setupViewGeneSymbol() {
    EKeys.load().then(eKeys => {
        const vcLinks = new VCLinks(eKeys);
        let links = vcLinks.generateLinks({gene_symbol: geneSymbolPage.geneSymbol});
        links = vcLinks.filterLinks(links, false, geneSymbolPage.quickLinksAllowList, geneSymbolPage.quickLinksBlockList);
        const geneAnchors = links.filter(link => link.geneLink).map(link => link.asAnchor("bootstrap"));
        $('#gene-quick-links').html($("<ul>", {class: "list-group", html: geneAnchors}));
    });

    if (geneSymbolPage.ideogram) {
        drawGeneSymbolIdeogram(geneSymbolPage.ideogram);
    }
    if (geneSymbolPage.showClassificationsHotspotGraph) {
        loadClassificationsHotspotGraph();
    }
    if (geneSymbolPage.showHotspotGraph) {
        loadHotspotGraph();
    }
    if (geneSymbolPage.geneCoverageUrl) {
        const gcContainer = $("#qc-gene-coverage-graph-container");
        gcContainer.html('<i class="fa fa-spinner"></i>');
        gcContainer.load(geneSymbolPage.geneCoverageUrl);
    }
    if ($(GENE_VARIANTS_TAG_COUNTS).length) {
        setupTagCountsSummary(GENE_VARIANTS_TAG_COUNTS, geneVariantsTagsChanged);
    }
    $("#gene-variants-show-clinvar").change(reloadGeneVariantsTable);
}

function drawGeneSymbolIdeogram(ideogramData) {
    const container = document.getElementById("gene-symbol-ideogram");
    if (!container) {
        return;
    }
    const ideogram = new Ideogram(container, {
        contig: ideogramData.contig,
        cytobands: ideogramData.cytobands,
        marker: ideogramData.marker,
        onBandClick: (cytoband) => {
            window.location.href = ideogramData.contigUrl + "?" + new URLSearchParams({band: cytoband.name}).toString();
        },
    });
    ideogram.draw();
}

function loadClassificationsHotspotGraph() {
    const url = Urls.classifications_gene_symbol_hotspot_graph(geneSymbolPage.genomeBuildName, geneSymbolPage.geneSymbol);
    $('#classifications-hotspot-graph-container').load(url, function() {
        $(".hotspot-graph", this).attr("hotspot_graph_click_func", "viewGeneClassificationsHotspotGraphClickFunc");
    });
}

function loadHotspotGraph() {
    const url = Urls.gene_symbol_hotspot_graph(geneSymbolPage.genomeBuildName, geneSymbolPage.geneSymbol);
    $('#hotspot-graph-container').load(url, function() {
        $(".hotspot-graph", this).attr("hotspot_graph_click_func", "viewGeneHotspotGraphClickFunc");
    });
}

function viewGeneHotspotGraphClickFunc(transcript_version_id, text, protein_position) {
    $("#gene-variants-grid-filtering-message").empty().html("Filtering to " + text).parent().show();
    geneVariantsGridParams["protein_position"] = protein_position;
    geneVariantsGridParams["protein_position_transcript_version_id"] = transcript_version_id;
    reloadGeneVariantsTable();
}

function viewGeneClassificationsHotspotGraphClickFunc(transcript_version_id, text, protein_position) {
    $("#gene-classifications-grid-filtering-message").empty().html("Filtering to " + text).parent().show();
    geneClassificationsGridParams["protein_position"] = protein_position;
    geneClassificationsGridParams["protein_position_transcript_version_id"] = transcript_version_id;
    classificationGroupingRedraw();
}

// Called by the classification groupings grid before each ajax call
function classificationGroupingFilterExtra(data) {
    data.gene_symbol = geneSymbolPage.geneSymbol;
    Object.assign(data, geneClassificationsGridParams);
}

function clearGeneClassificationsGrid() {
    geneClassificationsGridParams = {};
    classificationGroupingRedraw();
    $("#gene-classifications-grid-filtering-message").empty().parent().hide();
}

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

function geneListsDatatableFilter(data) {
    data.gene_symbol = geneSymbolPage.geneSymbol;
    data.page_id = "gene_symbol";
}
