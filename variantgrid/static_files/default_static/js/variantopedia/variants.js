// @ts-check
// variantopedia/templates/variantopedia/variants.html
// What the grid is showing - the controls only become this when Search is clicked. Set by initAllVariants
let gridExtraFilters = null;

function getAlignmentFiles() {
    return [];
}

// Called by DataTables before each ajax call (and so also by the CSV download, which is
// built from the table's live ajax params)
function allVariantsDatatableFilter(data) {
    if (gridExtraFilters) {
        data["extra_filters"] = JSON.stringify(gridExtraFilters);
    }
}

function filterGrid() {
    const table = $("#all-variants-datatable");
    if ($.fn.DataTable.isDataTable(table)) {
        table.DataTable().ajax.reload();
    }
}

function checkedValues(selector) {
    return $(selector).filter(":checked").map(function() { return $(this).val(); }).get();
}

function setChecked(input, checked) {
    // Bootstrap 4 button groups track state via the 'active' class on the wrapping label
    input.prop("checked", checked);
    input.closest("label").toggleClass("active", checked);
}

function readFilters() {
    const contigIds = checkedValues(".contig-filter").map(Number);
    return {
        contig_ids: contigIds,
        non_standard_contigs: $("#non-standard-contigs").is(":checked"),
        gene_symbols: $("#id_gene_symbols").val() || [],
        variant_types: checkedValues(".variant-type-filter"),
        min_count: Number($("#min_count").val()) || 0,
    };
}

function describeFilters(filters, numVariantTypes) {
    const parts = [];
    const contigNames = filters.contig_ids.map(function(contigId) {
        return $(".contig-filter[value='" + contigId + "']").data("contig-name");
    });
    if (filters.non_standard_contigs) {
        contigNames.push("non-standard contigs");
    }
    if (contigNames.length) {
        parts.push("chromosomes: " + contigNames.join(", "));
    }
    if (filters.gene_symbols.length) {
        parts.push("genes: " + filters.gene_symbols.join(", "));
    }
    if (filters.variant_types.length && filters.variant_types.length < numVariantTypes) {
        parts.push("types: " + filters.variant_types.join(", "));
    }
    if (filters.min_count) {
        parts.push("seen in at least " + filters.min_count + " sample(s)");
    }
    return parts.join("; ");
}

function saveFilters(genomeBuildName, filters) {
    $.ajax({
        type: "POST",
        url: Urls.set_all_variants_filter(genomeBuildName),
        contentType: "application/json",
        data: JSON.stringify(filters),
    });
}

function isSelective(filters) {
    // Variant type and min count alone still walk the whole variant table
    return Boolean(filters.contig_ids.length || filters.non_standard_contigs || filters.gene_symbols.length);
}

function showFilters(filters, numVariantTypes) {
    const description = describeFilters(filters, numVariantTypes);
    $("#all-variants-filter").text(description);
    $("#all-variants-filter-description").toggle(Boolean(description));
    $("#all-variants-unselective").toggle(!isSelective(filters));
}

function applyFallbacks(defaultContigId) {
    // Nothing selective would show no rows, and no variant type would show every type without
    // saying so - tick what's going to be searched, so the buttons describe the grid
    const filters = readFilters();
    if (!isSelective(filters) && defaultContigId !== null) {
        setChecked($(".contig-filter[value='" + defaultContigId + "']"), true);
    }
    if (!filters.variant_types.length) {
        setChecked($(".variant-type-filter"), true);
    }
    return readFilters();
}

function filtersMatchGrid() {
    const filters = readFilters();
    const asString = (f) => JSON.stringify([
        (f.contig_ids || []).slice().sort(), Boolean(f.non_standard_contigs),
        (f.gene_symbols || []).slice().sort(), (f.variant_types || []).slice().sort(),
        Number(f.min_count) || 0,
    ]);
    return asString(filters) === asString(gridExtraFilters);
}

function showPendingSearch() {
    // A solid button while the controls differ from what the grid is showing
    const pending = !filtersMatchGrid();
    $("#search-button").toggleClass("btn-primary", pending)
                       .toggleClass("btn-outline-primary", !pending);
}

function runSearch(options) {
    const filters = applyFallbacks(options.defaultContigId);
    gridExtraFilters = filters;
    showFilters(filters, options.numVariantTypes);
    showPendingSearch();
    saveFilters(options.genomeBuildName, filters);
    filterGrid();
}

function tickGeneSymbolContigs(genomeBuildName, geneSymbol) {
    // Tick the gene's chromosomes, otherwise a gene outside the current selection shows nothing
    const url = Urls.api_gene_symbol_detail(geneSymbol) + "?genome_build=" + genomeBuildName;
    return $.getJSON(url, function(data) {
        for (const gene of data.genes || []) {
            for (const version of gene.versions || []) {
                for (const contig of version.contigs || []) {
                    setChecked($(".contig-filter[value='" + contig.id + "']"), true);
                }
            }
        }
    });
}

function geneSymbolsChanged(genomeBuildName) {
    const geneSymbols = $("#id_gene_symbols").val() || [];
    $.when.apply($, geneSymbols.map((geneSymbol) => tickGeneSymbolContigs(genomeBuildName, geneSymbol))).always(showPendingSearch);
}

function resetFilters(genomeBuildName) {
    // An empty saved filter set means "use the defaults"
    $.ajax({
        type: "POST",
        url: Urls.set_all_variants_filter(genomeBuildName),
        contentType: "application/json",
        data: JSON.stringify({}),
        success: function() { window.location.reload(); },
    });
}

function setCheckedAll(selector, checked) {
    setChecked($(selector), checked);
    showPendingSearch();
}

/* options: initialFilters, genomeBuildName, defaultContigId, numVariantTypes */
function initAllVariants(options) {
    gridExtraFilters = options.initialFilters;
    $(document).ready(() => setupAllVariantsControls(options));
}

function setupAllVariantsControls(options) {
    $(".contig-filter, #non-standard-contigs, .variant-type-filter, #min_count").change(showPendingSearch);
    $("#min_count").keydown(function(event) {
        if (event.which === 13) {  // Enter searches, rather than doing nothing
            event.preventDefault();
            $("#search-button").click();
        }
    });
    $("#id_gene_symbols").change(() => geneSymbolsChanged(options.genomeBuildName));
    $("#contig-all").click(() => setCheckedAll(".contig-filter", true));
    $("#contig-none").click(() => setCheckedAll(".contig-filter", false));
    $("#variant-type-all").click(() => setCheckedAll(".variant-type-filter", true));
    $("#variant-type-none").click(() => setCheckedAll(".variant-type-filter", false));
    $("#search-button").click(() => runSearch(options));
    $("#all-variants-reset").click(() => resetFilters(options.genomeBuildName));
    showFilters(readFilters(), options.numVariantTypes);
    showPendingSearch();
}
