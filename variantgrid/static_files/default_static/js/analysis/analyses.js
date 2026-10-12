// @ts-check
// analysis/templates/analysis/analyses.html
/* global setupTagCountsSummary, getTagCountsSummarySelected */ // tag_counts_summary.js
const TAG_COUNTS_SUMMARY = "#analyses-tag-counts";
// Shared filter state - the build toggle and the tag pills both write to it
const gridParams = {};
let tagCountsRequest = null;

function renderAnalysisLink(data, type, row) {
    const dom = $('<span>');
    if (data.locked) {
        $('<i>', {class: 'fa fa-lock fa-lg mr-1', title: 'Locked'}).appendTo(dom);
    }
    $('<a>', {href: data.url, class: 'hover-link', text: data.text}).appendTo(dom);
    return dom.prop('outerHTML');
}

$(document).on('click', '.dt-analysis-settings', function() {
    $("#node-editor-container").load($(this).data('url'));
});

function renderAnalysisTags(data, type, row) {
    return VariantGridFormat.tagsGlobal(data, null, {id: row.id.text});
}

// Called by DataTables before each ajax call - lists go down as JSON, everything else as a plain param
function analysesDatatableFilter(data) {
    for (const [key, value] of Object.entries(gridParams)) {
        data[key] = Array.isArray(value) ? JSON.stringify(value) : value;
    }
}

function reloadAnalysesGrid() {
    $("#analyses-datatable").DataTable().ajax.reload();
}

/* Counting tags over every visible analysis is a group by across the whole tag table, so the
   pills are only fetched once the filter is opened, and recounted when the grid's contents change */
function loadTagCounts() {
    const container = $(TAG_COUNTS_SUMMARY);
    const selected = gridParams["tags"] || [];
    container.html("<i class='fas fa-spinner fa-spin'></i> Counting tags...");
    return $.ajax({
        url: Urls.analysis_list_tag_counts(),
        data: {genome_build_name: gridParams["genome_build_name"] || "",
               analysis_type: gridParams["analysis_type"] || "",
               tag: selected},
        traditional: true,
    }).done(function(html) {
        container.html(html);
        setupTagCountsSummary(container, tagsChanged);
        // A tag with nothing left to match isn't offered, so drop it rather than leave an invisible filter
        const stillSelected = getTagCountsSummarySelected(container);
        if (stillSelected.length !== selected.length) {
            tagsChanged(stillSelected);
        }
    });
}

function ensureTagCounts() {
    if (!tagCountsRequest) {
        tagCountsRequest = loadTagCounts();
    }
    return tagCountsRequest;
}

function refreshTagCounts() {
    if (tagCountsRequest) {
        tagCountsRequest = loadTagCounts();
    }
}

/* Multiple tags mean "carries any of these" */
function tagsChanged(tagIds) {
    if (tagIds.length) {
        gridParams["tags"] = tagIds;
    } else {
        delete gridParams["tags"];
    }
    reloadAnalysesGrid();
}

$(document).ready(() => {
    $('#id_analysis').change(function () {
        const analysisId = $("#id_analysis").val();
        window.location = Urls.analysis(analysisId);
    });

    $("#analyses-grid-filter input[name=genome_build_filter]").change(function() {
        const genomeBuildName = $(this).val();
        if (genomeBuildName) {
            gridParams["genome_build_name"] = genomeBuildName;
        } else {
            delete gridParams["genome_build_name"];
        }
        refreshTagCounts();
        reloadAnalysesGrid();
    });

    $("#analyses-grid-filter input[name=analysis_type_filter]").change(function() {
        const analysisType = $(this).val();
        if (analysisType) {
            gridParams["analysis_type"] = analysisType;
        } else {
            delete gridParams["analysis_type"];
        }
        refreshTagCounts();
        reloadAnalysesGrid();
    });

    $("#analyses-tag-filter").on("show.bs.collapse", ensureTagCounts);

    // Show Group Data changes which analyses the grid holds, so recount once it has reloaded
    $("#user-data-filter-analyses-datatable input[type=checkbox]").change(function() {
        $("#analyses-datatable").one("xhr.dt", refreshTagCounts);
    });

    // Clicking a tag on a row filters to it, the same as picking its pill
    $("#analyses-datatable").on("click", ".grid-tag", function() {
        const tag = $(this).attr("tag_id");
        $("#analyses-tag-filter").collapse("show");
        ensureTagCounts().done(function() {
            const container = $(TAG_COUNTS_SUMMARY);
            $(".summary-count[data-tag='" + tag + "']", container).addClass("selected");
            tagsChanged(getTagCountsSummarySelected(container));
        });
    });
});
