// @ts-check
// classification/templates/classification/tags/classification_groupings.html
/* global blankToNull */ // global.js
/* global classificationGroupingFilterExtra */ // optionally defined by the including page
/* global debounceClassificationGroupingRedraw:writable, classificationGroupingClinicalSignificance:writable */ // re-executed when the tag is in an AJAX-loaded fragment
function alleleOriginToggle() {
    classificationGroupingRedraw();
}

function classificationGroupingRedraw() {
    $('#vc-datatable').DataTable().ajax.reload();
    classificationGroupingSummaryLoad();
}
debounceClassificationGroupingRedraw = debounce(classificationGroupingRedraw);

classificationGroupingClinicalSignificance = null;

function classificationGroupingClinicalSignificanceSelect(clinicalSignificance) {
    classificationGroupingClinicalSignificance = clinicalSignificance;
    classificationGroupingRedraw();
}

function classificationGroupingSummaryLoad() {
    const summary = $('#classification-grouping-summary');
    if (!summary.length) {
        return;
    }
    const data = {};
    classificationGroupingFilter(data);
    $.getJSON(Urls.classification_grouping_counts(), data, (response) => {
        summary.empty().toggleClass('filtered', !!classificationGroupingClinicalSignificance);
        for (const entry of response.counts) {
            const clinicalSignificance = entry.clinical_significance || 'none';
            const selected = classificationGroupingClinicalSignificance === clinicalSignificance;
            $('<a>', {
                    class: 'summary-count' + (selected ? ' selected' : ''),
                    href: '#',
                    title: selected ? 'Show all classifications' : 'Only show ' + entry.label
                })
                .append($('<span>', {class: 'c-pill cs ' + entry.css_class, text: entry.label}))
                .append($('<span>', {class: 'count', text: entry.count.toLocaleString()}))
                .on('click', (event) => {
                    event.preventDefault();
                    classificationGroupingClinicalSignificanceSelect(selected ? null : clinicalSignificance);
                })
                .appendTo(summary);
        }
    });
}

$(document).ready(() => {
    $('.filter').on("change", function() {
        classificationGroupingRedraw();
    });
    $('#vc-datatable').on('init.dt', classificationGroupingSummaryLoad);
});

function classificationGroupingFilter(data) {
    if (classificationGroupingClinicalSignificance) {
        data.clinical_significance = classificationGroupingClinicalSignificance;
    }
    const allele_origin_filter_value = blankToNull($("input[name='allele-origin-toggle']:checked").val());
    if (allele_origin_filter_value) {
        data.allele_origin = allele_origin_filter_value;
        const germline_only = allele_origin_filter_value === "G";
        const column = $('#vc-datatable').DataTable().column(3);
        column.visible(!germline_only);
    }

    // allows implementation of classificationFilterExtra if surrounded code provides more filters
    if (typeof(classificationGroupingFilterExtra) === 'function') {
        classificationGroupingFilterExtra(data);
    }
}
function downloadAs(mode) {
    const data = {};
    classificationGroupingFilter(data);
    data["genome_build"] = $('#vc-datatable').attr('data-genome-build');
    const querystring = EncodeQueryData(data, true);
    let url = null;
    if (mode === 'redcap') {
        url = Urls.export_classifications_grid_redcap() + "?" + querystring;
    } else {
        url = Urls.export_classifications_grid() + "?" + querystring;
    }
    window.location = url;
    return false;
}
