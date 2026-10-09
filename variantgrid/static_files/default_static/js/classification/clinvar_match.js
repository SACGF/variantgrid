// @ts-check
// classification/templates/classification/clinvar_match.html
/* global severityIcon */ // global.js
// Fills in each [data-clinvar-row] placeholder in turn, two requests at a time
function loadClinVarMatches(clinvarKey) {
    loadNext(clinvarKey);
    loadNext(clinvarKey);
}

function loadNext(clinvarKey) {
    const first_data_dom = $('[data-clinvar-row]').first();
    if (first_data_dom.length == 0) {
        return;
    }
    const data = first_data_dom.attr('data-clinvar-row');
    first_data_dom.removeAttr('data-clinvar-row');
    const url = Urls.clinvar_match_detail(clinvarKey);
    first_data_dom.html($('<i class="fa fa-spinner"></i>'));

    $.ajax({
        type: "GET",
        url: url,
        async: true,
        data: {
            data_str: data
        },
        success: (results, textStatus, jqXHR) => {
            first_data_dom.replaceWith(results);
            loadNext(clinvarKey);
        },
        error: (call, status, text) => {
            // sometimes if the client browser goes to sleep we get errors and want to be able to retry
            // if it's going to take a long time
            first_data_dom.replaceWith(
                $('<a>', {
                    href: '#',
                    class: 'ajax-error',
                    html:[severityIcon('C'), "Error Loading Data, Retry?"],
                    'data-clinvar-row': data
                }).click(() => loadNext(clinvarKey)));
        },
        complete: (jqXHR, textStatus) => {}
    });
}
