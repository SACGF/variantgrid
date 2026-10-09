// @ts-check
// classification/templates/classification/clinvar_match.html
/* global severityIcon */ // global.js
const clinvarMatchData = readJsonData("clinvar-match-data");
function noOp() {

}
function loadNext() {
    const first_data_dom = $('[data-clinvar-row]').first();
    if (first_data_dom.length == 0) {
        return;
    }
    const data = first_data_dom.attr('data-clinvar-row');
    first_data_dom.removeAttr('data-clinvar-row');
    const url = Urls.clinvar_match_detail(clinvarMatchData.clinvar_key);
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
            loadNext();
        },
        error: (call, status, text) => {
            // sometimes if the client browser goes to sleep we get errors and want to be able to retry
            // if it's going to take a long time
            first_data_dom.replaceWith(
                $('<a>', {
                    onClick: "loadNext()",
                    href: '#',
                    class: 'ajax-error',
                    html:[severityIcon('C'), "Error Loading Data, Retry?"],
                    'data-clinvar-row': data
            }));
        },
        complete: (jqXHR, textStatus) => {}
    });
}
loadNext();
loadNext();
