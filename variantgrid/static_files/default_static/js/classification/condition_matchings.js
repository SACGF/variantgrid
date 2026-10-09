// @ts-check
// classification/templates/classification/condition_matchings.html
function idRenderer(data) {
    const text = data.normalized_text;
    const id = data.id;

    let textDom;
    if (!text) {
        textDom = $('<span>', {class: 'no-value', text:'<blank>'});
    } else {
        textDom = $('<span>', {text: text});
    }
    // FIXME regenerate URLs
    const aDom = $('<a>', {href:`/classification/condition_matching/${id}`, html: textDom, class: 'hover-link'});

    return aDom.prop('outerHTML');
}
function datatableFilter(data) {
    data.text_filter = $("#text-filter").val();
    data.filter_outstanding = $("#filter-outstanding").is(':checked');
}

function applyFilters() {
    $('#condition-text-datatable').DataTable().ajax.reload();
}

const debouncedApplyFilters = debounce(applyFilters, 250, {maxWait: 1000});
