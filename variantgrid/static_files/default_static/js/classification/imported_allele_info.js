// @ts-check
// classification/templates/classification/imported_allele_info.html
function datatableFilter(data) {
    $('#imported-allele-info-datatables-filter .form-check-input').each((index, element) => {
        const $element = $(element);
        data[$element.attr('id')] = $element.is(':checked');
    });
}
function applyFilters() {
    $('#imported-allele-info-datatables').DataTable().ajax.reload();
}
function render_validation(data, type, row) {
    const dom = $('<div>');
    const mainIncludeLine = $('<div>');
    if (data.include) {
        $('<i>', {class: 'fas fa-check-circle text-success'}).appendTo(dom);
    } else {
        $('<i>', {class: 'fas fa-times-circle text-danger'}).appendTo(dom);
    }
    for (const tag of data.tags) {
        const color = tag.severity == "E" ? "text-danger" : "text-warning";
        $('<div>', {text: tag.label, style: 'font-weight: bold', class: 'mt-1 ' + color}).appendTo(dom);
    }
    if (data.message) {
        $('<div>', {html: limitLengthSpan(data.message, 100), class: 'text-secondary'}).appendTo(dom);
    }
    return dom.prop('outerHTML');
}
