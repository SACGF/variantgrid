// @ts-check
// eventlog/templates/eventlog.html
// Custom renderer specified by datatable_config
function severityRenderer(data, type, row) {
    const domMaker = () => {
        switch (data) {
            case 'E': return $('<div>', {text: 'error', class: 'alert alert-danger text-center', role: 'alert'});
            case 'W': return $('<div>', {text: 'warning', class: 'alert alert-warning text-center', role: 'alert'});
            case 'I': return $('<div>', {text: 'info', class: 'alert alert-info text-center', role: 'alert'});
            case 'D': return $('<div>', {text: 'debug', class: 'alert alert-success text-center', role: 'alert'});
            default: return $('<span>', {text: data});
        }
    };
    const dom = domMaker();
    dom.addClass('severity');
    return dom.prop('outerHTML');
}

function datatableFilter(data) {
    data.filter = $("select#predefined-filters").val();
    data.exclude_admin = $("#exclude-admin").is(":checked");
}

function applyFilters() {
    $('#event-datatable').DataTable().ajax.reload();
}

$(document).ready(() => {
    $("select#predefined-filters").change(applyFilters);
    $("#exclude-admin").click(applyFilters);
});
