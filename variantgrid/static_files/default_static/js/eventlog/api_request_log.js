// @ts-check
// eventlog/templates/api_request_log.html
function apiRequestLogFilter(data) {
    data.view_name = $("select#api-view-name").val();
    data.errors_only = $("#api-errors-only").is(":checked");
    data.date_from = $("#api-date-from").val();
    data.date_to = $("#api-date-to").val();
}

$(document).ready(() => {
    const reload = () => $('#api-request-log-datatable').DataTable().ajax.reload();
    $("select#api-view-name").change(reload);
    $("#api-errors-only").click(reload);
    $("#api-date-from, #api-date-to").change(reload);
});
