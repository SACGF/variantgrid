// @ts-check
// classification/templates/classification/evidence_keys.html
$(document).ready(() => {
    $('#key-table').dataTable({
        fixedHeader: true,
        dom: 'frBtip',
        paginate: false,
        order: [[0, 'asc']],
        buttons: [
            // Uncomment the below line to put filters back in
            // 'searchPanes'
        ],
        language: {
            searchPanes: {
                collapse: 'Filter'
            }
        },
        columnDefs: [
            {searchPanes: {show: true}, targets: ['category', 'status', 'max-share-level']},
            {searchPanes: {show: false}, targets: '_all'},
        ]
    });
});
