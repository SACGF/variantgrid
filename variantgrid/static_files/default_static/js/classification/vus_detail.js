// @ts-check
// classification/templates/classification/vus_detail.html
$('#vus_detail_table').DataTable({
    paging: false,
    searching: true,
    dom: 'fp',
    language: {'searchPlaceholder': 'c.HGVS/gene symbol'},
    columns: [
        {visible: false, searchable:false}, // Allele ID
        {orderable: false, searchable: true}, // Lab (not actually searchable but we put data-search on it)
        {orderable: true, searchable: false}, // Score
        {orderable: true, searchable: false}, // Patients
        {orderable: false, searchable: false}, // PS4
        {orderable: false, searchable: false}, // PP1 BS4
        {orderable: false, searchable: false}, // PM6
        {orderable: false, searchable: false}, // PS2
        {orderable: false, searchable: false}, // PM3
        {orderable: false, searchable: false}, // PP4
        {orderable: false, searchable: false}, // PS3 / BS3
    ],
    rowGroup: {
        dataSrc: 0,
        endRender: null,
        startRender: function ( rows, group ) {
            const alleleId = rows.data()[0][0];
            return $(`#allele-${alleleId}-header`).clone().removeAttr('id');
        }
    },
    order: [[2, 'desc'], [3, 'desc']]
});
