// @ts-check
// classification/templates/classification/classification_grouping_export.html
const groupingExportData = readJsonData("classification-grouping-export-data");

function alleleOriginToggle(filterValue) {
    updateLink();
}

function updateLink() {
    const params = {};

    let labSelectionTab = $('#lab-selection .nav-link.active').attr('id');
    labSelectionTab = labSelectionTab.substring(1, labSelectionTab.length - 4); // drop off leading t and trailing -lab

    $("#custom-lab-warning").hide();
    $("#franklin-allele-origin-warning").hide();

    if (labSelectionTab === 'custom') {
        const labGroupNames = [];
        $('.custom-lab-selection:checked').each((index, value) => {
            labGroupNames.push($(value).attr('id'));
        });
        if (labGroupNames.length == 0) {
            $("#custom-lab-warning").show();
            $("#download-link").attr('href', '#').text('-');
            return;
        }
        labSelectionTab = labGroupNames.join(',');
    }
    params.labs = labSelectionTab;
    const rowsPerFile = $('#rows_per_file').val();
    if (rowsPerFile) {
        params.rows_per_file = `${rowsPerFile}`;
    }
    const since = $('#since').val();
    if (since) {
        params.since = `${since}`;
    }
    const allele_origin = $("input[name='allele-origin-toggle']:checked").val();
    if (allele_origin !== "A") {
        params.allele_origin = allele_origin;
    }

    let formatTab = $('#format .nav-link.active').attr('id');
    formatTab = formatTab.substring(1, formatTab.length - 4);
    params.format = formatTab;
    if (formatTab === 'vcf') {
        params.vcf_target_system = $('input[name=vcf_target_system]:checked').val();
        params.genome_build = $('input[name=vcf_genome_build]:checked').val();
    } else if (formatTab === 'franklin') {
        params.genome_build = $('input[name=franklin_genome_build]:checked').val();
        if (params.allele_origin !== "G") {
            $("#franklin-allele-origin-warning").show();
            $("#download-link").attr('href', '#').text('-');
            return;
        }
    }

    const paramParts = [];
    for (const [key, value] of Object.entries(params)) {
        paramParts.push(`${key}=${encodeURIComponent(value)}`);
    }
    const paramString = '?' + paramParts.join('&');

    const downloadLink = `${groupingExportData.base_url}${ paramString }`;
    $('#download-link').attr('href', downloadLink).text(downloadLink);
}

$(document).ready(() => {
    $('a[data-toggle="tab"]').on('shown.bs.tab', (e) => {
        updateLink();
    });
    $('input').change((e) => {
        updateLink();
    });
    $('input').keyup((e) => {
        updateLink();
    });
    updateLink();
});
