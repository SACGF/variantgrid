// @ts-check
// classification/templates/classification/classification_grouping_export.html
let groupingExportBaseUrl = null;  // set by initClassificationGroupingExport

function alleleOriginToggle(filterValue) {
    updateLink();
}

const dateRegex = /^\d{4}-(0?[1-9]|1[0-2])-(0?[1-9]|[12]\d|3[01])$/;
const positiveNum = /^[0-9]+$/;

function verifyDate(dateStr) {
    // returns a boolean to indicate if the value in "since" appears to be a valid
    // number for days ago, or a valid date in the past
    try {
        if (positiveNum.test(dateStr)) {
            if (parseInt(dateStr, 10) >= 0) {
                return true;
            }
        }
    } catch (e) {}
    try {
        if (dateRegex.test(dateStr)) {
            const [year, month, day] = dateStr.split('-');
            // month - 1 adjusts for JavaScript's 0-indexed months
            console.log(`dateStr is a date ${year} ${month} ${day}`);
            let sinceDate = new Date(year, month - 1, day);
            if (sinceDate > new Date()) {
                return false;
            } else {
                return true;
            }
        }
    } catch (e) {}
    return false;
}

function updateLink() {
    const params = {};

    let labSelectionTab = $('#lab-selection .nav-link.active').attr('id');
    labSelectionTab = labSelectionTab.substring(1, labSelectionTab.length - 4); // drop off leading t and trailing -lab

    $("#custom-lab-warning").hide();
    $("#franklin-allele-origin-warning").hide();
    $("#rows-per-file-warning").hide();
    $("#since-warning").hide();

    let error = false;

    if (labSelectionTab === 'custom') {
        const labGroupNames = [];
        $('.custom-lab-selection:checked').each((index, value) => {
            labGroupNames.push($(value).attr('id'));
        });
        if (labGroupNames.length == 0) {
            $("#custom-lab-warning").show();
            error = true;
        }
        labSelectionTab = labGroupNames.join(',');
    }
    params.labs = labSelectionTab;

    const rowsPerFile = $('#rows_per_file').val();
    if (rowsPerFile && rowsPerFile.trim().length > 0) {
        try {
            if (parseInt(rowsPerFile, 10) < 100) {
                $("#rows-per-file-warning").show();
                error = true;
            }
            params.rows_per_file = rowsPerFile;
        } catch (err) {
             $("#rows-per-file-warning").show();
             error = true;
        }
    }

    const since = $('#since').val();
    if (since && since.trim().length > 0) {
        if (!verifyDate(since)) {
            $("#since-warning").show();
            error = true;
        }
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
    const downloadLink = `${groupingExportBaseUrl}${ paramString }`;

    if (error) {
        $("#download-link").attr('href', '#').text('-');
    } else {
        $('#download-link').attr('href', downloadLink).text(downloadLink);
    }
}

function initClassificationGroupingExport(baseUrl) {
    groupingExportBaseUrl = baseUrl;
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
}
