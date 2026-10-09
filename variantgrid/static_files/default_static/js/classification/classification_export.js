// @ts-check
// classification/templates/classification/classification_export.html
const exportData = readJsonData("classification-export-data");
let baseUriApi = exportData.base_url;
let baseUriRedirect = exportData.base_url_redirect;
let paramString = '';
const errorString = null;
let lastFormat = null;
if (baseUriApi.indexOf('localhost') !== -1) {
    baseUriApi = Urls.classification_export_api();
    baseUriRedirect = Urls.classification_export_redirect();
}

function scrollToAlleleOrigin() {
    document.getElementById("allele-origin-toggle").scrollIntoView({ behavior: "smooth", block: "end"});
    $("#allele-origin-toggle").delay(400).fadeOut().fadeIn();
}

function generateUrl() {
    const params = {};
    let inputError = null;
    const format = $("input[name='format']:checked").val();
    if (format !== lastFormat) {
        $('.custom-option').hide();
        $(`[data-format=${format}]`).show('fast');
        lastFormat = format;
    }
    params.share_level = $("input[name='share_level']:checked").val();

    params.build = $("input[name='genome_build']:checked").val();

    const allele_origin = $("input[name='allele-origin-toggle']:checked").val();
    if (allele_origin !== "A") {
        params.allele_origin = allele_origin;
    }

    params.type = format;

    const rowsPerFile = $("#rows_per_file").val();
    if (rowsPerFile && rowsPerFile.length) {
        params.rows_per_file = rowsPerFile;
    }

    const rowLimit = $("#row_limit").val();
    if (rowLimit && rowLimit.length) {
        params.row_limit = rowLimit;
    }

    if (format === 'json') {
        const fullDetail = $("input[name='json_full_detail']:checked");
        if (fullDetail && fullDetail.length) {
            params.full_detail = 'true';
        }

    } else if (format == 'franklin') {
        if (!allele_origin || allele_origin !== "G") {
            inputError = 'Franklin download currently only supports Germline.<br/>Please set <a class="hover-link" onclick="scrollToAlleleOrigin()">Allele Origin</a> to <b>Germline</b> if you wish to export in Franklin format.';
        }
    } else if (format === 'vcf') {
        params.target_system = $("input[name='vcf_target_system']:checked").val();

    } else if (format === 'csv') {

        const valueFormat = $("input[name='value_format']:checked").val();
        if (valueFormat) {
            params.value_format = valueFormat;
        }

        const htmlHandling = $("input[name='html_handling']:checked").val();
        if (htmlHandling) {
            params.html_handling = htmlHandling;
        }

        const excludeTransient = $("input[name='exclude_transient']:checked");
        if (excludeTransient && excludeTransient.length) {
            params.exclude_transient = 'true';
        }

        const fullDetail = $("input[name='csv_full_detail']:checked");
        if (fullDetail && fullDetail.length) {
            params.full_detail = 'true';
        }
    }

    const labSelectionMode = $("input[name='lab_selection_mode']:checked").val();
    const selectedLabs = $("input[name='lab_selection']:checked");
    if (selectedLabs.length) {
        const labArray = selectedLabs.toArray().map(el => $(el).val());

        params[labSelectionMode === 'exclude' ? 'exclude_labs' : 'include_labs'] = labArray.join(',');
    } else if (labSelectionMode === 'include') {
        inputError = 'You must include at least 1 lab';
    }
    if (format === 'lab_compare'){
        if (selectedLabs.length !== 2 || labSelectionMode !== "include") {
            inputError = 'Lab Compare requires the <b>inclusion</b> of 2 labs.';
        }
    }

    const selectedOrgs = $("input[name='org_selection']:checked");
    if (selectedOrgs.length) {
        const orgArray = selectedOrgs.toArray().map(el => $(el).val());
        params.exclude_orgs = orgArray.join(',');
    }

    const since = $("input[name='since']").val().trim();
    if (since.length) {
        params.since = since;
    }
    const benchmark = $("input[name='benchmark']:checked");
    if (benchmark.length) {
        params.benchmark = 'true';
    }

    const allele = $("input[name='allele']").val();
    if (allele && allele.trim().length) {
        params.allele = allele.trim();
    }

    const paramParts = [];
    for (const [key, value] of Object.entries(params)) {
        paramParts.push(`${key}=${encodeURIComponent(value)}`);
    }
    paramString = '?' + paramParts.join('&');
    let url;
    const downloadMode = $("input[name='download_mode']:checked").val();
    let target = null;
    if (downloadMode === 'api') {
        url = baseUriApi;
    } else {
        url = baseUriRedirect;
        target = '_blank';
    }
    url += paramString;
    const a = $('#export-url');
    const error_box = $('#export-error');
    const error_text = $('#export-error-text');
    if (inputError) {
        error_text.html(inputError);
        a.hide();
        error_box.show();
    } else {
        error_box.hide();
        a.attr('href', url);
        a.attr('target', target);
        a.text(url);
        if (downloadMode === 'api') {
            a.addClass('download-link');
        } else {
            a.removeClass('download-link');
        }
        a.show();
    }
}

function alleleOriginToggle(filterValue) {
    this.generateUrl();
}

$(document).ready(() => {
    $('#export-fields input').change(() => {
       this.generateUrl();
    }).keyup(() => {
       this.generateUrl();
    });
    this.generateUrl();
});
