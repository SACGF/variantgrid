// @ts-check
// classification/templates/classification/clinvar_key_summary.html
/* global severityIcon */ // global.js
const clinvarKeySummaryData = readJsonData("clinvar-key-summary-data");
function records(data) {
    data.clinvar_key = clinvarKeySummaryData.clinvar_key;
}
function batches(data) {
    data.clinvar_key = clinvarKeySummaryData.clinvar_key;
}
function render_allele_origin_bucket(data) {
    return VCTable.allele_origin_bucket_label(data, null, "horizontal");
}
function renderId(data, type, row) {
    const allele_origin_bucket = data.allele_origin_bucket;
    const bucket_dom = VCTable.allele_origin_bucket_label(allele_origin_bucket);
    return $('<div>', {
        style: "display:flex; position: relative; top: -12px; left: 6px",
        html:[
        bucket_dom,
        "CE_" + data.id
    ]}).prop('outerHTML');
}
function batchId(data, type, row) {
    return $('<div>', {
        style: "display:flex; position: relative; top: -12px; left: 12px",
        html:[
        "CB_" + data
    ]}).prop('outerHTML');
}
function renderStatus(data, type, row) {
    let dom = null;
    if (data === 'D') {
        dom = severityIcon('success');
        dom.attr('title', 'Up to date');
    } else if (data === 'E') {
        dom = severityIcon('error');
        dom.attr('title', 'ClinVar conversion issues, unable to submit');
    } else if (data === 'N') {
        dom = $('<i class="fas fa-cloud-upload-alt text-success"></i>');
        dom.attr('title', 'New submission');
    } else if (data == 'C') {
        dom = $('<i class="fas fa-cloud-upload-alt text-primary"></i>');
        dom.attr('title', 'Changes to existing ClinVar record pending');
    } else if (data == 'X') {
        dom = $('<span/>', {html: [$('<i class="fas fa-globe-americas">'), $('</i><i class="fas fa-times"></i>')]});
        dom.attr('title', 'Classification has been specifically excluded from ClinVar Export');
    } else {
        dom = $('<span>', {text: data});
    }
    if (dom.attr('title')) {
        dom.attr('data-toggle', 'tooltip');
    }
    return dom.prop('outerHTML');
}
function renderReleaseStatus(data, type, row) {
    let dom = null;
    if (data === 'R') {
        dom = severityIcon('success');
        dom.attr('title', 'Release when ready');
        dom.attr('data-toggle', 'tooltip');
    } else if (data === 'H') {
        dom = $('<i class="fas fa-clock"></i>');
        dom.attr('title', 'On hold');
        dom.attr('data-toggle', 'tooltip');
    } else {
        dom = $('<span>', {text: data});
    }
    if (row.status === "E") {
        dom.css('opacity', 0.4);
    }
    return dom.prop('outerHTML');
}
function renderSCV(data, type, row) {
    if (!data) {
        return $('<span>', {text: '-', class:'no-value'}).prop('outerHTML');
    } else {
        return data;
    }
}
/*
function renderId(data, type, row) {
    let content = [$('<span/>', {text: data.genome_build}), $("<br/>"), $('<span/>', {text: data.c_hgvs})];
    // let id = data.id;
    // let elem = $('<a/>', {href: `/classification/clinvar_export/${id}`, html: content, class:'id-link'});
    let elem = $('<span>', {html: content});
    return elem.prop('outerHTML');
}
*/
function renderBatches(data, type, row) {
    if (!data) {
        return $('<span>', {text: '-', class:'no-value'}).prop('outerHTML');
    } else {
        return $('<span>', {text: data.join(", "), class:'text-monospace'}).prop('outerHTML');
    }
}
