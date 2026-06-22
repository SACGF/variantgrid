function renderGeneSymbol(geneSymbol, type, row) {
    let link = "";
    if (geneSymbol) {
        let linkObj = $('<a>', {
            href: Urls.view_gene_symbol(geneSymbol),
            class: 'hover-link',
            html: [
                $('<span>', {text: geneSymbol}),
            ]
        });
        link = linkObj.prop('outerHTML');
    }
    return link;
}

function renderAnalysisAuditLogSummary(summary, type, row) {
    let summaryText = summary['summary_text'];
    if (summaryText) {
        return "<span>" + summaryText + "</span>";
    }
    let changes = summary['changes'];
    if (changes) {
        const hideValues = new Set(['valid', 'status', 'version', 'shadow_color', 'appearance_version']);
        let changesSummary = "<table class='table'>"
        changesSummary += "<tr><th>field</th><th>old</th><th>new</th></tr>"
        for (const [key, value] of Object.entries(changes)) {
            if (!hideValues.has(key)) {
                changesSummary += `<tr><td>${key}</td><td>${value[0]}</td><td>${value[1]}</td>`;
            }
        }
        changesSummary += "</table>"
        return changesSummary
    }

}


function renderExpandAnalysisAuditLogEntry(x) {
    return formatJson({
        "changes": x["changes"],
        "additional_data": x["additional_data"],
    })
}
