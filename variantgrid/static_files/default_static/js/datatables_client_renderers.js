function renderGeneSymbol(geneSymbol, type, row) {
    let link = $('<a>', {
        href: Urls.view_gene_symbol(geneSymbol),
        class: 'hover-link',
        html: [
            $('<span>', {text: geneSymbol}),
        ]
    });
    return link.prop('outerHTML');
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
