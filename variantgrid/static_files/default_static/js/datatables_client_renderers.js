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