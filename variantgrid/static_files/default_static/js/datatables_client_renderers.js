function renderGeneSymbol(data, type, row) {
    let link = $('<a>', {
        href: Urls.view_gene_symbol(data.id),
        class: 'hover-link',
        html: [
            $('<span>', {text: data.id}),
        ]
    });
    return link.prop('outerHTML');
}