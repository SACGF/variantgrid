// @ts-check
// The variant page's overlapping-gene selector (#2524). A variant VEP gave no gene (an SV too long for VEP) has no
// transcript rows, which are what normally drive the gene card, so the select drives it instead: geneFunc is the
// page's gene card loader (variant_details.html), given the selected gene's id and symbol.
// Function declarations only - the analysis editor loads the variant page's scripts into one document repeatedly.

function initGeneOverlapSelect(geneFunc) {
    const select = $("#gene-overlap-select");
    const link = $("#gene-overlap-link");
    select.change(function () {
        const option = select.find("option:selected");
        const geneId = option.data("geneId");
        const geneSymbol = option.val();
        link.attr("href", geneSymbol ? Urls.view_gene_symbol(geneSymbol) : Urls.view_gene(geneId));
        geneFunc(geneId, geneSymbol);
    }).trigger("change");
}
