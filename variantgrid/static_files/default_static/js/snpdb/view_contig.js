// @ts-check
// snpdb/templates/snpdb/genomics/view_contig.html

let contigGenesContigId = null;

function initViewContig(contigId) {
    contigGenesContigId = contigId;
}

// data-datatable-data callback for the genes grid
function contigGenesFilter(data) {
    data.contig_id = contigGenesContigId;
}
