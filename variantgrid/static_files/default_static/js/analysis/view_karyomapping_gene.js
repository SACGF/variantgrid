// @ts-check
// analysis/templates/analysis/karyomapping/view_karyomapping_gene.html
function reverseArray(array) {
    return array.slice().reverse();
}

function showKaryotypeScatter() {
    const pageData = readJsonData("view-karyomapping-gene-data");
    const KARYOTYPE_BIN_SCATTER_DATA = pageData.karyotype_bin_scatter_data;
    const KARYOTYPE_BIN_LABELS = pageData.karyotype_bin_labels;
    const data = [];

    for(let i=0 ; i<KARYOTYPE_BIN_LABELS.length ; i++) {
        const k = KARYOTYPE_BIN_LABELS[i];
        const k_data = KARYOTYPE_BIN_SCATTER_DATA[k];
        // console.log("k: " + k);
        // console.log(k_data);
        let x_data = k_data['x'];
        if (!x_data.length) {
            x_data = [null]; // always show even if empty
        }

        const trace = {
          x: x_data,
          y: Array(x_data.length).fill(k),
          mode: 'markers',
          type: 'scatter',
          name: k,
          text: k_data['text'],
          marker: { size: 12 },
          visible: true,
        };
        data.push(trace);
    }

    const description = 'Karyomapping ' + escapeHtml(pageData.proband_sample);
    const geneDescription = escapeHtml(pageData.gene_symbol) + " ('" + escapeHtml(pageData.strand) + "' strand) Up: " + pageData.upstream_kb + "KB, Down: " + pageData.downstream_kb + ")KB";
    const coordinates = escapeHtml(pageData.iv);

    const layout = {
        title: {text: [description, geneDescription, coordinates].join('\n')},
        showlegend: false,
        type: 'category',
        xaxis : {
            // showgrid: false,
            showline: false,
        },
        yaxis : {
            showticklabels: true,
            showline: false,
            categoryorder: "array",
            categoryarray:  reverseArray(KARYOTYPE_BIN_LABELS),
        },

    };

    Plotly.newPlot('karyotype-graph', data, layout);
}

$(document).ready(function() {
    showKaryotypeScatter();
});
