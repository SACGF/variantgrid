// @ts-check
// analysis/templates/analysis/karyomapping/view_karyomapping_gene.html
function reverseArray(array) {
    return array.slice().reverse();
}

/* options: probandSample, geneSymbol, strand, upstreamKb, downstreamKb, iv - for the title */
function showKaryotypeScatter(scatterData, binLabels, options) {
    const data = [];

    for(let i=0 ; i<binLabels.length ; i++) {
        const k = binLabels[i];
        const k_data = scatterData[k];
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

    const description = 'Karyomapping ' + escapeHtml(options.probandSample);
    const geneDescription = escapeHtml(options.geneSymbol) + " ('" + escapeHtml(options.strand) + "' strand) Up: " + options.upstreamKb + "KB, Down: " + options.downstreamKb + "KB";
    const coordinates = escapeHtml(options.iv);

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
            categoryarray:  reverseArray(binLabels),
        },

    };

    Plotly.newPlot('karyotype-graph', data, layout);
}
