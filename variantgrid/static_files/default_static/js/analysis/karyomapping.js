// @ts-check
// analysis/templates/analysis/karyomapping/trio_karyomapping.html, analysis/templates/analysis/tags/trio_karyomapping_qc.html
function reverseArray(array) {
    return array.slice().reverse();
}

/* options: probandSample, geneSymbol, strand, upstreamKb, downstreamKb, iv - for the title,
            geneStart, geneEnd - shaded behind the points, variantIdLookup - point text to variant id for click through */
function showKaryotypeScatter(scatterData, binLabels, options) {
    const data = [];

    for(let i=0 ; i<binLabels.length ; i++) {
        const k = binLabels[i];
        const k_data = scatterData[k] || {x: [], text: []};
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
        title: {text: [description, geneDescription, coordinates].join('<br>')},
        showlegend: false,
        xaxis : {
            showline: false,
        },
        yaxis : {
            showticklabels: true,
            showline: false,
            categoryorder: "array",
            categoryarray:  reverseArray(binLabels),
        },
        shapes: [{
            type: 'rect', layer: 'below', xref: 'x', yref: 'paper',
            x0: options.geneStart, x1: options.geneEnd, y0: 0, y1: 1,
            fillcolor: '#6c757d', opacity: 0.15, line: {width: 0},
        }],
        annotations: [{
            xref: 'x', yref: 'paper', x: (options.geneStart + options.geneEnd) / 2, y: 1,
            yanchor: 'top', showarrow: false, text: escapeHtml(options.geneSymbol),
        }],
    };

    Plotly.newPlot('karyotype-graph', data, layout);
    const graph = /** @type {any} */ (document.getElementById('karyotype-graph'));
    graph.on('plotly_click', function (event) {
        const variantId = options.variantIdLookup[event.points[0].text];
        if (variantId) {
            window.open(Urls.view_variant(variantId), '_blank');
        }
    });
}

/* contigQC: list of {contig, father_percent, mother_percent, father_in_phase, father_out_of_phase,
                      mother_in_phase, mother_out_of_phase} in contig order */
function showKaryomappingQCGraph(selector, contigQC) {
    const contigs = contigQC.map(c => c.contig);
    const parents = [
        {name: 'Father', key: 'father', color: '#1f77b4'},
        {name: 'Mother', key: 'mother', color: '#d62728'},
    ];
    const data = parents.map(parent => ({
        x: contigs,
        y: contigQC.map(c => c[parent.key + '_percent']),
        customdata: contigQC.map(c => [c[parent.key + '_in_phase'], c[parent.key + '_out_of_phase']]),
        hovertemplate: parent.name + ' %{x}: %{y:.1f}% in phase (%{customdata[0]} in / %{customdata[1]} out)<extra></extra>',
        name: parent.name,
        type: 'scatter',
        mode: 'lines+markers',
        line: {color: parent.color},
    }));

    const layout = {
        title: {text: 'Proband in phase with each parent, per chromosome'},
        xaxis: {type: 'category', title: {text: 'Chromosome'}},
        yaxis: {range: [0, 100], title: {text: '% of parent het sites in phase'}},
        shapes: [{
            type: 'line', xref: 'paper', x0: 0, x1: 1, y0: 50, y1: 50,
            line: {dash: 'dash', color: '#6c757d', width: 1},
        }],
        legend: {orientation: 'h'},
        margin: {t: 40},
    };
    Plotly.newPlot(selector, data, layout, {responsive: true});
}
