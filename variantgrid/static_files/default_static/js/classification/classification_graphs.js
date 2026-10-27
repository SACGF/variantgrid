// @ts-check
// classification/templates/classification/classification_graphs.html
/* global Plotly */ // plotly-custom
/* global defaultLayout */ // plotly_helper.js
/* global Urls */ // django_js_reverse

const CLASSIFICATION_GRAPH_COLORS = ['#666666', // Other
                                     '#7777ff', // Benign
                                     '#aaaaff', // Likely benign
                                     '#cccccc', // VUS
                                     '#ffaaaa', // Likely pathogenic
                                     '#ff7777', // Pathogenic
                                    ];
const CLASSIFICATION_GRAPH_LABELS = ['Other', 'Benign', 'Likely Benign', 'VUS', 'Likely Pathogenic', 'Pathogenic'];
const UNCLASSIFIED_COLOR = '#1f77b4';


/**
 * Each graph loads separately so the page and the quicker graphs don't wait on the slowest one
 * @param {string} selector
 * @param {string} graph - key of CLASSIFICATION_GRAPHS in classification/views/views.py
 * @param {function(any): void} plot
 */
function loadClassificationGraph(selector, graph, plot) {
    const container = $(selector);
    container.html('<div class="text-center text-muted my-4">...loading...</div>');
    $.getJSON(Urls.classification_graph_data(graph))
        .done((data) => {
            container.empty();
            plot(data);
        })
        .fail(() => {
            container.html('<div class="text-danger my-4">Unable to load graph</div>');
        });
}


/**
 * @param {{showUnclassified: boolean}} options
 */
function initClassificationGraphs(options) {
    const colors = [...CLASSIFICATION_GRAPH_COLORS];
    const labels = [...CLASSIFICATION_GRAPH_LABELS];
    if (options.showUnclassified) {
        colors.push(UNCLASSIFIED_COLOR);
        labels.push('Unclassified');
    }

    loadClassificationGraph('#classifications-graph', 'allele_counts', (counts) => {
        const data = [{
            values: counts.values,
            labels: counts.labels,
            type: 'pie',
            sort: false,
            marker: {colors: [...colors].reverse()},
        }];
        Plotly.newPlot('classifications-graph', data, defaultLayout('Unique Classified Alleles'));
    });

    loadClassificationGraph('#gene-classifications-graph', 'gene_counts', (geneData) => {
        geneData.forEach((series, i) => {
            series.marker = {color: colors[i]};
        });
        const layout = defaultLayout('Unique Classified Alleles by Gene');
        layout.barmode = 'stack';
        Plotly.newPlot('gene-classifications-graph', geneData, layout);
    });

    loadClassificationGraph('#criteria-usage-container', 'criteria_counts', (acmgBySignificance) => {
        const criteriaContainer = $('#criteria-usage-container');
        const width = $('#container_table_content').width() * 0.7;
        for (let clinSig = labels.length - 1; clinSig >= 0; clinSig--) {
            const data = acmgBySignificance[`${clinSig}`];
            if (data) {
                const graphId = 'criteria-usage-' + clinSig;
                criteriaContainer.append($('<div />').attr({id: graphId}));
                const layout = defaultLayout(labels[clinSig], width, 400);
                layout.barmode = 'stack';
                Plotly.newPlot(graphId, data, layout);
            }
        }
    });
}
