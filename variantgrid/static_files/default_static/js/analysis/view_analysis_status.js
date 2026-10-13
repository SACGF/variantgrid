// @ts-check
// analysis/templates/analysis/view_analysis_status.html
/* global Plotly, setupNodeTypeSelect, blankToNull */

const SCHEDULER_SERIES = [
    ["Waiting", 'rgb(150, 150, 150)'],
    ["Backing off", 'rgb(240, 173, 78)'],
    ["Dispatched", 'rgb(70, 130, 180)'],
    ["Running", 'rgb(50, 171, 96)'],
    ["Lapsed", 'rgb(170, 50, 50)'],
];

const ALL_NODES_SERIES = [
    ["Queued", 'rgb(50, 50, 50)'],
    ["Running", 'rgb(50, 171, 96)'],
    ["User Error", 'rgb(255, 79, 0)'],
    ["Error", 'rgb(170, 50, 50)'],
    ["Finished", 'rgb(35, 110, 80)'],
];

/**
 * One stacked horizontal bar, a trace per series (the Pipeline Runs page style)
 * @param {string} selector element id
 * @param {Object<string, number>} summaryData count per series name
 * @param {Array<Array<string>>} seriesAndColors [name, color] in stacking order
 */
function showSummaryGraph(selector, summaryData, seriesAndColors) {
    const data = seriesAndColors.map(([name, color]) => ({
        x: [summaryData[name] || 0],
        y: ["counts"],
        orientation: 'h',
        name: name,
        marker: {color: color},
        type: 'bar',
    }));

    const layout = {
        width: 800,
        height: 100,
        barmode: 'stack',
        showlegend: true,
        legend: {orientation: 'h'},
        xaxis: {
            autorange: true,
            showgrid: false,
            zeroline: false,
            showline: false,
            ticks: '',
            showticklabels: false,
        },
        margin: {l: 100, r: 100, b: 50, t: 20},
    };
    Plotly.newPlot(selector, data, layout);
}

/**
 * @param {Object} options
 * @param {Object<string, number>|null} options.schedulerSummary null when the scheduler is idle (no graph drawn)
 * @param {Object<string, number>} options.nodesStatusSummary
 * @param {Object} options.nodeTypesDisplay
 */
function initViewAnalysisStatus(options) {
    if (options.schedulerSummary) {
        showSummaryGraph("analysis-scheduler-graph", options.schedulerSummary, SCHEDULER_SERIES);
    }
    showSummaryGraph("analysis-node-summary-graph", options.nodesStatusSummary, ALL_NODES_SERIES);

    setupNodeTypeSelect("#id_node_class", options.nodeTypesDisplay);

    const filterForm = $("#analysis-node-issues-filter");
    const reloadGrid = () => $('#analysis-node-issues-datatable').DataTable().ajax.reload();
    filterForm.find("select, input").on("change", reloadGrid);
    $("#clear-filters").on("click", () => {
        filterForm.find("select, input:not([type=hidden])").val("");
        filterForm.find("select").trigger("change.nodeTypeSelect");
        history.replaceState(null, "", window.location.pathname);  // drop the 'Reload matching' redirect's filters
        reloadGrid();
    });
    filterForm.on("submit", (event) => {
        const numNodes = $('#analysis-node-issues-datatable').DataTable().page.info().recordsDisplay;
        if (!confirm(`Reload ${numNodes} failed nodes (and the other failed nodes in their analyses)?`)) {
            event.preventDefault();
        }
    });
}
