// @ts-check
// annotation/templates/annotation/variant_annotation_runs.html
const filter_data = {"status": "", "variant_annotation_version_id": ""};
function initVariantAnnotationRuns(genomeBuildSummary) {
    function showSummaryGraph(buildVersionName, buildData) {
        const selector = "build-" + buildVersionName + "-summary-graph";
        const data = [];
        const summary_and_colors = [
            ["Queued", 'rgb(50, 50, 50)'],
            ["Running", 'rgb(50,171, 96)'],
            ["External", 'rgb(70, 130, 180)'],
            ["Error", 'rgb(170, 50, 50)'],
            ["Finished", 'rgb(35, 110, 80)'],
        ];

        for (let i = 0; i < summary_and_colors.length; ++i) {
            const name = summary_and_colors[i][0];
            const color = summary_and_colors[i][1];
            const c = buildData[name] || 0;
            const trace = {
                x: [c],
                y: ["Counts: "],
                orientation: 'h',
                name: name,
                marker: {
                    color: color
                },
                type: 'bar'
            };
            data.push(trace);
        }

        const layout = {
            width: 600,
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
                showticklabels: false
            },
            margin: {
                l: 100,
                r: 100,
                b: 50,
                t: 20,
                //pad: 4
            },
        };
        Plotly.newPlot(selector, data, layout);
    }

    for (const buildName in genomeBuildSummary) {
        const buildData = genomeBuildSummary[buildName];
        for (const vavId in buildData) {
            const vavData = buildData[vavId];  // { pipelineType: {summaryState: count} }
            for (const pipelineType in vavData) {
                const selectorKey = buildName + "-" + vavId + "-" + pipelineType;
                showSummaryGraph(selectorKey, vavData[pipelineType]);
            }
        }
    }

    $('.table-filter').change(() => {
        $('#annotation-runs-table').DataTable().ajax.reload();
    });
}
function idRenderer(data, type, row) {
    return $('<a>', {href:Urls.view_annotation_run(data), text:data}).prop('outerHTML');
}
function datatableFilter(data) {
    data.status = $('#status-filter').val();
    data.variant_annotation_version_id = $('input[name="variant_annotation_version_filter"]:checked').val();
}
