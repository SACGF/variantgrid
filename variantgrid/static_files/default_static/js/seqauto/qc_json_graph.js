// @ts-check
// seqauto/templates/seqauto/json_graphs/qc_json_graph.html
// Loaded by ajax into a container (named in the data) - sets up each graph whose data hasn't been used yet
function setupQCJsonGraph(jsonDataId) {
    const graphData = readJsonData(jsonDataId);
    const containerSelector = $("#" + graphData.container_name);
    const qc_data = graphData.qc_data;
    const current_label = graphData.current_label;
    const label_column = graphData.label_column;
    const gold_column = graphData.gold_column;

    function createBar(name, color) {
        return {x : [], y: [], type: 'bar', name: name, marker: { color: color} };
    }

    const plotColumn = function() {
        const column = $('option:selected', this).attr('value');
        const labels = qc_data[label_column];
        const is_gold = qc_data[gold_column];
        const column_data = qc_data[column];

        const gold = createBar("gold", "rgb(255,215,0)");
        const active = createBar("current", 'rgba(222,45,38,0.8)');
        const other = createBar("other", 'rgba(204,204,204,0.8)');

        for (let i=0 ; i<labels.length ; ++i) {
            let label = labels[i];
            // JSON can't handle newlines, so had to escape them
            label = label.replace("__newline__", "\n");

            /*  Issue #701 - Sort bar by date.
                We need to put labels/data (NaN is invisible) in all categories so that they are sorted together.
                See https://community.plot.ly/t/how-to-sort-bars-in-a-bar-chart/274/2 */

            gold['x'].push(label);
            active['x'].push(label);
            other['x'].push(label);

            let gold_data = NaN;
            let active_data = NaN;
            let other_data = NaN;

            if (label == current_label) {
                active_data = column_data[i];
            } else if (is_gold[i]) {
                gold_data = column_data[i];
            } else {
                other_data = column_data[i];
            }

            gold['y'].push(gold_data);
            active['y'].push(active_data);
            other['y'].push(other_data);
        }

        const data = [gold, active, other];

        const layout = {
            title: {text: column},
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(0,0,0,0)'
        };

        const dom = $('.plotly-graph-container', containerSelector)[0];
        Plotly.newPlot(dom, data, layout);
    };

    $(document).ready(function() {
        const runStatsColumn = $("select#id_column", containerSelector);
        runStatsColumn.change(plotColumn);
        runStatsColumn.each(plotColumn); // initial plot
    });
}

$('script[id^="qc-json-graph-data-"]').not("[data-graph-setup]").each(function() {
    $(this).attr("data-graph-setup", "true");
    setupQCJsonGraph(this.id);
});
