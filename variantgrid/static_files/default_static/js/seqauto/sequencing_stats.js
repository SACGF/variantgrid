// @ts-check
// seqauto/templates/seqauto/sequencing_stats.html
function showGraph(elementId, title, data) {
    const labels = [];
    const values = [];
    for(let i=0 ; i<data.length ; ++i) {
        labels.push(data[i][0]);
        values.push(data[i][1]);
    }

    const plotData = [{
      labels: labels,
      values: values,
      type: 'pie'
    }];

    const width = 500;
    const height = width;

    const layout = {
        title: {text: title},
        width: width,
        height: height,
    };

    $("#" + elementId).empty();
    Plotly.newPlot(elementId, plotData, layout);
}

$(document).ready(function() {
    const data = readJsonData("sequencing-stats-data");
    const sequencingRunInfo = data.sequencing_run_info;
    const sequencingSampleInfo = data.sequencing_sample_info;

    showGraph('sequencing-run-model', 'Sequencer Model', sequencingRunInfo['sequencer_model']);
    showGraph('sequencing-run-sequencer', 'Sequencer', sequencingRunInfo['sequencer']);
    showGraph('sequencing-run-enrichment_kit', 'Sequencer EnrichmentKit', sequencingRunInfo['enrichment_kit']);

    showGraph('samples-model', 'Sequencer Model', sequencingSampleInfo['sequencer_model']);
    showGraph('samples-sequencer', 'Sequencer', sequencingSampleInfo['sequencer']);
    showGraph('samples-enrichment_kit', 'EnrichmentKit', sequencingSampleInfo['enrichment_kit']);
});
