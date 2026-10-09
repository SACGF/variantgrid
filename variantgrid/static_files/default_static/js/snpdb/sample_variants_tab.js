// @ts-check
// snpdb/templates/snpdb/data/sample_variants_tab.html
/* global setupNodeGrid */ // grid.js
/* global AnalysisDownloadTracker */ // analysis_downloads.js
/* global AnalysisMessagePoller */ // analysis_updates.js
/* global analysisDownloadTracker:writable */ // read via getAnalysisWindow() in grid.js
function gridLoadError() {
    console.log("gridLoadError");
}

function on_error_function() {
    console.log("on_error_function");
}

function showNodeGrid(analysisId, analysisVersion, nodeStatus) {
    const nodeId = nodeStatus.id;
    const version = nodeStatus.version;

    const config_url = Urls.node_grid_config(analysisId, analysisVersion, nodeId, version, "default");
    const handler_url = Urls.node_grid_handler(analysisId);
    const unique_code = "sample-variants-grid";  // ok to be constant on this page
    function noOp() {};
    setupNodeGrid(config_url, handler_url, analysisId, nodeId, version, unique_code,
                  noOp, gridLoadError, on_error_function);
}

function initSampleVariantsTab(sampleId, analysisId, analysisVersion) {
    // Every grid on this page is this sample's, so its pills are read against them
    nodeProbandSampleId = sampleId;
    // The CSV/VCF exports run as Celery jobs then hand back a file - see analysis_downloads.js
    analysisDownloadTracker = new AnalysisDownloadTracker(analysisId, "#sample-variants-downloads");

    const messagePoller = new AnalysisMessagePoller(Urls.nodes_status(analysisId));
    const gridContainer = $("#node-data-grid");

    $("#id_node").change(function() {
        const nodeId = $(this).val();
        gridContainer.empty();
        if (nodeId) {
            $("<table/>").attr({
                id: 'grid-' + nodeId,
                class: 'grid'
            }).appendTo(gridContainer);
            messagePoller.observe_node(nodeId, "ready", (nodeStatus) => showNodeGrid(analysisId, analysisVersion, nodeStatus));
        }
    });
    messagePoller.update_loop();
}
