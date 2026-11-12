// @ts-check
// analysis/templates/analysis/analysis_settings.html
/* global IN_ANALYSIS:writable, ANALYSIS_NODE_ZOOM_PERCENT, setNodeZoomPercent */
function reloadNodes(analysisId, onlyErrors) {
    $.ajax({
        type: "POST",
        url: Urls.analysis_reload(analysisId),
        data: 'only_errors=' + encodeURIComponent(JSON.stringify(!!onlyErrors)),
        success: function(data) {
            if (IN_ANALYSIS) {
                checkAndMarkDirtyNodes();
            }
        }
    });
}

function lockAnalysis(analysisId, lock) {
    if (typeof(lock) === 'undefined') {
        lock = true;
    }
    $.ajax({
        type: "POST",
        url: Urls.analysis_settings_lock(analysisId),
        data: 'lock=' + encodeURIComponent(JSON.stringify(lock)),
        success: function(data) {
            // force a reload of the page - either analysis or analyses listing
            window.location.reload();
        }
    });
}

// The slider changes the graph as it moves and saves to the user's settings when let go.
// Default clears the user's value, so the zoom goes back to the lab / organisation one
function initNodeZoomSlider() {
    const slider = $("#node-zoom-slider");
    const label = $("#node-zoom-value");

    function showZoom(percent) {
        slider.val(percent);
        label.text(percent + "%");
        setNodeZoomPercent(percent);
    }

    function saveZoom(percent) {
        $.ajax({
            type: "POST",
            url: Urls.set_analysis_node_zoom(),
            data: {zoom: percent},
            success: function(data) { showZoom(data.zoom); },
        });
    }

    showZoom(ANALYSIS_NODE_ZOOM_PERCENT);
    slider.on("input", function() { showZoom(Number(slider.val())); });
    slider.on("change", function() { saveZoom(slider.val()); });
    $("#node-zoom-default").click(function() { saveZoom(""); });
}

function initAnalysisSettings(analysisId) {
    IN_ANALYSIS = $("#analysis-and-toolbar-container").length > 0;

    $('button#close-analysis-settings').click(function() {
        $("#analysis-settings-container").parent().empty();
    });

    if (IN_ANALYSIS) {
        initNodeZoomSlider();
    } else {
        $("#force-reload-nodes-button").hide();
        $("#force-reload-error-nodes-button").hide();
        $("#node-zoom-row").hide();
    }

    $('#force-reload-nodes-button').click(function() { reloadNodes(analysisId, false); });
    $('#force-reload-error-nodes-button').click(function() { reloadNodes(analysisId, true); });
    $("button#lock-analysis-button").click(function() { lockAnalysis(analysisId); });
    $("button#unlock-analysis-button").click(function() { lockAnalysis(analysisId, false); });
    // TODO: Lock history...
}
