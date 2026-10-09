// @ts-check
// analysis/templates/analysis/analysis_settings.html
/* global IN_ANALYSIS:writable */
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

function initAnalysisSettings(analysisId) {
    IN_ANALYSIS = $("#analysis-and-toolbar-container").length > 0;

    $('button#close-analysis-settings').click(function() {
        $("#analysis-settings-container").parent().empty();
    });

    if (!IN_ANALYSIS) {
        $("#force-reload-nodes-button").hide();
        $("#force-reload-error-nodes-button").hide();
    }

    $('#force-reload-nodes-button').click(function() { reloadNodes(analysisId, false); });
    $('#force-reload-error-nodes-button').click(function() { reloadNodes(analysisId, true); });
    $("button#lock-analysis-button").click(function() { lockAnalysis(analysisId); });
    $("button#unlock-analysis-button").click(function() { lockAnalysis(analysisId, false); });
    // TODO: Lock history...
}
