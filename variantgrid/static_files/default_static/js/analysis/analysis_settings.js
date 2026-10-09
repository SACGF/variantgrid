// @ts-check
// analysis/templates/analysis/analysis_settings.html
/* global IN_ANALYSIS:writable, ANALYSIS_SETTINGS:writable */
IN_ANALYSIS = $("#analysis-and-toolbar-container").length > 0;

function reloadNodes(onlyErrors) {
    const data = readJsonData("analysis-settings-data");
    $.ajax({
        type: "POST",
        url: Urls.analysis_reload(data.analysis_id),
        data: 'only_errors=' + encodeURIComponent(JSON.stringify(!!onlyErrors)),
        success: function(data) {
            if (IN_ANALYSIS) {
                checkAndMarkDirtyNodes();
            }
        }
    });
}

function lockAnalysis(lock) {
    if (typeof(lock) === 'undefined') {
        lock = true;
    }
    const data = readJsonData("analysis-settings-data");
    $.ajax({
        type: "POST",
        url: Urls.analysis_settings_lock(data.analysis_id),
        data: 'lock=' + encodeURIComponent(JSON.stringify(lock)),
        success: function(data) {
            // force a reload of the page - either analysis or analyses listing
            window.location.reload();
        }
    });
}

$(document).ready(function() {
    ANALYSIS_SETTINGS = readJsonData("analysis-settings-data").new_analysis_settings;

    $('button#close-analysis-settings').click(function() {
        $("#analysis-settings-container").parent().empty();
    });

    if (!IN_ANALYSIS) {
        $("#force-reload-nodes-button").hide();
        $("#force-reload-error-nodes-button").hide();
    }

    $('#force-reload-nodes-button').click(function() { reloadNodes(false); });
    $('#force-reload-error-nodes-button').click(function() { reloadNodes(true); });
    $("button#lock-analysis-button").click(function() { lockAnalysis(); });
    $("button#unlock-analysis-button").click(function() { lockAnalysis(false); });
    // TODO: Lock history...

});
