// @ts-check
// analysis/templates/analysis/tags/base_related_analyses.html
// This is called from analysis_output_node_downloads tag. The export runs under Celery, so poll
// in place - this page has no analysis toolbar, but its icons live as long as the page does.
function downloadNodeCSV(link, analysisId, nodeId, nodeVersion) {
    const anchor = $(link);
    const state = anchor.data("downloadState");
    if (state === "pending") {
        return;  // Already being generated
    }
    if (state === "done") {
        window.location.href = anchor.data("downloadUrl");  // Grab it again without regenerating
        return;
    }
    const icon = $("i.csv-icon", anchor);
    const gridParam = {
        node_id: nodeId,
        version_id: nodeVersion,
        rows: 0,
        export_type: 'csv',
        use_canonical_transcripts: true, // if available
    };
    const url = Urls.node_grid_export(analysisId) + "?" + EncodeQueryData(gridParam);

    function setIcon(state, iconClass, title) {
        anchor.data("downloadState", state);
        icon.attr("class", iconClass);
        anchor.attr("title", title);
    }

    setIcon("pending", "icon fas fa-spinner fa-spin", "Preparing download...");
    $.getJSON(url, function (data) {
        poll_cached_generated_file(Urls.cached_generated_file_check(data.cgf_id),
            function (d) {
                anchor.data("downloadUrl", d.url);
                setIcon("done", "icon fas fa-download", "Download ready");
                window.location.href = d.url;
            },
            function (d) {
                setIcon(null, "icon csv-icon", "Download failed: " + ((d && d.exception) || ""));
            },
            function (progress) {
                anchor.attr("title", `Preparing download - ${Math.floor(100 * (progress || 0))}%`);
            });
    }).fail(function () {
        setIcon(null, "icon csv-icon", "Could not start the download");
    });
}
