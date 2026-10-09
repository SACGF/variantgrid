// @ts-check
// analysis/templates/analysis/node_editors/grid_editor_debug_tab.html
/* global hljs */ // js/lib/highlight/highlight.pack.js
function initGridEditorDebugTab(highlightJsUrl, analysisId, nodeId) {
    // This tab is the only user of highlight.js, and it arrives by ajax - fetch the
    // library on first open rather than on every analysis page
    function highlightSql() {
        $('pre code.sql').each(function(i, block) {
            hljs.highlightBlock(block);
        });
    }
    if (window.hljs) {
        highlightSql();
    } else {
        $.ajax({url: highlightJsUrl, dataType: "script", cache: true}).done(highlightSql);
    }

    $('input#show-grid-columns').click(function() {
        const checked = $(this).prop('checked');
        const min_sql = $('#min-sql');
        const grid_sql = $('#grid-sql');
        if (checked) {
            min_sql.hide();
            grid_sql.show();
        } else {
            grid_sql.hide();
            min_sql.show();
        }
    });

    $("button#node-populate-clingen-alleles").click(function() {
        const btn = $(this);
        if (!btn.hasClass("disabled")) {
            $.ajax({
                type: "POST",
                url: Urls.node_populate_clingen_alleles(analysisId, nodeId),
                success: function(data) {
                    btn.addClass("disabled");
                },
            });
        }
    });
}
