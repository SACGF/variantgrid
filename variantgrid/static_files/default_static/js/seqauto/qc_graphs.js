// @ts-check
// seqauto/templates/seqauto/qc_graphs.html
const QC_TYPES_TOTALS = readJsonData("qc-graphs-data").qc_type_totals;
console.log("QC");
console.log(QC_TYPES_TOTALS);

function changeGraph() {
    const qc_type = $("#id_qc_type").val();
    const totalField = QC_TYPES_TOTALS[qc_type];
    const percent_selector = $("input#percent");
    const percent_container = percent_selector.closest("div.form-group");
    if (totalField) {
        percent_container.show();
    } else {
        percent_container.hide();
        percent_selector.prop("checked", false);
    }
}

function loadGraph() {
    console.log("loadGraph");
    const graph_selector = $('#qc-column-graph');
    const qc_column = $('#id_qc_column').val();
    const percent = $("input#percent").is(":checked");

    const error_messages = [];
    if (!qc_column) {
        error_messages.push("No QC Column selected");
    }

    const errorContainer = $("#error-container");
    if (error_messages.length) {
        errorContainer.html(error_messages.join(', '));
    } else {
        errorContainer.empty();
        graph_selector.empty();
        graph_selector.addClass('graph-loading');
        const JSON_GRAPH_URL = Urls.qc_column_graph(qc_column, percent);
        graph_selector.load(JSON_GRAPH_URL, function() { graph_selector.removeClass("graph-loading"); });
    }
}

$(document).ready(() => {
    const qc_type = $('#id_qc_type');
    const qc_column = $('#id_qc_column');

    qc_type.change(function() {
        clearAutocompleteChoice(qc_column);
    });
    qc_column.change(changeGraph);
    changeGraph(); // initially hide percent

    const graph_selector = $('#qc-column-graph');
    graph_selector.html("<div id='initial-message'>Select a column to graph</div>");
    $("button#load-graph").click(loadGraph);
});
