/* This is include with AnalysisNode form media */
function analysisVariableNodeFieldSetup() {
    const nodeId = $(this).parents("#node-editor-wrapper").attr("node_id");
    const field = $(this).attr("field");
    $(".add-analysis-variable-button", this).click(function() {
        analysisVariable(nodeId, field, 'add', function() {
            addAnalysisVariableButton(nodeId, field);
            lockNodeField(nodeId, field, true);
        });
    });

    const fieldSet = analysisNodeVariables[nodeId];
    if (fieldSet) {
        if (fieldSet.has(field)) {
            lockNodeField(nodeId, field, true);
        }
    }
}

$(document).ready(function() {
    $(".analysis-variable-node-field-wrapper").each(analysisVariableNodeFieldSetup);

});