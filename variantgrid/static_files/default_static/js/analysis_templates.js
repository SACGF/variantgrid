/* This is include with AnalysisNode form media */
function analysisVariableNodeFieldSetup() {
    let nodeId = $(this).parents("#node-editor-wrapper").attr("node_id");
    let field = $(this).attr("field")
    $(".add-analysis-variable-button", this).click(function() {
        analysisVariable(nodeId, field, 'add', function() {
            addAnalysisVariableButton(nodeId, field);
            lockNodeField(nodeId, field, true);
        });
    });

    let fieldSet = analysisNodeVariables[nodeId];
    if (fieldSet) {
        if (fieldSet.has(field)) {
            lockNodeField(nodeId, field, true);
        }
    }
}

$(document).ready(function() {
    $(".analysis-variable-node-field-wrapper").each(analysisVariableNodeFieldSetup);

});