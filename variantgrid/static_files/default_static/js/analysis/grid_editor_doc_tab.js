// @ts-check
// analysis/templates/analysis/node_editors/grid_editor_doc_tab.html
/* global retrieveAndUpdateNodeAppearances */ // analysis_nodes.js
function initGridEditorDocTab(nodeId) {
    // Exit Node name edit on enter
    $('#id_name').keypress(function(e) {
      if(e.keyCode == 13) {
        $(this).blur();
      }
    });

    const nodeDocForm = $("form#node-doc-form");
    const options = {
        target: $("#node-doc"),
        success: function () {
            retrieveAndUpdateNodeAppearances([nodeId]);
        },
    };
    nodeDocForm.ajaxForm(options);
    $("#id_name", nodeDocForm).on('input', function() {
        $("#id_auto_node_name", nodeDocForm).prop("checked", false);
    });
    $("#id_auto_node_name", nodeDocForm).change(function() {
        if ($(this).is(":checked")) {
            const autoNodeName = $("#id_auto_name", nodeDocForm).val();
            $("#id_name", nodeDocForm).val(autoNodeName);
        }
    });
}
