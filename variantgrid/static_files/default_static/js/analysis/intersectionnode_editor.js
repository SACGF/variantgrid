// @ts-check
// analysis/templates/analysis/node_editors/intersectionnode_editor.html
/* global accordionForm, ajaxForm */ // analysis/templates/analysis/node_editors/base_editor.html

function initIntersectionNodeEditor() {
    $(document).ready(function () {
        const intersectionNodeForm = $('form#intersection-node-form');
        accordionForm(intersectionNodeForm);
        ajaxForm(intersectionNodeForm);
    });
}
