// @ts-check
// analysis/templates/analysis/node_editors/builtinfilternode_editor.html
/* global setupSlider */ // analysis.js
/* global ajaxForm */ // analysis/templates/analysis/node_editors/base_editor.html
function setupClinVarStars(bifSelect) {
    // StarsWidget renders the zero option itself - we only show/hide the whole thing
    const starsRow = $("#clinvar-stars-row");

    function checkClinVarStarsVisibility() {
        starsRow.toggle(bifSelect.val() == 'C');
    }

    bifSelect.change(checkClinVarStarsVisibility);
    checkClinVarStarsVisibility();
}

function setupCOSMICCount(bifSelect) {
    function checkCOSMICCount() {
        $("#cosmic-count-widget").toggle(bifSelect.val() === 'M');
    }

    bifSelect.change(checkCOSMICCount);
    checkCOSMICCount();
}

$(document).ready(function() {
    const bifForm = $("form#built-in-filter");
    const bifSelect = $("#id_built_in_filter", bifForm);
    setupClinVarStars(bifSelect);
    setupCOSMICCount(bifSelect);
    setupSlider($("#id_cosmic_count_min"), $("#cosmic-count-slider"));
    ajaxForm(bifForm);
});
