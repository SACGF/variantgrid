// @ts-check
// analysis/templates/analysis/node_editors/conservationnode_editor.html
/* global setupSlider, setVisibleSliderValue */ // analysis.js
/* global ajaxForm */ // analysis/templates/analysis/node_editors/base_editor.html

function showIndividualConservationScores() {
    $(".individual-conservation-scores").show();
    $("#show-individual-conservation-scores").hide();
}

function setSliderEnabled(sliderSelector, enabled) {
    sliderSelector.prop("disabled", !enabled);
}

function setupConservationNodeSlider(inputSelector, sliderSelector, enableInput) {
    setupSlider(inputSelector, sliderSelector);
    const setEnabledIfChecked = function () {
        setSliderEnabled(sliderSelector, enableInput.is(":checked"));
    };
    setEnabledIfChecked();
    enableInput.change(setEnabledIfChecked);
}

/**
 * @param {{individualFields: string[]}} options - the conservation scores the annotation version populates
 */
function initConservationNodeEditor(options) {
    function setIndividualSlidersFromMaster(masterValue) {
        if ($("#id_use_individual_sliders").is(":checked")) {
            return;
        }
        const fraction = parseFloat(masterValue);
        for (const fieldName of options.individualFields) {
            const inputSelector = $("#id_" + fieldName);
            const sliderSelector = $("#" + fieldName + "-slider");
            const min = parseFloat(inputSelector.attr("min"));
            const max = parseFloat(inputSelector.attr("max"));
            const newSliderValue = min + (max - min) * fraction;
            sliderSelector.val(newSliderValue);
            inputSelector.val(newSliderValue);
            setVisibleSliderValue(inputSelector, sliderSelector, newSliderValue);
        }
    }

    $(document).ready(function() {
        ajaxForm($("form#conservation-filter"));

        const useIndividualSliders = $("#id_use_individual_sliders");
        const masterInput = $("#id_any_scaled_min");
        const masterSlider = $("#any_scaled_min-slider");
        setupSlider(masterInput, masterSlider);

        const setEnabledIfNotChecked = function () {
            setSliderEnabled(masterSlider, !useIndividualSliders.is(":checked"));
        };
        useIndividualSliders.change(setEnabledIfNotChecked);
        setEnabledIfNotChecked();
        masterSlider.on('change', function () {
            setIndividualSlidersFromMaster(this.value);
        });

        for (const fieldName of options.individualFields) {
            setupConservationNodeSlider($("#id_" + fieldName), $("#" + fieldName + "-slider"), useIndividualSliders);
        }
        $("#show-individual-conservation-scores a").click(function (event) {
            event.preventDefault();
            showIndividualConservationScores();
        });

        // after all setup, now set from initial master value if applicable
        setIndividualSlidersFromMaster(masterInput.val());
        if (useIndividualSliders.is(":checked")) {
            showIndividualConservationScores();
        }
    });
}
