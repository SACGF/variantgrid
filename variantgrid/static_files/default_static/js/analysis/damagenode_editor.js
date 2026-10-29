// @ts-check
// analysis/templates/analysis/node_editors/damagenode_editor.html
/* global setupSlider, setVisibleSliderValue */ // analysis.js
/* global ajaxForm */ // analysis/templates/analysis/node_editors/base_editor.html

function toggleRequiredFilter(show) {
    if (show) {
        $(".filter-required").show();
    } else {
        $(".filter-required").hide();
        // clear them so they don't apply if invisible
        $("input", ".filter-required").prop("checked", false);
    }
}

function setupVariantClassGroups() {
    // Group checkboxes are unnamed - they drive the individual (named) checkboxes that get submitted
    $(".variant-class-group").each(function() {
        const groupCheckbox = $(this);
        updateVariantClassGroup(groupCheckbox);
        groupCheckbox.change(function() {
            variantClassGroupCheckboxes(groupCheckbox).prop("checked", groupCheckbox.prop("checked"));
        });
        variantClassGroupCheckboxes(groupCheckbox).change(function() {
            updateVariantClassGroup(groupCheckbox);
        });
    });
}

function variantClassGroupCheckboxes(groupCheckbox) {
    return $("#variant-class-individual").find(".variant-class-group-classes[data-group='" +
             groupCheckbox.data("group") + "'] input[type=checkbox]");
}

function updateVariantClassGroup(groupCheckbox) {
    const boxes = variantClassGroupCheckboxes(groupCheckbox);
    const numChecked = boxes.filter(":checked").length;
    groupCheckbox.prop("checked", numChecked > 0 && numChecked === boxes.length);
    groupCheckbox.prop("indeterminate", numChecked > 0 && numChecked < boxes.length);
}

function showIndividualPathogenicPredictions(columnsVersion) {
    // From v4 the rankscores share this class but are split into their own legacy group, so the
    // main toggle reveals only the non-rankscore (raw-score) rows; rankscores use showRankscores().
    let rows = $(".individual-pathogenic-prediction");
    if (columnsVersion >= 4) {
        rows = rows.not(".rankscore-prediction");
    }
    rows.show();
    $("#show-individual-pathogenic-prediction").hide();
}

function showRankscores() {
    $(".rankscore-prediction").show();
    $("#show-rankscore-prediction").hide();
}

// Draw the calibrated PP3/BP4 cutoff bands over a raw-score slider track. The thresholds come
// from data-pathogenic-min / data-benign-max (set in DamageNodeForm from the per-tool TOOLS table);
// sliders without them (rankscores, splice, cosmic) get no bands.
function addCutoffBands(inputSelector, sliderSelector) {
    const sliderMin = Number(inputSelector.attr("min"));
    const sliderMax = Number(inputSelector.attr("max"));
    const range = sliderMax - sliderMin;
    if (!range) {
        return;
    }
    const benignMax = inputSelector.attr("data-benign-max");
    if (typeof benignMax !== 'undefined') {
        const pct = (Number(benignMax) - sliderMin) / range * 100;
        $('<div class="cutoff-band cutoff-band-benign"></div>')
            .css({left: 0, width: pct + "%"})
            .attr("title", "BP4-supporting (benign) ≤ " + benignMax)
            .prependTo(sliderSelector.parent(".slider-track"));
    }
    const pathMin = inputSelector.attr("data-pathogenic-min");
    if (typeof pathMin !== 'undefined') {
        const pct = (Number(pathMin) - sliderMin) / range * 100;
        $('<div class="cutoff-band cutoff-band-pathogenic"></div>')
            .css({left: pct + "%", width: (100 - pct) + "%"})
            .attr("title", "PP3-supporting (pathogenic) ≥ " + pathMin)
            .prependTo(sliderSelector.parent(".slider-track"));
    }
}

// Colour the chosen threshold readout by which calibrated band it falls in. No-op for sliders
// that have no calibrated thresholds, so splice/cosmic/rankscore values keep their default colour.
function colourSliderValue(inputSelector, sliderSelector) {
    const container = sliderSelector.parents(".slider-container");
    const sliderValue = $(".slider-value", container);
    sliderValue.removeClass("in-benign in-intermediate in-pathogenic");
    const benignMax = inputSelector.attr("data-benign-max");
    const pathMin = inputSelector.attr("data-pathogenic-min");
    if (typeof benignMax === 'undefined' && typeof pathMin === 'undefined') {
        return;
    }
    const value = inputSelector.val();
    if (value === '') {
        return;
    }
    const v = Number(value);
    if (typeof pathMin !== 'undefined' && v >= Number(pathMin)) {
        sliderValue.addClass("in-pathogenic");
    } else if (typeof benignMax !== 'undefined' && v <= Number(benignMax)) {
        sliderValue.addClass("in-benign");
    } else {
        sliderValue.addClass("in-intermediate");
    }
}

function setupDamageNodeSlider(inputSelector, sliderSelector, enableInput) {
    setupSlider(inputSelector, sliderSelector);
    addCutoffBands(inputSelector, sliderSelector);
    colourSliderValue(inputSelector, sliderSelector);
    // setupSlider's own input handler updates inputSelector.val() first, so recolour from it here
    sliderSelector.on("input", function () {
        colourSliderValue(inputSelector, sliderSelector);
    });
    const sliderVal = inputSelector.val();

    if (enableInput.length) {
        // A row can hold more than one slider (Splicing + AlphaGenome splicing) sharing its Required /
        // Allow Null, which stay enabled while any of the row's sliders is
        const row = enableInput.parents("tr");
        const frInputs = $("input", $(".filter-required", row));
        const updateRequiredInputs = function () {
            frInputs.prop("disabled", $("input.enable-slider:checked", row).length === 0);
        };

        enableInput.change(function () {
            if ($(this).is(":checked")) {
                sliderSelector.prop("disabled", false);
                // A slider with a suggested start (data-default) begins there rather than at its minimum
                const defaultValue = inputSelector.attr("data-default");
                if (inputSelector.val() === '' && typeof defaultValue !== 'undefined') {
                    inputSelector.val(defaultValue);
                    sliderSelector.val(defaultValue);
                    setVisibleSliderValue(inputSelector, sliderSelector, defaultValue);
                }
                colourSliderValue(inputSelector, sliderSelector);
            } else {
                setVisibleSliderValue(inputSelector, sliderSelector, "");
                sliderSelector.val(sliderSelector.attr("min"));
                sliderSelector.prop("disabled", true);
                inputSelector.val("");
                colourSliderValue(inputSelector, sliderSelector);
            }
            updateRequiredInputs();
        });
        enableInput.prop("checked", sliderVal !== '');
        updateRequiredInputs();
    }

    sliderSelector.prop("disabled", sliderVal === '');
}

/**
 * @param {{columnsVersion: number, showRequiredFilter: boolean, hasIndividualPathogenicPredictions: boolean,
 *          hasRankscorePredictions: boolean, sliders: {input: string, slider: string, enable: string}[]}} options
 */
function initDamageNodeEditor(options) {
    let showRequiredFilter = options.showRequiredFilter;

    $(document).ready(function() {
        ajaxForm($("form#damage-filter"));
        setupVariantClassGroups();

        for (const s of options.sliders) {
            setupDamageNodeSlider($(s.input), $(s.slider), $(s.enable));
        }

        $("#toggle-required-filter").click(function (event) {
            event.preventDefault();
            showRequiredFilter = !showRequiredFilter;
            toggleRequiredFilter(showRequiredFilter);
        });
        $("#toggle-variant-classes").click(function (event) {
            event.preventDefault();
            $("#variant-class-individual").toggle();
        });
        $("#show-individual-pathogenic-prediction a").click(function (event) {
            event.preventDefault();
            showIndividualPathogenicPredictions(options.columnsVersion);
        });
        $("#show-rankscore-prediction a").click(function (event) {
            event.preventDefault();
            showRankscores();
        });

        if (showRequiredFilter) {
            $(".filter-required").show();
        }
        if (options.hasIndividualPathogenicPredictions) {
            showIndividualPathogenicPredictions(options.columnsVersion);
        }
        if (options.hasRankscorePredictions) {
            showRankscores();
        }
    });
}
