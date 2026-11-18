// @ts-check
// analysis/templates/analysis/node_editors/damagenode_editor.html
/* global setupSlider */ // analysis.js
/* global ajaxForm, nodeEditorHooks */ // analysis/templates/analysis/node_editors/base_editor.html

/**
 * One filter row (analysis/models/nodes/filters/damage_node_filters.py:EffectFilter). A filter is in "all" (All of)
 * when it has no requiredField or that field is ticked, otherwise in "any" (Any of).
 * @typedef {{key: string, label: string, category: string, valueFields: string[], words: string[],
 *            requiredField: ?string, allowNullField: ?string, nullWords: string,
 *            legacyNullField: ?string, legacyNullWords: string,
 *            defaultRequired: ?boolean, defaultAllowNull: ?boolean}} EffectFilter
 */

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

// The readout shows the threshold as the sentence does (setVisibleSliderValue rounds to 2 places, so 0.644 read 0.64)
function showEffectSliderValue(inputSelector, sliderSelector) {
    const value = inputSelector.val();
    $(".slider-value", sliderSelector.parents(".slider-container")).text(value === "" ? "" : formatEffectNumber(value));
    colourSliderValue(inputSelector, sliderSelector);
}

// Set a slider's form field (blank = not set) and move its handle to match
function setEffectSliderValue(inputSelector, sliderSelector, value) {
    const blank = value === null || value === undefined || value === "";
    inputSelector.val(blank ? "" : value);
    sliderSelector.val(blank ? sliderSelector.attr("min") : value);
    showEffectSliderValue(inputSelector, sliderSelector);
}

// Where a slider starts when its filter is added: its published cutoff, else the tool's PP3 cutoff, else its minimum
function effectSliderStart(inputSelector) {
    for (const attr of ["data-default", "data-pathogenic-min", "min"]) {
        const value = inputSelector.attr(attr);
        if (typeof value !== 'undefined' && value !== "") {
            return value;
        }
    }
    return 0;
}

// A number as the server's f"{value:g}" writes it (damage_node_filters.py), so the two sentences match
function formatEffectNumber(value) {
    return String(Number(Number(value).toPrecision(6)));
}

/**
 * The sentence under the groups - damage_node_filters.py:effect_filter_sentence builds the same for the node's
 * summary and hover card, from the same SENTENCE strings.
 * @param {string[]} anyOf
 * @param {string[]} allOf
 * @param {Object<string, string>} sentence
 */
function effectFilterSentence(anyOf, allOf, sentence) {
    const out = $("<span/>");
    if (!anyOf.length && !allOf.length) {
        return out.text(sentence.none);
    }
    const group = function (terms, name) {
        if (terms.length > 1) {
            out.append($("<b/>", {text: sentence[name]}), ": ");
        }
        terms.forEach((term, i) => {
            if (i) {
                out.append($("<span/>", {class: "effect-sentence-separator", text: sentence.separator}));
            }
            out.append($("<span/>", {class: "effect-sentence-term", text: term}));
        });
    };
    out.append(document.createTextNode(sentence.lead));
    if (anyOf.length) {
        group(anyOf, "any_of");
    }
    if (allOf.length) {
        if (anyOf.length) {
            out.append(", ", $("<b/>", {text: sentence.and}), " ");
        }
        group(allOf, "all_of");
    }
    return out.append(".");
}

/**
 * @param {{filters: EffectFilter[], sentence: Object<string, string>, canEdit: boolean}} options
 */
function initDamageNodeEditor(options) {
    $(document).ready(function() {
        const form = $("form#damage-filter");
        ajaxForm(form);

        const editor = $("#effect-filter-editor", form);
        const pool = $(".effect-filter-pool", editor);
        const filters = options.filters;
        /** @type {Object<string, EffectFilter>} */
        const filtersByKey = Object.fromEntries(filters.map(f => [f.key, f]));
        // Rows added this visit that have no value yet (a Variant type row before anything is ticked)
        const pending = new Set();
        let openPalette = null;

        const fieldInputs = name => $("[name='" + name + "']", form);
        const filterRow = f => $(".effect-filter[data-filter='" + f.key + "']", editor);

        /** @type {Object<string, {input: JQuery, slider: JQuery}>} */
        const sliders = {};
        $("input.effect-slider", editor).each(function () {
            const slider = $(this);
            const input = fieldInputs(slider.attr("data-field"));
            setupSlider(input, slider);
            addCutoffBands(input, slider);
            showEffectSliderValue(input, slider);
            // setupSlider's own input handler writes the form field first, so the readout is redone from it here
            slider.on("input", function () {
                showEffectSliderValue(input, slider);
            });
            sliders[slider.attr("data-field")] = {input, slider};
        });
        // A row in a group always has a value - "not set" is removing the row
        $(".effect-filter select option[value='']", editor).prop("hidden", true);

        function fieldHasValue(name) {
            if (sliders[name]) {
                return sliders[name].input.val() !== "";
            }
            const inputs = fieldInputs(name);
            if (inputs.is("select")) {
                return Boolean(inputs.val());
            }
            return inputs.filter(":checked").length > 0;
        }

        function setFieldStart(name) {
            if (sliders[name]) {
                setEffectSliderValue(sliders[name].input, sliders[name].slider, effectSliderStart(sliders[name].input));
                return;
            }
            const inputs = fieldInputs(name);
            if (inputs.is("select")) {
                const start = inputs.attr("data-default");
                const hasStart = typeof start !== 'undefined' && $("option", inputs).filter((i, o) => o.value === start).length;
                inputs.val(hasStart ? start : $("option", inputs).filter((i, o) => o.value !== "").first().val());
            } else if (inputs.length === 1) {
                inputs.prop("checked", true);  // a flag: the row being there is the filter being on
            }  // Variant type starts with nothing ticked
        }

        function clearField(name) {
            if (sliders[name]) {
                setEffectSliderValue(sliders[name].input, sliders[name].slider, "");
                return;
            }
            const inputs = fieldInputs(name);
            if (inputs.is("select")) {
                inputs.val("");
            } else {
                inputs.prop("checked", false);
            }
        }

        // A value as the filter's words show it, or null when it isn't set
        function fieldWords(name) {
            if (!fieldHasValue(name)) {
                return null;
            }
            if (sliders[name]) {
                return formatEffectNumber(sliders[name].input.val());
            }
            const inputs = fieldInputs(name);
            if (inputs.is("select")) {
                return $("option:selected", inputs).text().trim();
            }
            if (inputs.length === 1) {
                return "";
            }
            // Variant type: nothing or everything ticked is no restriction (snpdb/variant_filters.py)
            const ticked = inputs.filter(":checked");
            if (ticked.length === inputs.length) {
                return null;
            }
            return ticked.map((i, input) => $(input).closest("label").text().trim()).get().join(" or ");
        }

        // A row drawing several sliders (Splicing with AlphaGenome) switches each on and off by its own box, and
        // only joins them with "or" after one that is on
        function syncSliderEnable(enable) {
            const name = enable.attr("data-field");
            enable.prop("checked", fieldHasValue(name));
            sliders[name].slider.prop("disabled", !fieldHasValue(name));
        }

        function syncValueConjunctions(row) {
            let earlierOn = false;
            $(".effect-filter-option-group", row).each(function () {
                $(".effect-value-conj", this).toggle(earlierOn);
                const enable = $("input.effect-slider-enable", this);
                earlierOn = earlierOn || (enable.length > 0 && fieldHasValue(enable.attr("data-field")));
            });
        }

        // A node saved under an earlier allow null rule keeps it until switched (splicing, #2836)
        const usesLegacyNull = f => Boolean(f.legacyNullField) && fieldInputs(f.legacyNullField).prop("checked");

        function setLegacyNull(f, on) {
            if (f.legacyNullField) {
                fieldInputs(f.legacyNullField).prop("checked", on);
            }
        }

        const hasValue = f => f.valueFields.some(fieldHasValue);
        const groupOf = f => (!f.requiredField || fieldInputs(f.requiredField).prop("checked")) ? "all" : "any";
        const shown = f => hasValue(f) || pending.has(f.key);

        function filterTerm(f) {
            const parts = [];
            f.valueFields.forEach((name, i) => {
                const value = fieldWords(name);
                if (value !== null) {
                    parts.push(f.words[i].replace("{}", value));
                }
            });
            if (!parts.length) {
                return null;
            }
            let term = parts.join(" or ");
            if (f.allowNullField && groupOf(f) === "all" && fieldInputs(f.allowNullField).prop("checked")) {
                term += " (" + (usesLegacyNull(f) ? f.legacyNullWords : f.nullWords) + ")";
            }
            return term;
        }

        function updateSentence() {
            const terms = {any: [], all: []};
            for (const f of filters) {
                const term = filterTerm(f);
                if (term !== null) {
                    terms[groupOf(f)].push(term);
                }
            }
            $(".effect-sentence", editor).empty().append(effectFilterSentence(terms.any, terms.all, options.sentence));
        }

        function renderAddControls(group, used) {
            const add = $(".effect-group[data-group='" + group + "'] .effect-add", editor).empty();
            $("<button/>", {type: "button", class: "effect-add-button", "data-palette": group,
                            "aria-expanded": String(openPalette === group), text: "+ Add filter"}).appendTo(add);

            const palette = $(".effect-group[data-group='" + group + "'] .effect-palette", editor).empty();
            palette.prop("hidden", openPalette !== group);
            if (openPalette !== group) {
                return;
            }
            $("<div/>", {class: "effect-palette-head"}).append(
                $("<span/>", {text: "Add to " + (group === "any" ? "Any of" : "All of")}),
                $("<button/>", {type: "button", class: "effect-action effect-palette-done", text: "Done"})).appendTo(palette);
            const categories = [...new Set(filters.map(f => f.category))];
            for (const category of categories) {
                const items = $("<span/>", {class: "effect-palette-items"});
                for (const f of filters.filter(x => x.category === category)) {
                    const fixed = group === "any" && !f.requiredField;
                    const title = fixed ? "Always narrows, so it is added under All of" : (used.has(f.key) ? "Already in use" : "");
                    $("<button/>", {type: "button", "data-add": f.key, "data-group": group, text: f.label,
                                    disabled: used.has(f.key) || fixed, title}).appendTo(items);
                }
                $("<div/>", {class: "effect-palette-category"}).append(
                    $("<span/>", {class: "effect-palette-category-name", text: category}), items).appendTo(palette);
            }
        }

        // Place every row from the form's values in one pass, so a saved node and a just-applied preset look the same
        function layout() {
            const groups = {any: [], all: []};
            for (const f of filters) {
                const row = filterRow(f);
                if (shown(f)) {
                    groups[groupOf(f)].push(row);
                } else {
                    row.appendTo(pool);
                }
            }
            const used = new Set(filters.filter(shown).map(f => f.key));
            for (const [group, rows] of Object.entries(groups)) {
                const groupDiv = $(".effect-group[data-group='" + group + "']", editor);
                const rowsDiv = $(".effect-group-rows", groupDiv);
                rows.forEach((row, i) => {
                    row.appendTo(rowsDiv);
                    $(".effect-filter-conj", row).text(i === 0 ? "" : (group === "any" ? "or" : "and"));
                    // Allow null only means something in All of
                    $(".effect-filter-null", row).toggle(group === "all");
                    const legacy = usesLegacyNull(filtersByKey[row.attr("data-filter")]);
                    $(".effect-null-current", row).toggle(!legacy);
                    $(".effect-null-legacy", row).toggle(legacy);
                    $(".effect-filter-legacy-note", row).toggle(legacy && group === "all");
                    $(".effect-move", row)
                        .text(group === "any" ? "↓ All of" : "↑ Any of")
                        .attr("title", group === "any" ? "Every variant must pass this filter"
                                                       : "One of the alternatives a variant can match");
                });
                let empty = "";
                if (!rows.length) {
                    empty = group === "any" ? "Nothing here, so this group lets every variant through."
                                            : "Nothing here, so this group removes nothing.";
                }
                $(".effect-group-empty", groupDiv).text(empty).toggle(!rows.length);
                let hint = "a variant must pass every one, so adding one can only narrow the results";
                if (group === "any") {
                    hint = "a variant needs at least one of these: the first filter narrows the results, "
                         + "each one after that widens them";
                }
                $(".effect-group-hint", groupDiv).text(hint);
                if (options.canEdit) {
                    renderAddControls(group, used);
                }
            }
            $("input.effect-slider-enable", editor).each(function () {
                syncSliderEnable($(this));
            });
            $(".effect-filter", editor).each(function () {
                syncValueConjunctions($(this));
            });
            updateSentence();
        }

        function setRequired(f, required) {
            if (f.requiredField) {
                fieldInputs(f.requiredField).prop("checked", required);
            }
        }

        function addFilter(f, group) {
            setRequired(f, group === "all");
            if (f.allowNullField) {
                fieldInputs(f.allowNullField).prop("checked", Boolean(f.defaultAllowNull));
            }
            setLegacyNull(f, false);
            setFieldStart(f.valueFields[0]);  // a second slider on the row (AlphaGenome) starts off
            pending.add(f.key);
        }

        function removeFilter(f) {
            f.valueFields.forEach(clearField);
            setRequired(f, Boolean(f.defaultRequired));
            if (f.allowNullField) {
                fieldInputs(f.allowNullField).prop("checked", Boolean(f.defaultAllowNull));
            }
            setLegacyNull(f, false);
            pending.delete(f.key);
        }

        editor.on("click", "button", function () {
            const button = $(this);
            const rowKey = button.closest(".effect-filter").attr("data-filter");
            if (button.attr("data-palette")) {
                openPalette = openPalette === button.attr("data-palette") ? null : button.attr("data-palette");
                layout();
                return;
            }
            if (button.hasClass("effect-palette-done")) {
                openPalette = null;
                layout();
                return;
            }
            if (button.attr("data-add")) {
                addFilter(filtersByKey[button.attr("data-add")], button.attr("data-group"));
            } else if (button.hasClass("effect-move")) {
                const f = filtersByKey[rowKey];
                setRequired(f, groupOf(f) === "any");
            } else if (button.hasClass("effect-remove")) {
                removeFilter(filtersByKey[rowKey]);
            } else if (button.hasClass("effect-null-switch")) {
                setLegacyNull(filtersByKey[rowKey], false);
            } else {
                return;
            }
            layout();
            form.trigger("change");  // these set fields without the browser's change event, which node_presets.js watches
        });

        // A second slider on a row (AlphaGenome splicing) switches on and off by itself
        editor.on("change", "input.effect-slider-enable", function () {
            const name = $(this).attr("data-field");
            const row = $(this).closest(".effect-filter");
            const f = filtersByKey[row.attr("data-filter")];
            if (this.checked) {
                setFieldStart(name);
            } else {
                clearField(name);
                if (!hasValue(f)) {
                    removeFilter(f);  // every slider off is the row removed
                    layout();
                    return;
                }
            }
            syncSliderEnable($(this));
            syncValueConjunctions(row);
            updateSentence();
        });
        form.on("input change", updateSentence);

        // ACMG presets (node_presets.js) fill every field: the sliders are ours to move, then the layout follows
        nodeEditorHooks.applyPresetField = function (name, value) {
            if (!sliders[name]) {
                return false;
            }
            setEffectSliderValue(sliders[name].input, sliders[name].slider, value);
            return true;
        };
        nodeEditorHooks.presetApplied = function () {
            pending.clear();
            openPalette = null;
            layout();
        };

        layout();
    });
}
