// @ts-check
// analysis/templates/analysis/node_editors/_node_presets.html
/* global nodeEditorHooks */ // analysis/templates/analysis/node_editors/base_editor.html

// Fill one form field from a preset value. An editor with a control the generic rules can't drive (the Effect
// node's sliders) registers nodeEditorHooks.applyPresetField(name, value) and returns true for the fields it owns.
function setPresetField(form, name, value) {
    if (nodeEditorHooks.applyPresetField && nodeEditorHooks.applyPresetField(name, value)) {
        return;
    }
    const inputs = $("[name='" + name + "']", form);
    if (!inputs.length) {
        return; // a form field this annotation version's editor doesn't draw
    }
    const type = inputs.attr("type");
    if (type === "checkbox") {
        if (Array.isArray(value)) {
            const values = value.map(String);
            inputs.each(function () {
                this.checked = values.includes(this.value);
            });
        } else {
            inputs.prop("checked", Boolean(value));
        }
        inputs.trigger("change");
    } else if (type === "radio") {
        inputs.each(function () {
            this.checked = this.value === String(value);
            // Bootstrap 4 data-toggle="buttons" shows the checked radio through its label's active class
            $(this).closest("label.btn").toggleClass("active", this.checked);
        });
        inputs.filter(":checked").trigger("change");
    } else {
        inputs.val(value === null ? "" : value).trigger("change");
    }
}

function applyNodePreset(form, values) {
    for (const [name, value] of Object.entries(values)) {
        setPresetField(form, name, value);
    }
    if (nodeEditorHooks.presetApplied) {
        nodeEditorHooks.presetApplied(values);
    }
}

/**
 * @param {{select: string, form: string,
 *          presets: {key: string, label: string, description: string, values: Object<string, any>}[]}} options
 */
function initNodePresets(options) {
    const select = $(options.select);
    const description = select.siblings(".node-preset-description");
    select.on("change", function () {
        const preset = options.presets.find(p => p.key === this.value);
        if (!preset) {
            return;
        }
        applyNodePreset($(options.form), preset.values);
        description.text(preset.description);
        this.value = ""; // the dropdown is an action - what the form now holds is what the node will filter on
    });
}
