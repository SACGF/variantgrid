// @ts-check
// uicore/templates/uicore/widgets/optional_range_widget.html
/* Ticking the box shows the slider and has it submitted. Unticked, the slider is disabled so the field posts
   nothing and saves as null */
registerDomProcessor('[data-optional-range]', (container) => {
    const setBox = container.find(".optional-range-set");
    const range = container.find(".optional-range-input");
    const valueLabel = container.find(".optional-range-value");
    const suffix = container.attr("data-value-suffix") || "";

    function update() {
        const isSet = setBox.is(":checked");
        range.prop("disabled", !isSet || setBox.prop("disabled")).prop("hidden", !isSet);
        valueLabel.text(isSet ? range.val() + suffix : "");
    }

    setBox.on("change", update);
    range.on("input", update);
    update();
});
