// jQuery plugin: $span.age(options) to create, $span.age('value', val) to get/set

class AgeWidget {

    constructor(element, options) {
        this.element = $(element);
        this.options = $.extend({value: null, updated: null}, options);

        this.element.addClass('custom-age');

        this.entryText = $('<input>', {class: 'custom-age-number'});

        this.unitSelect = $('<select>', {class: 'custom-age-unit custom-age-select', html: [
                $('<option>', {text: "years", value: ""}),
                $('<option>', {text: "months", value: "months"}),
                $('<option>', {text: "weeks gestation", value: "weeks_gestation"})
        ]});

        this.element.append(this.entryText);
        this.element.append(this.unitSelect);

        this.entryText.on('keyup', () => {this.fieldsUpdated();});
        this.unitSelect.on('change', () => {this.fieldsUpdated();});
        this.unitSelect.chosen({width: '160px'});

        this.value(this.options.value);
    }

    value(value) {
        if (value === undefined) {
            return this.options.value;
        }
        const parts = {num: null, unit: null};

        if (value === null) {
            parts.num = '';
            parts.unit = this.unitSelect.val();
        } else if (typeof value == 'number') {
            parts.num = 'y';
            parts.unit = value;
        } else {
            const unitsM = /^(.*?)(months|weeks_gestation)?$/.exec(`${value}`);
            parts.num = unitsM[1] || '';
            parts.unit = unitsM[2] || '';
        }
        this.options.value = `${parts.num}${parts.unit || 'y'}`;
        this.refreshView(parts);
    }

    fieldsUpdated() {
        const parts = {
            num: this.entryText.val().trim(),
            unit: this.unitSelect.val() || ''
        };
        this.options.value = `${parts.num}${parts.unit}`;
        this.refreshView(parts);

        if (this.options.updated) {
            this.options.updated(null, {value: this.options.value});
        }
    }

    refreshView(parts) {
        this.entryText.val(parts.num);
        this.unitSelect.val( parts.unit );
        this.unitSelect.trigger("chosen:updated");
    }
}

$.fn.age = function(optionsOrMethod, ...args) {
    let methodResult;
    this.each(function() {
        const element = $(this);
        const widget = element.data('ageWidget');
        if (typeof optionsOrMethod === 'string') {
            if (widget) {
                methodResult = widget[optionsOrMethod](...args);
            }
        } else {
            element.data('ageWidget', new AgeWidget(this, optionsOrMethod));
        }
    });
    return methodResult === undefined ? this : methodResult;
};
