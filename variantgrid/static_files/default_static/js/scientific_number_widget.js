function isBlankOrNull(value) {
    if (value == null || (typeof(value) == 'string' && value.trim().length == 0)) {
        return true;
    }
    return false;
}

function toFixedString(value) {
    let strValue = value;
    if (isBlankOrNull(value)) {
        return null;
    }
    value = Number(value);

    if (!isNaN(value)) {
        if (Math.floor(value) === value) {
            // note large number will still be turned scientific but we shouldn't see those
            strValue = value.toFixed(0);
        } else {
            strValue = value.toFixed(10);
            // strip trailing 0s, if Math.foor value != value then we should be guarenteed to have a decimal
            // but just double check
            if (strValue.indexOf('.') !== -1) {
                while (strValue[strValue.length - 1] == '0') {
                    strValue = strValue.slice(0, -1);
                }
                // trailing number was over factionDigits
                // e.g. 1.00000000000000000302 which would then set to
                // 1.0000000000000000, after stripping 0s it was 1. so now strip the "."
                if (strValue[strValue.length - 1] == '.') {
                    strValue = strValue.slice(0, -1);
                }
            }
        }
    }
    return strValue;
}
function toPercent(value, multiplier) {
    const numValue = Number(value);
    if (!isBlankOrNull(value) && !isNaN(numValue)) {
        const percent = numValue * (multiplier || 100);
        if (percent >= 10) {
            return percent.toFixed(1) + '%';
        } else if (percent >= 1) {
            return percent.toFixed(2) + '%';
        } else if (percent === 0) {
            return '0%';
        } else if (percent < 0.000001) {
            return '<0.000001%';
        } else {
            let precise = percent.toPrecision(2);
            if (precise.indexOf('.') !== -1 && precise.endsWith('0')) {
                precise = precise.substring(0, precise.length-1);
            }
            return precise + '%';
        }
    } else {
        return null;
    }
}

// jQuery plugin: $input.scientific(options) wraps an input with a live % representation

class ScientificWidget {

    constructor(element, options) {
        this.element = $(element);
        this.options = $.extend({placeholder: null, placeholder_short: null, tooltip: null, multiplier: null}, options);

        this.wrapper = $('<span>', {class: 'custom-scientific'});
        this.wrapper.insertAfter(this.element);
        this.element.appendTo(this.wrapper);

        this.notation = $('<span>', {class: 'notation', text: '', title: this.options.tooltip});
        this.notation.appendTo(this.wrapper);

        this.element.attr('placeholder', this.options.placeholder_short);
        this.element.attr('title', this.options.placeholder);
        this.element.on('keyup', () => {this.refreshScientificNote();});
        this.value(this.element.val());
        this.element.on('blur', () => {
            const oldValue = this.element.val();
            this.value(this.element.val());
            const neatValue = this.element.val();
            if (oldValue != neatValue) {
                this.element.trigger('change');
            }
        });

        this.element.attr('customPopulate', true);
        this.element.on('onpopulate', (event, val) => {
            this.value(val);
        });
    }

    value(value) {
        this.element.val(toFixedString(value));
        this.refreshScientificNote();
    }

    refreshScientificNote() {
        let scientificValue = toPercent(this.element.val(), this.options.multiplier);
        if (scientificValue === null) {
            scientificValue = '';
        }
        this.notation.text(`${scientificValue}`);
    }
}

$.fn.scientific = function(optionsOrMethod, ...args) {
    let methodResult;
    this.each(function() {
        const element = $(this);
        const widget = element.data('scientificWidget');
        if (typeof optionsOrMethod === 'string') {
            if (widget) {
                methodResult = widget[optionsOrMethod](...args);
            }
        } else {
            element.data('scientificWidget', new ScientificWidget(this, optionsOrMethod));
        }
    });
    return methodResult === undefined ? this : methodResult;
};
