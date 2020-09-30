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
    let numValue = Number(value);
    if (!isBlankOrNull(value) && !isNaN(numValue)) {
        let percent = numValue * (multiplier || 100);
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

$.widget('custom.scientific', {

    options: {
        placeholder: null,
        tooltip: null,
        multiplier: null
    },

    _create: function() {
        this.wrapper = $('<span>', {class: 'custom-scientific'});
        this.wrapper.insertAfter(this.element);
        this.element.appendTo(this.wrapper);

        this.notation = $('<span>', {class: 'notation', text: '', title: this.options.tooltip});
        //this.element.insertAfter(this.notation);
        this.notation.appendTo(this.wrapper);
        // save this so we can unbind it later
        this.updateBinding = () => {
            this.refreshScientificNote();
        };

        this.element.attr('placeholder', this.options.placeholder_short);
        this.element.attr('title', this.options.placeholder);
        this.element.keyup(this.updateBinding);
        this.value(this.element.val());
        this.element.bind('blur', () => {
            let oldValue = this.element.val();
            this.value(this.element.val());
            let neatValue = this.element.val();
            if (oldValue != neatValue) {
                this.element.trigger('change');
            }
        });

        this.element.attr('customPopulate', true);
        this.element.bind('onpopulate', (event, val) => {
            this.value(val);
        });
    },

    value: function(value) {
        if (value === undefined) {
            return this.value;
        }
        let strValue = toFixedString(value);
        this.element.val(strValue);
        this.refreshScientificNote();
    },

    fieldsUpdated: function() {
        this.refreshScientificNote();
        //this._trigger("updated", null, {value: this.element.val()});
    },

    refreshScientificNote: function() {
        let scientificValue = toPercent(this.element.val(), this.options.multiplier);
        if (scientificValue === null) {
            scientificValue = ''; // this.options.placeholder;
        }
        this.notation.text(`${scientificValue}`);
    },

    _destroy: function() {
        this.element.unbind('keyup', this.updateBinding);
        this.element.insertBefore(this.wrapper);
        this.wrapper.remove();
    }
});