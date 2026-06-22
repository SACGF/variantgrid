

$.widget('custom.age', {

    options: {
        value: null,
    },

    value: function(value) {
        if (value === undefined) {
            return this.options.value;
        }
        let parts = {num: null, unit: null};

        if (value === null) {
            parts.num = '';
            parts.unit = this.unitSelect.val();
        } else if (typeof value == 'number') {
            parts.num = 'y';
            parts.unit = value;
        } else {
            let unitsM = /^(.*?)(months|weeks_gestation)?$/.exec(`${value}`);
            parts.num = unitsM[1] || '';
            parts.unit = unitsM[2] || '';
        }
        this.options.value = `${parts.num}${parts.unit || 'y'}`;
        this._refreshView(parts);
    },

    _create: function() {
        this.element.addClass('custom-age');

        this.entryText = $('<input>', {class: 'custom-age-number'});

        this.unitSelect = $('<select>', {class: 'custom-age-unit custom-age-select', html: [
                $('<option>', {text: "years", value: ""}),
                $('<option>', {text: "months", value: "months"}),
                $('<option>', {text: "weeks gestation", value: "weeks_gestation"})
        ]});
        $('<option>', {text: "years"}).appendTo();

        this.element.append(this.entryText);
        this.element.append(this.unitSelect);

        $(this.entryText).keyup(() => {this.fieldsUpdated()});
        $(this.unitSelect).change(() => {this.fieldsUpdated()});
        $(this.unitSelect).chosen({width: '160px'});

        this.value(this.options.value);
    },

    fieldsUpdated: function() {
        let parts = {
            num: this.entryText.val().trim(),
            unit: this.unitSelect.val() || ''
        };
        let value = `${parts.num}${parts.unit}`;
        this.options.value = `${parts.num}${parts.unit}`;
        this._refreshView(parts);

        this._trigger("updated", null, {value: this.options.value});
    },

    _refreshView: function(parts) {
        this.entryText.val(parts.num);
        this.unitSelect.val( parts.unit );
        this.unitSelect.trigger("chosen:updated");
    },

    _destroy: function() {
        this.rangeSelect.remove();
        this.element.removeClass('custom-age');
    }
});