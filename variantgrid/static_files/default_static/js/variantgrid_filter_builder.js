/* Column filter builder for the variant grids.

   Emits the rule JSON the grid engine already speaks:

     {"groupOp": "AND", "rules": [{"field": "locus__position", "op": "lt", "data": "2500"}]}

   which is what JqGrid.get_q turns into a Q object, what a grid request carries as
   '_search=true&filters=<json>', and what FilterNode persists as FilterNodeItem rows. Fields,
   their types and their choices all come from the definition JSON's filterBuilder block
   (@see datatable_filter_fields_from_colmodels), so the list is whatever the grid can filter on.

   Two callers:
     - DataTableDefinition puts one behind the "Filter grid..." toolbar button on any adapter served
       grid, and sends the rules with the table's ajax params.
     - The FilterNode editor mounts one directly and POSTs the rules to the node.
*/

const VariantGridFilterBuilder = (function() {
    "use strict";

    const NUMERIC_TYPES = ["int", "float"];

    /* params:
         container    where to render (jQuery)
         filterBuilder  the definition's filterBuilder block: {fields: [...], operations: [...]}
         onApply      called with the rule JSON when the user presses Find
         onReset      called when the user clears the filters (defaults to onApply with no rules)
         readOnly     render the rules without the controls that change them
         applyLabel   text for the apply button
    */
    const VariantGridFilterBuilder = function(params) {
        this.container = params.container;
        this.fields = (params.filterBuilder || {}).fields || [];
        this.operations = (params.filterBuilder || {}).operations || [];
        this.onApply = params.onApply;
        this.onReset = params.onReset;
        this.onChange = params.onChange;
        this.readOnly = Boolean(params.readOnly);
        this.applyLabel = params.applyLabel || "Find";
        this.groupOp = "AND";
        this.rulesContainer = null;
        this.groupOpSelect = null;
    };

    VariantGridFilterBuilder.prototype = {

        fieldByName: function(fieldName) {
            for (const field of this.fields) {
                if (field.field === fieldName) {
                    return field;
                }
            }
            return null;
        },

        /* Operations that make sense for the field's type - ordering/contains operations on a
           select field would compare against the stored code rather than the label */
        operationsForField: function(field) {
            const type = (field && field.type) || "text";
            if (type === "select") {
                return this.operations.filter(o => ["eq", "ne", "in", "ni", "nu", "nn"].indexOf(o.op) !== -1);
            }
            if (NUMERIC_TYPES.indexOf(type) !== -1 || type === "date") {
                return this.operations.filter(o => ["cn", "nc", "bw", "bn", "ew", "en"].indexOf(o.op) === -1);
            }
            return this.operations;
        },

        buildFieldSelect: function(selectedField) {
            const select = $("<select>", {class: "form-control form-control-sm filter-rule-field"});
            if (selectedField && !this.fieldByName(selectedField)) {
                // A saved rule on a column this grid no longer offers - keep it selectable so
                // loading and re-saving a FilterNode doesn't quietly drop the rule
                $("<option>", {value: selectedField, text: selectedField}).appendTo(select);
            }
            for (const field of this.fields) {
                $("<option>", {value: field.field, text: field.label}).appendTo(select);
            }
            if (selectedField) {
                select.val(selectedField);
            }
            return select;
        },

        buildOperationSelect: function(field, selectedOp) {
            const select = $("<select>", {class: "form-control form-control-sm filter-rule-op"});
            for (const operation of this.operationsForField(field)) {
                $("<option>", {value: operation.op, text: operation.label}).appendTo(select);
            }
            if (selectedOp) {
                select.val(selectedOp);
            }
            return select;
        },

        buildDataInput: function(field, op, value) {
            const type = (field && field.type) || "text";
            let input;
            if (type === "select" && field.choices && ["in", "ni"].indexOf(op) === -1) {
                input = $("<select>", {class: "form-control form-control-sm filter-rule-data"});
                $("<option>", {value: "", text: ""}).appendTo(input);
                for (const key of Object.keys(field.choices)) {
                    $("<option>", {value: key, text: field.choices[key]}).appendTo(input);
                }
            } else {
                const inputType = (NUMERIC_TYPES.indexOf(type) !== -1 && ["in", "ni"].indexOf(op) === -1)
                    ? "number" : "text";
                input = $("<input>", {type: inputType, class: "form-control form-control-sm filter-rule-data"});
                if (["in", "ni"].indexOf(op) !== -1) {
                    input.attr("placeholder", "comma separated");
                }
            }
            if (value !== undefined && value !== null) {
                input.val(value);
            }
            return input;
        },

        takesData: function(op) {
            for (const operation of this.operations) {
                if (operation.op === op) {
                    return operation.takesData;
                }
            }
            return true;
        },

        /* Field and operation both decide what the value input is, so it is rebuilt rather than
           reused when either changes */
        refreshRuleData: function(row) {
            const field = this.fieldByName($(".filter-rule-field", row).val());
            const op = $(".filter-rule-op", row).val();
            const cell = $(".filter-rule-data-cell", row);
            const existing = $(".filter-rule-data", cell).val();
            cell.empty();
            if (this.takesData(op)) {
                cell.append(this.buildDataInput(field, op, existing));
            }
        },

        addRule: function(rule) {
            const self = this;
            const row = $("<div>", {class: "form-row filter-rule mb-1"}).appendTo(this.rulesContainer);
            const ruleField = rule && rule.field;
            const field = this.fieldByName(ruleField) || (ruleField ? null : this.fields[0]);

            const fieldSelect = this.buildFieldSelect(ruleField || (field && field.field));
            const opSelect = this.buildOperationSelect(field, rule && rule.op);

            $("<div>", {class: "col-5"}).append(fieldSelect).appendTo(row);
            $("<div>", {class: "col-3"}).append(opSelect).appendTo(row);
            $("<div>", {class: "col-3 filter-rule-data-cell"}).appendTo(row);
            const actions = $("<div>", {class: "col-1"}).appendTo(row);

            if (this.readOnly) {
                $("select, input", row).prop("disabled", true);
            } else {
                $("<button>", {type: "button", class: "btn btn-sm btn-outline-secondary filter-rule-remove",
                               title: "Remove this filter", html: '<i class="fas fa-minus"></i>'})
                    .appendTo(actions)
                    .click(function() {
                        row.remove();
                        self.changed();
                    });

                fieldSelect.change(function() {
                    const newField = self.fieldByName($(this).val());
                    opSelect.replaceWith(self.buildOperationSelect(newField, opSelect.val()));
                    self.refreshRuleData(row);
                    self.changed();
                });
                row.on("change", ".filter-rule-op", function() {
                    self.refreshRuleData(row);
                    self.changed();
                });
                row.on("change", ".filter-rule-data", function() {
                    self.changed();
                });
            }

            this.refreshRuleData(row);
            if (rule && rule.data !== undefined) {
                $(".filter-rule-data", row).val(rule.data);
            }
            if (this.readOnly) {
                $("select, input", row).prop("disabled", true);
            }
            return row;
        },

        changed: function() {
            if (this.onChange) {
                this.onChange(this.getFilters());
            }
        },

        getFilters: function() {
            const rules = [];
            const self = this;
            $(".filter-rule", this.rulesContainer).each(function() {
                const field = $(".filter-rule-field", this).val();
                const op = $(".filter-rule-op", this).val();
                if (!field || !op) {
                    return;
                }
                const dataInput = $(".filter-rule-data", this);
                const data = self.takesData(op) ? (dataInput.val() || "") : "";
                rules.push({field: field, op: op, data: String(data)});
            });
            return {groupOp: this.groupOp, rules: rules};
        },

        setFilters: function(filters) {
            this.rulesContainer.empty();
            this.groupOp = (filters && filters.groupOp) || "AND";
            if (this.groupOpSelect) {
                this.groupOpSelect.val(this.groupOp);
            }
            const rules = (filters && filters.rules) || [];
            for (const rule of rules) {
                this.addRule(rule);
            }
            if (!rules.length && !this.readOnly) {
                this.addRule(null);
            }
        },

        hasRules: function() {
            return this.getFilters().rules.length > 0;
        },

        render: function() {
            const self = this;
            const container = this.container;
            container.empty().addClass("variantgrid-filter-builder");

            if (!this.fields.length) {
                container.text("This grid has no filterable columns.");
                return this;
            }

            const header = $("<div>", {class: "form-inline mb-2"}).appendTo(container);
            $("<label>", {class: "mr-2", text: "Match"}).appendTo(header);
            this.groupOpSelect = $("<select>", {class: "form-control form-control-sm mr-2"})
                .append($("<option>", {value: "AND", text: "all"}))
                .append($("<option>", {value: "OR", text: "any"}))
                .appendTo(header);
            $("<span>", {text: "of the following:"}).appendTo(header);

            this.rulesContainer = $("<div>", {class: "filter-rules"}).appendTo(container);

            if (this.readOnly) {
                this.groupOpSelect.prop("disabled", true);
                this.setFilters(null);
                return this;
            }

            this.groupOpSelect.change(function() {
                self.groupOp = $(this).val();
                self.changed();
            });

            const buttons = $("<div>", {class: "mt-2"}).appendTo(container);
            $("<button>", {type: "button", class: "btn btn-sm btn-outline-secondary mr-2",
                           html: '<i class="fas fa-plus"></i> Add filter'})
                .appendTo(buttons)
                .click(function() {
                    self.addRule(null);
                    self.changed();
                });
            $("<button>", {type: "button", class: "btn btn-sm btn-primary mr-2", text: this.applyLabel})
                .appendTo(buttons)
                .click(function() {
                    if (self.onApply) {
                        self.onApply(self.getFilters());
                    }
                });
            $("<button>", {type: "button", class: "btn btn-sm btn-outline-secondary", text: "Reset"})
                .appendTo(buttons)
                .click(function() {
                    self.setFilters(null);
                    const filters = self.getFilters();
                    if (self.onReset) {
                        self.onReset(filters);
                    } else if (self.onApply) {
                        self.onApply(filters);
                    }
                });

            this.setFilters(null);
            return this;
        }
    };

    return VariantGridFilterBuilder;
})();
