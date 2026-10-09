// @ts-check
// snpdb/templates/snpdb/tags/model_fields_version_diff.html
function getFieldDiffClass(a_val, b_val) {
    let fieldDiffClass = null;
    if (a_val != b_val) {
        if (!a_val) {
            if (b_val) {
                fieldDiffClass = "addition";
            }
        } else if (!b_val) {
            if (a_val) {
                fieldDiffClass = "deletion";
            }
        } else {
            fieldDiffClass = "modification";
        }
    }
    return fieldDiffClass;
}

function create_select_from_tuples(selector, list_of_tuples) {
    const sel = $('<select>').appendTo(selector);
    $(list_of_tuples).each(function() {
        sel.append($("<option>").attr('value', this[0]).text(this[1]));
    });
    return sel;
}

const models_by_version = readJsonData("model-fields-version-diff-data").model_dicts_by_version;
const modelVersionIds = Object.keys(models_by_version).map(Number);
const versions = readJsonData("model-fields-version-diff-data").versions;

function showModelDiff(aId, bId) {
    const table = $("#version-comparison");
    $("tr.field-row", table).remove();

    // test if they're the same...
    const sm = $("#select-messages");
    if (aId == bId) {
        sm.text("Can't diff the same versions.");
        sm.show();
        return;
    } else {
        sm.hide();
    }

    const a = models_by_version[aId];
    const b = models_by_version[bId];

    const model_fields = Object.keys(a).sort();
    for (let i=0 ; i<model_fields.length ; ++i ) {
        const f = model_fields[i];
        const a_val = a[f];
        const b_val = b[f];
        const fieldDiffClass = getFieldDiffClass(a_val, b_val);
        const row = $("<tr>", {"id" : f + "-row", class: "field-row " + fieldDiffClass});
        row.append($("<th>", {"class" : "field-name"}).text(f));
        row.append($("<td>", {"id" : f + "-" + aId}).text(a_val));
        row.append($("<td>", {"id" : f + "-" + bId}).text(b_val));
        table.append(row);
    }
}


function add_selects() {
    const a_sel = create_select_from_tuples($("#a-version-header"), versions);
    const b_sel = create_select_from_tuples($("#b-version-header"), versions);

    const selectChanged = function() {
        showModelDiff(a_sel.val(), b_sel.val());
    };

    const firstId = Math.min.apply(null, modelVersionIds);
    const lastId = Math.max.apply(null, modelVersionIds);

    a_sel.val(firstId);
    b_sel.val(lastId);

    a_sel.change(selectChanged);
    b_sel.change(selectChanged);
    selectChanged();

}


$(document).ready(function() {
    add_selects();
});
