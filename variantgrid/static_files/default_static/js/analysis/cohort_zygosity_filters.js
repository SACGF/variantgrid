// @ts-check
// analysis/templates/analysis/node_editors/cohort_zygosity_filters.html
function checkIfAllChecked(toggle, array) {
    let numChecked = 0;
    array.each(function() { numChecked += $(this).is(":checked"); });
    const allChecked = numChecked == array.length;
    $(toggle).prop('checked', allChecked);
}


function checkZygConfigWarning() {
    const zygosityTable = $(".zygosity-filter-table");
    let samplesWithAllZyg = 0;
    let samplesWithSomeZyg = 0;
    $(".cohort-node-zygosity-filter", zygosityTable).each(function() {
        const zygFilters = $(".zygosity_ref, .zygosity_het, .zygosity_hom, .zygosity_none", this);
        const numZyg = $("input:checked", zygFilters).length;
        if (numZyg) {
            if (numZyg === 4) {
                samplesWithAllZyg += 1;
            } else {
                samplesWithSomeZyg += 1;
            }
        }
    });

    if (samplesWithAllZyg || samplesWithSomeZyg) {
        $("#zyg-config-no-samples-warning").hide();
        if (samplesWithAllZyg && !samplesWithSomeZyg) {
            $("#zyg-config-unintuitive-no-zyg-call-warning").show();
        } else {
            $("#zyg-config-unintuitive-no-zyg-call-warning").hide();
        }

        const showInGrid = $(".show_in_grid", zygosityTable);
        if (!$("input:checked", showInGrid).length) {
            $("#zyg-config-no-visible-samples-warning").show();
        } else {
            $("#zyg-config-no-visible-samples-warning").hide();
        }


    } else {
        $("#zyg-config-no-samples-warning").show();
    }
}

function checkboxSetup(addClickHandler) {
    $(".toggle-row").each(function() {
        const rowInputs = $("input[type=checkbox][class!='toggle-row']", $(this).parents("tr"));
        checkIfAllChecked(this, rowInputs);

        if (addClickHandler) {
            $(this).click(function() {
                rowInputs.prop('checked', $(this).is(":checked"));
            });
        }
    });

    $(".toggle-column").each(function() {
        const columnName = $(this).attr("column_name");
        const columnInputs = $("input[type=checkbox]", "td." + columnName);
        checkIfAllChecked(this, columnInputs);

        if (addClickHandler) {
            $(this).click(function() {
                columnInputs.prop('checked', $(this).is(":checked"));
            });
        }
    });
}

$(document).ready(function() {
    checkboxSetup(true);
    checkZygConfigWarning();

    $("input[type=checkbox]", "table.zygosity-filter-table").click(function() {
        checkboxSetup(false);
        checkZygConfigWarning();
    });
});
