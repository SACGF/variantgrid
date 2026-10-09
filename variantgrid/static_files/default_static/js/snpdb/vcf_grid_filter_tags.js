// @ts-check
// snpdb/templates/snpdb/tags/vcf_grid_filter_tags.html
/* global ignoreFilterTable:writable, vcfGridParams:writable, tagContainer:writable */ // this file
// Plain assignments (not let/const) as may be reloaded multiple times in tabs
ignoreFilterTable = false;
vcfGridParams = {};
// The filter tags (one per tab) that this load added
tagContainer = $(".vcf-grid-filter-tags").not("[data-filter-setup]").attr("data-filter-setup", "true");

function vcfShowAll() {
    ignoreFilterTable = true;
    $("input:radio[name=genome_build_filter]:first", tagContainer).click();  // reset radio to "All"
    vcfClearAutoCompletes();
    ignoreFilterTable = false;
    vcfGridParams = {};
    filterTable();
}

function vcfClearAutoCompletes() {
    $("select[data-autocomplete-light-function=select2]", tagContainer).each(function() {
        clearAutocompleteChoice(this);
    });
}

// This is called by DataTables before each ajax call
function vcfGridDatatableFilter(data) {
    for (const [key, value] of Object.entries(vcfGridParams)) {
        data[key] = Array.isArray(value) ? JSON.stringify(value) : value;
    }
}

function filterTable() {
    if (ignoreFilterTable) {
        return;
    }

    const descriptionContainer = $("#vcf-grid-filter-description", tagContainer);
    if (!$.isEmptyObject(vcfGridParams)) {
        const descriptions = [];
        const genomeBuildName = vcfGridParams["genome_build_name"];
        if (genomeBuildName) {
            descriptions.push("Genome Build = " + genomeBuildName);
        }

        const project = vcfGridParams["project"];
        if (project) {
            descriptions.push("Project = " + project);
        }

        const variantsType = vcfGridParams["variants_type"];
        if (variantsType) {
            const VARIANTS_TYPE_LABELS = readJsonData("vcf-grid-filter-tags-data-" + tagContainer.attr("data-table-id")).variants_type_labels;
            const variantsTypeLabels = [];
            for (let i=0 ; i<variantsType.length ; i++) {
                const vt = variantsType[i];
                variantsTypeLabels.push(VARIANTS_TYPE_LABELS[vt]);
            }
            let variantsTypeDescriptions = variantsTypeLabels.join(", ");
            if (!variantsTypeDescriptions) {
                variantsTypeDescriptions = "(None selected)";
            }
            descriptions.push("Variants Type = " + variantsTypeDescriptions);
        }

        if (descriptions) {
            $("#vcf-grid-filter", tagContainer).html(descriptions.join(", "));
            descriptionContainer.show();
        }
    } else {
        descriptionContainer.hide();
    }

    $("#" + tagContainer.attr("data-table-id")).DataTable().ajax.reload();
}

$(document).ready(function() {
    $('#id_project_' + tagContainer.attr("data-table-id"), tagContainer).change(function() {
        vcfGridParams["project"] = $(this).val();
        filterTable();
    });

    $("input[name=genome_build_filter]", tagContainer).change(function() {
        const genomeBuildName = $(this).val();
        if (genomeBuildName) {
            vcfGridParams["genome_build_name"] = genomeBuildName;
        } else {
            delete vcfGridParams["genome_build_name"];
        }
        filterTable();
    });

    $("input[type=checkbox]", "#id_variants_type").change(function() {
        const variantsTypeList = [];
        $("input:checked", "#id_variants_type").each(function() {
            variantsTypeList.push($(this).val());
        });
        vcfGridParams["variants_type"] = variantsTypeList;
        filterTable();
    });

    $("#vcf-grid-filter-description", tagContainer).hide();
});
