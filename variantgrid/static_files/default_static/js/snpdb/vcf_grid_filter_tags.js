// @ts-check
// snpdb/templates/snpdb/tags/vcf_grid_filter_tags.html
/* global ignoreFilterTable:writable, vcfGridParams:writable, tagContainer:writable */ // this file
/* global vcfGridTableId:writable, vcfGridVariantsTypeLabels:writable */ // this file
// Plain assignments (not let/const) as the tag may be reloaded multiple times in tabs - the last one set up
// is the one the filter functions act on
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
            const variantsTypeLabels = [];
            for (let i=0 ; i<variantsType.length ; i++) {
                const vt = variantsType[i];
                variantsTypeLabels.push(vcfGridVariantsTypeLabels[vt]);
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

    $("#" + vcfGridTableId).DataTable().ajax.reload();
}

function setupVcfGridFilterTags(tableId, variantsTypeLabels) {
    ignoreFilterTable = false;
    vcfGridParams = {};
    vcfGridTableId = tableId;
    vcfGridVariantsTypeLabels = variantsTypeLabels;
    tagContainer = $(".vcf-grid-filter-tags-" + tableId);

    // Called from inside the tag's own container, before the controls below it are parsed
    $(document).ready(function() {
        $('#id_project_' + tableId, tagContainer).change(function() {
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
}
