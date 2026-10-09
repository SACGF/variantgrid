// @ts-check
// snpdb/templates/snpdb/tags/vcf_grid_filter_tags.html
// The tag is on both the Samples and VCF tabs of the Data page, so each instance keeps its state in a closure,
// reached through its own container (keyed by table_id) rather than page globals

function vcfGridFilterContainer(tableId) {
    return $(".vcf-grid-filter-tags-" + tableId);
}

// The DataTables ajax `data` hook for tableId - data-datatable-data="vcfGridDatatableFilter('samples-datatable')"
function vcfGridDatatableFilter(tableId) {
    return function(data) {
        const datatableFilter = vcfGridFilterContainer(tableId).data("vcfGridDatatableFilter");
        if (datatableFilter) {
            datatableFilter(data);
        }
    };
}

function initVcfGridFilter(tableId, options) {
    const tagContainer = vcfGridFilterContainer(tableId);
    const variantsTypeLabels = options.variantsTypeLabels;
    let ignoreFilterTable = false;
    let vcfGridParams = {};

    function filterTable() {
        if (ignoreFilterTable) {
            return;
        }

        const descriptionContainer = $(".vcf-grid-filter-description", tagContainer);
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
                const variantsTypeDescriptions = variantsType.map(vt => variantsTypeLabels[vt]).join(", ");
                descriptions.push("Variants Type = " + (variantsTypeDescriptions || "(None selected)"));
            }

            $(".vcf-grid-filter", tagContainer).html(descriptions.join(", "));
            descriptionContainer.show();
        } else {
            descriptionContainer.hide();
        }

        $("#" + tableId).DataTable().ajax.reload();
    }

    function showAll() {
        ignoreFilterTable = true;
        $("input:radio[name=genome_build_filter]:first", tagContainer).click();  // reset radio to "All"
        $("select[data-autocomplete-light-function=select2]", tagContainer).each(function() {
            clearAutocompleteChoice(this);
        });
        ignoreFilterTable = false;
        vcfGridParams = {};
        filterTable();
    }

    tagContainer.data("vcfGridDatatableFilter", function(data) {
        for (const [key, value] of Object.entries(vcfGridParams)) {
            data[key] = Array.isArray(value) ? JSON.stringify(value) : value;
        }
    });

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

        const variantsTypeCheckboxes = $("#id_variants_type input[type=checkbox]", tagContainer);
        variantsTypeCheckboxes.change(function() {
            vcfGridParams["variants_type"] = variantsTypeCheckboxes.filter(":checked").map(function() {
                return $(this).val();
            }).get();
            filterTable();
        });

        $(".vcf-grid-filter-show-all", tagContainer).click(function(event) {
            event.preventDefault();
            showAll();
        });

        $(".vcf-grid-filter-description", tagContainer).hide();
    });
}
