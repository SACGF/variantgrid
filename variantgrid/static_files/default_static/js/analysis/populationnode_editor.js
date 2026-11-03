// @ts-check
// analysis/templates/analysis/node_editors/populationnode_editor.html
/* global ajaxForm, nodeEditorHooks */ // analysis/templates/analysis/node_editors/base_editor.html

function showGnomadPopulations() {
    $("#id_gnomad_af").prop("checked", false);
    $("#id_gnomad_popmax_af").prop("checked", false);
    $("#show-gnomad-populations-link").fadeOut();
    $("#gnomad-populations").slideDown();
}

function initPopulationNodeEditor() {
    $(document).ready(function () {
        ajaxForm($("form#population"));

        const gnomadPopulationSelector = $("input[type='checkbox']", "#id_gnomad_population");

        function deselectGnomadPopulations() {
            gnomadPopulationSelector.prop("checked", false);
            $("#show-gnomad-populations-link").fadeIn();
            $("#gnomad-populations").slideUp();
        }

        function selectedGnomadPopulations() {
            return gnomadPopulationSelector.filter(function () {
                return $(this).is(":checked");
            });
        }

        $("#id_gnomad_af").click(function () {
            if ($(this).is(":checked")) {
                deselectGnomadPopulations();
                $("#id_gnomad_popmax_af").prop("checked", false);
            }
        });

        $("#id_gnomad_popmax_af").click(function () {
            if ($(this).is(":checked")) {
                deselectGnomadPopulations();
                $("#id_gnomad_af").prop("checked", false);
            }
        });

        function changeUseInternalCounts() {
            const internalPop = $("#internal-population-container");
            if ($("#id_use_internal_counts").is(":checked")) {
                internalPop.slideDown();
            } else {
                internalPop.slideUp();
            }
        }

        // Show hidden things if set
        if (selectedGnomadPopulations().length) {
            $("#show-gnomad-populations-link").hide();
            $("#gnomad-populations").show();
        }

        $("#id_use_internal_counts").change(changeUseInternalCounts);
        changeUseInternalCounts();

        // A preset unticks every individual population, so put the picker away as the AF / popmax clicks do
        nodeEditorHooks.presetApplied = function () {
            if (!selectedGnomadPopulations().length) {
                deselectGnomadPopulations();
            }
        };
    });
}
