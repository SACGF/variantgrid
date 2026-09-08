/*
 * What the samples ticked on the VCF / Cohort page can be turned into (issue private#933).
 *
 * One bar for the lot: how many samples are selected, and only the action that many can actually
 * launch - a Duo/Trio/Quad analysis, or a cohort of the selection. The hint says what to select
 * to get the others.
 */
const SampleSelectionActions = (function () {
    let config = null;

    function selectedSampleIds() {
        return $("input.sample-select:checked").map(function () {
            return Number($(this).attr("sample_id"));
        }).get();
    }

    function familyAnalysisUrl(button, sampleIds) {
        return Urls[$(button).data("url-name")].apply(null, [config.cohortId].concat(sampleIds));
    }

    function createSubCohort(sampleIds) {
        $.ajax({
            type: "POST",
            url: config.subCohortUrl,
            data: {sample_id_list: JSON.stringify(sampleIds)},
            success: (data) => window.location = Urls.view_cohort(data["cohort_id"]),
        });
    }

    /** "2 for a Duo, 3 for a Trio or any 2 or more for a cohort" - everything but what's selected now */
    function hintText(numSelected) {
        const options = [];
        $(".family-analysis-action").each(function () {
            const numSamples = Number($(this).data("num-samples"));
            if (numSamples !== numSelected) {
                options.push(`${numSamples} for a ${$(this).data("label")}`);
            }
        });
        if (numSelected < 2) {
            options.push("any 2 or more for a cohort");
        }
        if (!options.length) {
            return "";
        }
        const last = options.pop();
        return "Select " + (options.length ? options.join(", ") + " or " + last : last);
    }

    function update() {
        const numSelected = selectedSampleIds().length;
        $("#sample-actions .selected-number").text(numSelected);
        $("#sample-actions .selected-label").text(numSelected === 1 ? "sample selected" : "samples selected");
        $(".family-analysis-action").each(function () {
            $(this).toggle(Number($(this).data("num-samples")) === numSelected);
        });
        $("#create-sub-cohort-action").toggle(numSelected >= 2);
        $("#sample-actions .sample-actions-hint").text(hintText(numSelected));
    }

    function init(cfg) {
        config = cfg;
        $("#sample-actions").on("click", ".family-analysis-action", function () {
            window.location = familyAnalysisUrl(this, selectedSampleIds());
        });
        $("#create-sub-cohort-action").click(() => createSubCohort(selectedSampleIds()));
        // The cohort membership table redraws itself, so watch the checkboxes from the document
        $(document).on("change", "input.sample-select", update);
        $(document).on("change", "input.sample-select-all", function () {
            $("input.sample-select").prop("checked", $(this).is(":checked"));
            update();
        });
        update();
    }

    return {
        init: init,
        update: update,
    };
})();
