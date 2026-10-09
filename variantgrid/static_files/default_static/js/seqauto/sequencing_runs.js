// @ts-check
// seqauto/templates/seqauto/sequencing_runs.html
/* global ignoreAutcompleteChangeEvent:writable */ // never declared - an implicit global, only ever written here
let sequencingRunFilters = {};

function showAll() {
    clearAutoCompletes();
    sequencingRunFilters = {};
    filterTable();
}

function clearAutoCompletes(leave_autocomplete_term_type) {
    // Clear the autcompletes (temp disabling change handlers)
    const AUTOCOMPLETE_IDS = {
        'enrichment_kit': 'id_enrichment_kit',
    };

    ignoreAutcompleteChangeEvent = true;
    for (const t in AUTOCOMPLETE_IDS) {
        if (t == leave_autocomplete_term_type) {
            continue; // leave this one
        }
        const autocomplete_id = AUTOCOMPLETE_IDS[t];
        clearAutocompleteChoice("#" + autocomplete_id);
    }
    ignoreAutcompleteChangeEvent = false;
}

// This is called by DataTables before each ajax call
function sequencingRunDatatableFilter(data) {
    for (const [key, value] of Object.entries(sequencingRunFilters)) {
        data[key] = value;
    }
}

function filterTable() {
    const descriptionContainer = $("#sequencing-runs-filter-description");
    if (!$.isEmptyObject(sequencingRunFilters)) {
        const descriptions = [];
        if (sequencingRunFilters["enrichment_kit_id"]) {
            descriptions.push("Enrichment Kit");
        }
        if (descriptions) {
            $("#sequencing-runs-filter").html(descriptions.join(" "));
            descriptionContainer.show();
        }
    } else {
        descriptionContainer.hide();
    }

    $("#sequencing-runs-datatable").DataTable().ajax.reload();
}

$(document).ready(function() {
    $('#id_enrichment_kit').change(function() {
        sequencingRunFilters["enrichment_kit_id"] = $(this).val();
        filterTable();
    });
    $('#id_sequencing_run').change(function () {
        const sequencingRunId = $("#id_sequencing_run").val();
        window.location = Urls.view_sequencing_run(sequencingRunId);
    });
});
