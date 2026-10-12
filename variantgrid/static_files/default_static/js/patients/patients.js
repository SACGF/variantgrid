// @ts-check
// patients/templates/patients/patients.html
/* global setupPatientSamplenode */ // patient_samplenode.js
const PATIENT_TABLE = "#patient-datatable";
let ignoreAutcompleteChangeEvent = false;
let termFilter = null;
let filterFields = {};

function patientDatatableFilter(data) {
    if (termFilter) {
        data.term_type = termFilter.term_type;
        data.term = termFilter.term;
    }
    Object.assign(data, filterFields);
}

function reloadPatientTable() {
    $(PATIENT_TABLE).DataTable().ajax.reload();
}

function changeFilter(field, value) {
    if (value) {
        filterFields[field] = value;
    } else {
        delete filterFields[field];
    }
    reloadPatientTable();
}

function showCreatePatientForm() {
    $("#new-patient-link").hide();
    $("#new-patient-form-container").slideDown();
}

function resetGrid() {
    $("#grid-filter-message").hide();
    clearAutoCompletes();
    filterFields = {};
    termFilter = null;
    reloadPatientTable();
}

function clearAutoCompletes(leave_autcomplete_term_type) {
    // Clear the autcompletes (temp disabling change handlers)
    const AUTOCOMPLETE_IDS = {
        'HP': 'id_hpo',
        'OMIM': 'id_omim',
        'HGNC': 'id_gene_symbol',
        'MONDO': 'id_mondo',
    };

    ignoreAutcompleteChangeEvent = true;
    for (const t in AUTOCOMPLETE_IDS) {
        if (t === leave_autcomplete_term_type) {
            continue; // leave this one
        }
        const autocomplete_id = AUTOCOMPLETE_IDS[t];
        clearAutocompleteChoice("#" + autocomplete_id);
    }
    ignoreAutcompleteChangeEvent = false;
}

function filterPatientGrid(term_type, value, leave_autcomplete_term_type) {
    const msgBox = $("#grid-filter-message");
    if (value) {
        const url = Urls.ontology_term_text(term_type, value);
        msgBox.html("Filtering grid to <span class='" + term_type.toLowerCase() + "'><a href='" + url + "' target='_blank'>" + value + "</a></span>... <a href='javascript:resetGrid()'>show all</a>");
        msgBox.show();
        termFilter = {term_type: term_type, term: value};
        reloadPatientTable();
    } else {
        msgBox.html("Error: could not search on <span class='" + term_type.toLowerCase() + "'> value='" + value + "'</span>... <a href='javascript:resetGrid()'>show all</a>");
        msgBox.show();
    }

    clearAutoCompletes(leave_autcomplete_term_type);
}

function filterPatientGridFromBarClick(term_type, data) {
    const barClicked = data.points[0];
    const value = barClicked.y;
    filterPatientGrid(term_type, value);
}

function plotlyHPOClickHandler(data) {
    filterPatientGridFromBarClick('HP', data);
}

function plotlyOMIMClickHandler(data) {
    filterPatientGridFromBarClick('OMIM', data);
}

$(document).ready(function() {
    // The term links and collapsed groups are made as the table draws, so delegate from the document
    $(document).on('click', '.ontology-terms-container .collapsed-term', expandCollapsedOntologyTerm);
    $(document).on('click', '.grid-term-link:not(.collapsed-term)', function() {
        filterPatientGrid($(this).attr("term_type"), $(this).attr("term"));
    });

    function loadPatientOnSearchSelect() {
        $('#id_patient').change(function () {
            const patientId = $(this).val();
            if (patientId) {
                window.location = Urls.view_patient(patientId);
            }
        });
    }

    function filterGridOnTermSearchSelect() {
        function searchForText(term_type, regex, selector) {
            if (ignoreAutcompleteChangeEvent) {
                return;
            }

            if (!$(selector).val()) { // empty - show all on grid
                resetGrid();
                return;
            }

            const raw_text = $("option:selected", selector).text();
            const match = new RegExp(regex).exec(raw_text);
            if (match) {
                const value = match[1].trim();
                // pass in 3rd param to tell it we're from autocomplete box (and leave it)
                filterPatientGrid(term_type, value, term_type);
            }
        }

        $('#id_hpo').change(function () {
            const regex = "HP:\\d{7} ([^\\(]*)";
            searchForText('HP', regex, this);
        });
        $('#id_omim').change(function () {
            const regex = "OMIM:\\d+ ([^\\(]*)";
            searchForText('OMIM', regex, this);
        });
        $('#id_mondo').change(function () {
            const regex = "MONDO:\\d+ ([^\\(]*)";
            searchForText('MONDO', regex, this);
        });
        $('#id_hgnc').change(function () {
            const regex = "HGNC:\\d+ ([^\\(]*)";
            searchForText('HGNC', regex, this);
        });
    }

    const form = $("form#patient-form");
    $("#id_first_name", form).change(function () {
        changeFilter('first_name', $(this).val());
    });
    $("#id_last_name", form).change(function () {
        changeFilter('last_name', $(this).val());
    });

    $("#id_sex", form).change(function () {
        let data = $(this).val();
        if (data === 'U') {
            data = null; // Use blank to not search when unknown
        }
        changeFilter('sex', data);
    });
    $("button#reset-form").click(function () {
        $("#patient-success").remove();
        filterFields = {};
        reloadPatientTable();
    });

    $('form#create-patient-form').ajaxForm({});
    setupPatientSamplenode({name: ''});
    loadPatientOnSearchSelect();
    filterGridOnTermSearchSelect();

    $("button#create-new-patient").click(showCreatePatientForm);
});
