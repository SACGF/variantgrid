// @ts-check
/* global $, Urls, displayPhenotypeMatches, phenotypeMatchesToGridData, clearAutocompleteChoice */

/*
 * Phenotype text entry (patients/templatetags/phenotype_tags.py:phenotype_entry): highlights the ontology terms
 * matched in the text, lists them in a grid, and in edit mode re-matches as the user types or picks terms
 * from the autocompletes. Called once per tag on the page, everything scoped to its container.
 */
function initPhenotypeEntry(options) {
    const DONE_TYPING_INTERVAL = 2000;  // ms
    const AUTOCOMPLETE_IDS = ['id_hpo', 'id_omim', 'id_mondo', 'id_gene_symbol'];
    const phenotypeEditorTag = $("#phenotype-entry-" + options.flattenedUuid);
    const excludeString = options.phenotypeExcludeString || null;
    let typingTimer;
    let ignoreAutocompleteChangeEvent = false;

    function getPhenotypeInput() {
        return $("#" + options.entryFieldId, phenotypeEditorTag);
    }

    function showGrid() {
        $("#phenotype-entry-terms-grid-container", phenotypeEditorTag).show();
    }

    function ontologyTermFormatter(termText) {
        const m = termText.match(/(.*):(\d+) /);
        if (m) {
            const url = Urls.ontology_term(m[1] + "_" + m[2]);
            return $("<a />", {href: url, target: "_blank"}).text(termText).prop("outerHTML");
        }
        return $("<span />").text(termText).html();
    }

    // The rows come from the phenotype matcher rather than an endpoint, so this is a client side table
    function createTermGrid(gridData) {
        const selector = $("#phenotype-entry-terms-grid-" + options.flattenedUuid);

        if ($.fn.dataTable.isDataTable(selector)) {
            selector.DataTable().clear().rows.add(gridData).draw();
        } else {
            selector.DataTable({
                data: gridData,
                columns: [
                    {data: 'accession', title: 'Accession', render: ontologyTermFormatter},
                    {data: 'name', title: 'Name', render: $.fn.dataTable.render.text()},
                    {data: 'gene_symbols', title: 'Gene Symbols', render: $.fn.dataTable.render.text()},
                ],
                pageLength: 20,
                dom: '<"toolbar">rtBip',
                buttons: [{extend: 'csvHtml5', text: 'Download as CSV', filename: 'phenotypes'}],
            });
        }

        if (options.showGrid && gridData.length) {
            showGrid();
        }
    }

    function doneTyping() {
        const phenotypeText = getPhenotypeInput().val();
        const descriptionBox = $("#description-highlighted", phenotypeEditorTag);
        descriptionBox.siblings(".phenotype-messages").remove();
        descriptionBox.html('<i class="fa fa-spinner"></i>  Querying text for known terms');

        $.ajax({
            type: "POST",
            data: 'phenotype_text=' + encodeURIComponent(phenotypeText),
            url: Urls.phenotypes_matches(),
            success: function(phenotypeMatches) {
                displayPhenotypeMatches(descriptionBox, phenotypeText, phenotypeMatches, excludeString);
                const gridData = phenotypeMatchesToGridData(phenotypeMatches);
                if (gridData.length) {
                    $(".has-phenotypes-icon").show();
                    createTermGrid(gridData);
                } else {
                    $(".has-phenotypes-icon").hide();
                    $("#phenotype-entry-terms-grid-container", phenotypeEditorTag).hide();
                }
            },
        });
    }

    function startDescriptionUpdateTimer() {
        clearTimeout(typingTimer);
        if (getPhenotypeInput().val()) {
            typingTimer = setTimeout(doneTyping, DONE_TYPING_INTERVAL);
        }
    }

    function addToDescription(text) {
        const description = getPhenotypeInput();
        let existingText = description.val();

        if (existingText) {
            if ($("#new-terms-only", phenotypeEditorTag).is(":checked")) {
                if (existingText.indexOf(text) !== -1) {
                    return;
                }
            }
            existingText += '\n';
        }
        description.val(existingText + text);
        startDescriptionUpdateTimer();
    }

    function addTerm(value) {
        addToDescription(value);

        ignoreAutocompleteChangeEvent = true;
        for (const autocompleteId of AUTOCOMPLETE_IDS) {
            clearAutocompleteChoice($("#" + autocompleteId, phenotypeEditorTag));
        }
        ignoreAutocompleteChangeEvent = false;
    }

    function searchForText(regex, selector) {
        if (ignoreAutocompleteChangeEvent || !$(selector).val()) {
            return;
        }

        const rawText = $("option:selected", selector).text();
        const match = new RegExp(regex).exec(rawText);
        if (match) {
            addTerm(match[1].trim());
        }
    }

    function refreshExcludeWarning() {
        const warning = $("#phenotype-exclude-warning-" + options.flattenedUuid);
        const text = getPhenotypeInput().val() || "";
        warning.toggle(text.indexOf(excludeString) !== -1);
    }

    $(document).ready(function() {
        $(".show-phenotype-grid", phenotypeEditorTag).click(function(e) {
            e.preventDefault();
            $("#phenotype-entry-terms-grid-container", phenotypeEditorTag).slideDown();
        });

        const phenotypeMatches = options.phenotypeMatches;
        if (phenotypeMatches) {
            const descriptionBox = $("#description-highlighted", phenotypeEditorTag);
            displayPhenotypeMatches(descriptionBox, options.initialPhenotypeText, phenotypeMatches, excludeString);
            createTermGrid(phenotypeMatchesToGridData(phenotypeMatches));
            $(".has-phenotypes-icon").show();
        } else if (excludeString) {
            // The exclude string stops anything being matched and stored, but show the potential matches
            // (as if they removed it) by matching the text now
            const currentText = getPhenotypeInput().val() || "";
            if (currentText.indexOf(excludeString) !== -1) {
                doneTyping();
            }
        }

        if (options.edit) {
            getPhenotypeInput().keyup(startDescriptionUpdateTimer);

            $('#id_hpo', phenotypeEditorTag).change(function() {
                searchForText("HP:\\d{7} ([^\\(]*)", this);
            });
            $('#id_omim', phenotypeEditorTag).change(function() {
                searchForText("OMIM:\\d+ ([^\\(]*)", this);
            });
            $('#id_mondo', phenotypeEditorTag).change(function() {
                searchForText("MONDO:\\d+ ([^\\(]*)", this);
            });
            $('#id_gene_symbol', phenotypeEditorTag).change(function() {
                searchForText("(.*)", this);
            });
        }

        if (options.showGrid && options.hasPatientResults) {
            showGrid();
        }

        if (excludeString) {
            getPhenotypeInput().on("keyup change", refreshExcludeWarning);
            refreshExcludeWarning();
        }
    });
}
