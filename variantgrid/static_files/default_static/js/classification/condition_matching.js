// @ts-check
// classification/templates/classification/condition_matching.html
/* global severityIcon */ // global.js
const conditionMatchingData = readJsonData("condition-matching-data");
const gene_symbol = conditionMatchingData.gene_symbol;
const basicHighlight = ["autosomal", "x-linked", "recessive", "dominant", gene_symbol];
const basicHighlightStr = basicHighlight.join("|");
const highlightRegex = new RegExp(`(${basicHighlightStr})`, 'ig');

function submitUserChoices() {
    const id = $('#mp-id').val();

    const selectedText = $('#mp-selected').val();
    const {terms, __} = parseTermText(selectedText);

    let finalText = terms.join(", ");
    let mode = 'N';
    if (terms.length >= 2) {
        mode = $('[name=mp-multimode]:checked').val();
        switch (mode) {
            case "U": finalText += "; uncertain"; break;
            case "C": finalText += "; co-occurring"; break;
            default: finalText += "; uncertain/co-occurring"; break;
        }
    }
    applyChange(parseInt(id, 10), terms, mode);
}

function parseTermText(text) {
    let termText = null;
    let join = 'N';

    const parts = text.split(";");
    if (parts.length === 2) {
        termText = parts[0].trim();
        join = parts[1].trim().toLowerCase();
        if (join.includes('/')) {
            join = 'N';
        } else if (join.includes('un')) {
            join = 'U';
        } else if (join.includes('co')) {
            join = 'C';
        }
    } else {
        termText = text;
    }
    const terms = termText.split(",").map(tt => tt.trim());

    return {
        terms,
        join
    };
}

function worstMessage(suggestion) {
    if (suggestion.messages) {
        const severity = {
            "i": 1,
            "s": 2,
            "w": 3,
            "e": 4
        };
        const maxMessage = suggestion.messages.reduce((best, message) => {
            const rank = severity[message.severity.toLowerCase()[0]] || 0;
            return (best === null || rank > best.rank) ? {rank, message} : best;
        }, null);
        return maxMessage ? maxMessage.message.severity : null;
    }
}

function toggleTerm(fullText, term) {
    const termArray = fullText.split(",").map(tt => tt.trim()).filter(tt => tt.length);
    const termDict = {};
    for (const term of termArray) {
        termDict[term] = true;
    }
    if (termDict[term]) {
        delete termDict[term];
    } else {
        termDict[term] = true;
    }
    const termList = Object.keys(termDict);
    updateMulti(termList.length);
    if (termList.length > 0) {
        termList.sort();
        return termList.join(", ");
    } else {
        return '';
    }
}

function highlightDescription(description) {
    if (!description) {
        return description;
    }
    return description.replace(highlightRegex, "<b>$1</b>");
}

function idSafe(id) {
    let idStr = id.replaceAll(/\W/g,"-");
    if (idStr.length == 0 || idStr.match(/^[\W].*$/)) {
        idStr = "x" + idStr;
    }
    return idStr;
}

function ontologyTr(onto) {
    const row = $('<tr>', {id: `result-${idSafe(onto.id)}`, class: 'onto-row'});
    const toggleButton = $('<button>', {class:'btn btn-outline-primary', text:'Toggle'}).click(() => {
        const mpSelected = $('#mp-selected');
        const fullText = mpSelected.val();
        const newText = toggleTerm(fullText, onto.id);
        mpSelected.val(newText);
        updateOntoRowSelections();
    });

    $('<td>', {class: 'toggle-cell', html:toggleButton}).appendTo(row);

     $('<td>', {class: 'term-cell', html:[
        $('<a>', {target:'blank', href:onto.url, text:onto.id}),
        $('<br>'),
        onto.title
    ]}).appendTo(row);

    const contexts = $('<td>', {class: 'relevance-cell'});
    contexts.append(renderContext(onto));

    contexts.appendTo(row);
    const description = $('<td>', {class:'description-cell text-muted' }).appendTo(row);
    description.append(highlightDescription(onto.definition));
    if (onto.sibling_count || onto.children_count) {
        $('<div>', {html:[
            $('<a>', {text: `Siblings: ${onto.sibling_count}, Children: ${onto.children_count}`})
        ]}).appendTo(description);
    }

    return row;
}

function updateOntoRowSelections() {
    $('.onto-row').removeClass('selected');
    const allSelected = $('#mp-selected').val().split(',').map(t => t.trim()).filter(t => t.length > 0);
    for (const selected of allSelected) {
        $(`#result-${idSafe(selected)}`).addClass('selected');
    }
}

function renderOntologies(ontos) {
    const tbody = $('#mp-results tbody');
    tbody.empty();

    const errors = ontos.errors;
    for (const error of errors) {
        tbody.append($('<tr>', {html:
                $('<td>', {'colspan': 4, html: [severityIcon('error'), $('<span>', {text:error})]})
        }));
    }

    const terms = ontos.terms;
    for (const onto of terms) {
        ontologyTr(onto).appendTo(tbody);
    }
}

function updateMulti(count) {
    if (count > 1) {
        $('#mp-multimode').show();
        $('#mp-singlemode').hide();
    } else {
        $('#mp-multimode').hide();
        $('#mp-singlemode').show();
    }
}

function searchMondo(data) {
    const searchText = data.searchText;
    const geneSymbol = data.geneSymbol;
    const selected = data.selected;
    updateMulti(selected ? selected.split(",").length : 0);

    const url = Urls.api_mondo_search();
    const resultsDom = $("#mp-results");
    resultsDom.LoadingOverlay('show');

    $.ajax({
        headers: {
            'Accept': 'application/json',
            'Content-Type': 'application/json'
        },
        data: {
            "search_term": searchText,
            "gene_symbol": geneSymbol,
            "selected": selected
        },
        url: url,
        type: 'GET',
        error: (call, status, text) => {
            resultsDom.LoadingOverlay('hide');
            resultsDom.find('tbody').empty().append($('<tr>', {html:
                    $('<td>', {'colspan': 4,  html: [severityIcon('critical'), $('<span>', {text:"Error retrieving MONDO suggestions."})]})
            }));
            console.log(status);
            console.log(text);
        },
        success: async (results) => {
            resultsDom.LoadingOverlay('hide');
            renderOntologies(results);
            updateOntoRowSelections();
        }
    });
}

function renderContext(onto) {
    const allContexts = $('<div>');
    if (onto.direct_reference) {
        $('<div>', {html: [
                $(`<i class="fas fa-circle text-dark"></i>`),
                $('<span>', {text: "Referenced in Condition Text" })
        ]}).appendTo(allContexts);
    }
    if (onto.text_search) {
         $('<div>', {html: [
                $(`<i class="fas fa-circle text-secondary"></i>`),
                $('<span>', {text: "Text search" })
        ]}).appendTo(allContexts);
    }
    if (onto.gene_relationships && onto.gene_relationships.length) {
        const html = [
                $(`<i class="fas fa-circle text-success"></i>`),
                $('<span>', {text: "Established gene relationship"}),
                " "
            ];
        let first = true;
        for (const relationship of onto.gene_relationships) {
            if (!first) {
                html.push(", ");
            } else {
                first = false;
            }
            if (relationship.relation === "panelappau") {
                const extra = relationship.extra || {};
                const phenoEvidences = extra.phenotypes_and_evidence || [];
                const allEvidencesSet = new Set();
                for (const pe of phenoEvidences) {
                    const evidences = pe.evidence || [];
                    if (evidences.indexOf("Expert Review Green") !== -1) {
                        for (const evidence of (pe.evidence || [])) {
                            if (!evidence.startsWith("Expert Review")) {
                                allEvidencesSet.add(evidence);
                            }
                        }
                    }
                }
                const allEvidences = Array.from(allEvidencesSet);
                allEvidences.sort();

                let tooltip = `Expert Review Green Sources: ${allEvidences.join(', ')}`;
                if (relationship.via) {
                    tooltip += `<br/>Via related term: ${relationship.via}`;
                }
                html.push( $('<span>', {class: 'text-muted', text: 'PanelApp AU', 'data-toggle':'popover', 'title': 'Gene Relationship', 'data-content': tooltip}) );
            } else if (relationship.source == "gencc_file") {
                const extra = relationship.extra || {};
                const sources = extra.sources || [];
                let tooltip = "";
                for (const source of extra.sources) {
                    if (tooltip.length) {
                        tooltip += "<br/><br/>";
                    }
                    tooltip += `${source.submitter || 'Unknown source'} : ${source.mode_of_inheritance || '-'} : ${source.gencc_classification}`;
                }
                html.push( $('<span>', {class: 'text-muted', text: 'GenCC', 'data-toggle':'popover', 'title': 'Gene Relationship', 'data-content': tooltip}) );
            } else if (relationship.source === "hpo_disease") {
                const tooltip = `Via related term: ${relationship.via}`;
                html.push( $('<span>', {class: 'text-muted', text: 'DEPRECATED', 'data-toggle':'popover', 'title': 'Gene Relationship', 'data-content': "This relationship is from a deprecated file. Other relationships listed are fine."}) );
            } else if (relationship.source === "mondo_file") {
                const tooltip = `Relationship: ${relationship.relation}`;
                html.push( $('<span>', {class: 'text-muted', text: "MONDO", 'data-toggle':'popover', 'title': 'Gene Relationship', 'data-content': tooltip}) );
            } else {
                const tooltip = `Relationship: ${relationship.relation}`;
                html.push( $('<span>', {class: 'text-muted', text: relationship.source, 'data-toggle':'popover', 'title': 'Gene Relationship', 'data-content': tooltip}) );
            }
        }
        $('<div>', {
            html: html
        }).appendTo(allContexts);
    }
    return allContexts;
}

function triggerSearchMondo() {
    const searchText = $('#mp-search-text').val();
    let geneSymbol = $('#mp-gene-symbol').text();
    const selected = $('#mp-selected').val();
    if (geneSymbol === "-") {
        geneSymbol = "";
    }

    searchMondo({
        searchText,
        geneSymbol,
        selected
    });
}

function canSuggest(suggestion) {
    if (!suggestion.is_applied && !suggestion.info_only && suggestion.terms.length >= 1) {
        for (const message of suggestion.messages) {
            if (message.severity.toLowerCase().startsWith("e")) {
                return false;
            }
        }
        return true;
    }
    return false;
}

function renderSuggestion(suggestion, id) {
    const dom = $('<div>');
    if (suggestion.terms && suggestion.terms.length) {
        dom.addClass('suggestion');
        if (suggestion.info_only) {
            dom.addClass('info-only');
            dom.append($('<div>', {class: 'text-muted', text:'Text match to inform gene level suggestions:'}));
        }
        let termStyle = '';

        if (canSuggest(suggestion)) {

            let suggest = true;
            const worst = worstMessage(suggestion);
            if (worst != null && ["s", "i"].indexOf(worst[0].toLowerCase()) == -1) {
                // we have a message and it's not success or info
                suggest = false;
            }
            const input = $('<input>', {
                type: 'checkbox',
                class: 'form-check-input suggestion-checkbox',
                click: () => {
                    checkSuggestionCount();
                }
            });
            if (suggest) {
                input.attr('checked', 'checked');
            }

            $('<label>', {
                style: 'margin-left:1.25rem',
                    class: 'form-check-label',
                    html: [
                    input,
                    'Selected',
                ]
            }).appendTo(dom);
            termStyle = 'opacity:0.5';
        }
        for (const term of suggestion.terms) {
            // TODO use generated URL
            const termDom = $('<div>', {
                style: termStyle, html: [
                    $('<a>', {
                        href: `/ontology/term/${term.id.replace(':', '_')}`,
                        text: term.id,
                        target: '_blank',
                        title: term.definition,
                        'data-toggle': 'tooltip'
                    }),
                    ' ',
                    term.name
                ]
            });
            dom.append(termDom);
        }
        if (suggestion.terms.length > 1) {
            let joiner_text = 'Combination type undecided';
            switch (suggestion.joiner) {
                case 'U':
                    joiner_text = 'Uncertain';
                    break;
                case 'C':
                    joiner_text = 'Co-occurring';
                    break;
            }
            $('<div>', {text: joiner_text, class:'font-italic'}).appendTo(dom);
        }
        if (suggestion.user) {
            $('<div>', {text: `Set by ${suggestion.user.username}`, class:'text-secondary'}).appendTo(dom);
        }
    }
    if (suggestion.messages && suggestion.messages.length) {
        if (!suggestion.is_applied) {
            const maxSev = worstMessage(suggestion);
            dom.addClass(maxSev);
        }
        dom.addClass('suggestion');
        for (const message of suggestion.messages) {
            $('<div>', {html:[
                    severityIcon(message.severity),
                    ' ',
                    message.text
            ]}).appendTo(dom);
        }
    }
    return dom;
}

function checkSuggestionCount() {
    const checkedCount = $('.suggestion-checkbox:checked').length;
    const button = $('#approve-suggestions');
    button.text(`Apply ${checkedCount} selected suggestion${checkedCount === 1 ? '' : 's'}`);
    if (checkedCount === 0) {
        button.hide();
    } else {
        button.show();
    }
}

function approveSuggestions() {
    const checked = $('.suggestion-checkbox:checked').parents('.condition-match-row');
    const changes = [];
    for (let placeholder of checked) {
        placeholder = $(placeholder);
        const cid = parseInt(placeholder.attr('data-id'));
        const terms = placeholder.attr('data-terms').split(",");
        const joiner = placeholder.attr('data-join');
        changes.push({
            'ctm_id': cid,
            'terms': terms,
            'joiner': joiner
        });
    }
    if (changes.length) {
        console.log(changes);
        applyChanges(changes);
    }
}

function applyChange(id, terms, joiner) {
    applyChanges([{
        'ctm_id': id,
        'terms': terms,
        'joiner': joiner
    }]);
}

function applyChanges(changes) {
    $('#condition-list').LoadingOverlay('show');
    $.ajax({
        headers: {
            'Accept': 'application/json',
            'Content-Type': 'application/json'
        },
        data: JSON.stringify({"changes": changes}),
        url: Urls.condition_text_matching_api(conditionMatchingData.condition_text_id),
        type: 'POST',
        error: (call, status, text) => {
            $('#condition-list').LoadingOverlay('hide');
            alert('There was an error updating this term');
            console.log(status);
            console.log(text);
        },
        success: async (results) => {
            applySuggestionUpdates(results);
            $('#condition-list').LoadingOverlay('hide');
        }
    });
}

function updateSuggestions() {
    $('#condition-list').LoadingOverlay('show');
    $.ajax({
        headers: {
            'Accept': 'application/json',
            'Content-Type': 'application/json'
        },
        url: Urls.condition_text_matching_api(conditionMatchingData.condition_text_id),
        type: 'GET',
        error: (call, status, text) => {
            $('#condition-list').LoadingOverlay('hide');
            alert('There was an error retrieving terms');
            console.log(status);
            console.log(text);
        },
        success: async (results) => {
            applySuggestionUpdates(results, true);
            $('#condition-list').LoadingOverlay('hide');
        }
    });
}

function parentTerms(rowDom) {
    const parentId = rowDom.attr('data-parent-id');
    if (parentId) {
        const parent = $(`#condition-match-${parentId}`);
        if (parent.is("[data-applied]")) {
            const terms = parent.attr('data-terms').split(",");
            if (terms.length) {
                return terms.map(t => t.trim()).join(", ");
            }
        }
        return parentTerms(parent);
    }
    return null;
}

function applySuggestionUpdates(results, complete) {
    const errors = results.errors;
    if (errors && errors.length) {
        window.alert(errors.join("\n") + "\n- rejecting change.");
    }

    const suggestions = results.suggestions;
    $("#count_outstanding").text(results.count_outstanding);
    if (results.count_outstanding == 0) {
        $("#no_outstanding").html(severityIcon("S"));
    } else {
        $("#no_outstanding").html(severityIcon("W"));
    }

    if (complete) {
        $('.condition-match-row').removeAttr('data-terms').removeAttr('data-join').removeAttr('data-applied').removeAttr('data-info');
    }
    for (const suggestion of suggestions) {
        const wrapper = $(`#condition-match-${suggestion.id}`);
        wrapper.removeAttr('data-terms').removeAttr('data-join').removeAttr('data-applied').removeAttr('data-info');

        const valuesDom = wrapper.find('.condition-match-values').empty();
        const suggestionDom = wrapper.find('.condition-match-suggestion').empty();

        wrapper.attr('data-terms', suggestion.terms.map(term => term.id).join(","));
        wrapper.attr('data-join', suggestion.joiner);
        if (suggestion.is_applied) {
            valuesDom.html(renderSuggestion(suggestion));
            if (suggestion.terms && suggestion.terms.length) {
                wrapper.attr('data-applied', 'true');
            }
        } else {
            if (suggestion.info_only) {
                wrapper.attr('data-info', 'true');
            }
            suggestionDom.html(renderSuggestion(suggestion));
        }
    }

    const inheritings = $('.condition-match-row:not([data-applied])');
    $('.inherit').remove();
    for (let inheriting of inheritings) {
        inheriting = $(inheriting);
        const inheritedTerms = parentTerms(inheriting);
        const valuesDom = inheriting.find('.condition-match-values');
        const parentId = inheriting.attr('data-parent-id');
        if (parentId) {
            const parent = $(`#condition-match-${parentId}`).attr('data-label');
            if (inheritedTerms) {
                valuesDom.prepend($('<span>', {class:'inherit text-muted font-italic', text:inheritedTerms, title:`Inherits from ${parent}`}));
            } else {
                valuesDom.prepend($('<span>', {class:'inherit text-muted font-italic', html:'none', title:`Inherits from ${parent}`}));
            }
        } else {
            valuesDom.prepend($('<span>', {class:'inherit text-muted font-italic', text:`-`, title:`No top level terms provided`}));
        }
    }
    const topLevel = $('.condition-match-row:not([data-parent-id])');
    const topLevelDataTerms = topLevel.attr('data-terms');
    const button = topLevel.find('.mondo-picker');
    // if we already have data-terms on top level, don't provide warning
    if ((topLevel.attr('data-terms') || '').trim().length && !topLevel.attr('data-info')) {
        button.css('opacity','1.0');
        button.removeAttr('title');
        button.tooltip('dispose');
    } else {
        button.css('opacity','0.3');
        button.attr('title', 'It is recommended you choose gene level terms using the buttons below');
        button.tooltip({html:true, trigger : 'hover'});
    }

    checkSuggestionCount();
}

$(document).ready(() => {
    updateSuggestions();
    checkSuggestionCount();

    // open MODAL
    // TODO : Make it so you can't click this button until the server
    $('.mondo-picker').click(function () {
        const mondoButton = $(this);
        const mondoRow = mondoButton.parents('.condition-match-row');
        const geneSymbol = mondoRow.attr('data-gene-symbol');
        const label = mondoRow.attr('data-label');
        const inheritance = mondoRow.attr('data-mode-of-inheritance');
        if (geneSymbol === '') {
            $('#mp-gene-symbol').text('-').addClass('no-value');
        } else {
            $('#mp-gene-symbol').text(geneSymbol).removeClass('no-value');
        }
        const recordId = mondoRow.attr('data-id');
        const selected = mondoRow.attr('data-terms');
        const joiner = mondoRow.attr('data-join');
        // let parsedMatch = parseTermText(selected);
        $(`[name=mp-multimode][value=${joiner}]`).trigger('click');

        $('#mp-id').val(recordId);
        $('#mp-search-text').val(conditionMatchingData.normalized_text.trim());
        $('#mp-selection-title').text(label);
        $('#mp-selected').val(selected);
        if (inheritance === 'N/A') {
            $('#mp-mode-of-inheritance').html($('<span>', {class: 'no-value', text: 'N/A'}));
        } else if (inheritance === 'None') {
            $('#mp-mode-of-inheritance').html($('<span>', {class: 'no-value', text: 'Not Specified'}));
        } else {
            $('#mp-mode-of-inheritance').text(inheritance);
        }
        $('#mondoModal').modal('show');

        triggerSearchMondo();
        return false;
    });

    $('#mp-search-button').click(() => {
        triggerSearchMondo();
    });

    $('#mp-selected').change(() => {
        updateOntoRowSelections();
    });
});
