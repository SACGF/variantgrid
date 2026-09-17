/*
 * One span per ontology per matched stretch of text - "osteogenesis imperfecta" matches 24 OMIM terms, which
 * were 24 nested spans. The tooltip lists every term the span stands for.
 */
function groupPhenotypeMatchesByOntology(phenotypeMatches) {
    const groups = new Map();
    for (const pm of phenotypeMatches) {
        const key = `${pm.ontology_service}:${pm.offset_start}:${pm.offset_end}`;
        if (!groups.has(key)) {
            groups.set(key, {ontology_service: pm.ontology_service, offset_start: pm.offset_start,
                             offset_end: pm.offset_end, matches: []});
        }
        groups.get(key).matches.push(pm.match);
    }
    return Array.from(groups.values());
}

function phenotypeMatchGroupTitle(group) {
    if (group.matches.length === 1) {
        return group.matches[0];
    }
    return [`${group.matches.length} ${group.ontology_service} terms:`, ...group.matches].join("\n");
}

const PHENOTYPE_EXCLUDE_MARKER_TOOLTIP = "This string prevents the phenotype term from being officially matched and used. Remove it to signal human approval";

/* Each occurrence of the exclude string, as a range drawn alongside the matches */
function phenotypeExcludeMarkerRanges(phenotypeText, excludeString) {
    const ranges = [];
    if (excludeString) {
        let start = phenotypeText.indexOf(excludeString);
        while (start !== -1) {
            ranges.push({excludeMarker: true, offset_start: start, offset_end: start + excludeString.length});
            start = phenotypeText.indexOf(excludeString, start + excludeString.length);
        }
    }
    return ranges;
}

/* A warning whose detail (a list of terms) is too long to show until asked for */
function phenotypeWarningWithDetails(message, detailLabel, detailItems) {
    const listElement = $("<li/>").addClass("warning").text(message + " ");
    const details = $("<ul/>").addClass("ambiguous-candidates").hide();
    for (const item of detailItems) {
        $("<li/>").text(item).appendTo(details);
    }
    const toggle = $("<a/>", {href: "#"}).text(`Show ${detailLabel}`).click(function () {
        details.toggle();
        $(this).text(`${details.is(":visible") ? "Hide" : "Show"} ${detailLabel}`);
        return false;
    });
    return listElement.append(toggle, details);
}

function displayPhenotypeMatches(descriptionBox, phenotypeText, phenotypeMatches, excludeString) {
    const ambiguousAcronymCandidates = {};  // acronym -> [{accession, name}, ...]
    const ambiguous = {};
    const realMatches = [];
    for (const pm of phenotypeMatches) {
        const acronym = pm.ambiguous_alias;
        if (acronym) {
            if (!(acronym in ambiguousAcronymCandidates)) {
                ambiguousAcronymCandidates[acronym] = pm.ambiguous_alias_candidates || [];
            }
            continue;  // warning-only; not a real match - skip highlighting/grid
        }
        if (pm.ambiguous) {
            const key = `${pm.ontology_service}:${pm.ambiguous}`;
            if (!(key in ambiguous)) {
                ambiguous[key] = {text: pm.ambiguous, ontologyService: pm.ontology_service, accessions: new Set()};
            }
            ambiguous[key].accessions.add(pm.accession);
        }
        realMatches.push(pm);
    }

    // Longest first at the same start, so a match inside another opens inside it
    const excludeMarkers = phenotypeExcludeMarkerRanges(phenotypeText, excludeString);
    const matchGroups = groupPhenotypeMatchesByOntology(realMatches).concat(excludeMarkers).sort(
        (a, b) => (a.offset_start - b.offset_start) || (b.offset_end - a.offset_end));
    const openGroups = [];
    let nextGroup = 0;
    let phenotypeHTML = '';

    function closeGroupsEndingBy(offset) {
        for (let j = openGroups.length - 1; j >= 0; --j) {
            if (openGroups[j].offset_end <= offset) {
                phenotypeHTML += "</span>";
                openGroups.splice(j, 1);
            }
        }
    }

    for (let i = 0; i < phenotypeText.length; i++) {
        closeGroupsEndingBy(i);
        while (nextGroup < matchGroups.length && matchGroups[nextGroup].offset_start <= i) {
            const group = matchGroups[nextGroup++];
            if (group.excludeMarker) {
                phenotypeHTML += `<span class="phenotype-exclude-marker" title="${PHENOTYPE_EXCLUDE_MARKER_TOOLTIP}">`;
            } else {
                const serviceClass = group.ontology_service.toLowerCase();
                phenotypeHTML += `<span title="${escapeHtml(phenotypeMatchGroupTitle(group))}" class="ontology-service ${serviceClass}">`;
            }
            openGroups.push(group);
        }

        const char = phenotypeText.charAt(i);
        if (char === '\n') {
            phenotypeHTML += "<br />";
        } else {
            phenotypeHTML += escapeHtml(char);
        }
    }
    closeGroupsEndingBy(Infinity);

    descriptionBox.html(phenotypeHTML);

    const ambiguousAcronyms = Object.keys(ambiguousAcronymCandidates);
    if (Object.keys(ambiguous).length > 0 || ambiguousAcronyms.length > 0) {
        const phenoMessages = $("<div/>").addClass("phenotype-messages");
        const messageContainer = $("<ul/>").addClass("messages");
        phenoMessages.append(messageContainer);

        for (const {text, ontologyService, accessions} of Object.values(ambiguous)) {
            const msg = `'${text}' matched ${accessions.size} different ${ontologyService} terms, and genes from all of them ` +
                "are used - some may not apply. Please make it more specific (the full term name, or an ID).";
            messageContainer.append(phenotypeWarningWithDetails(msg, `${accessions.size} terms`, accessions));
        }
        for (const acronym of ambiguousAcronyms) {
            const candidates = ambiguousAcronymCandidates[acronym];
            const msg = `'${acronym}' is an ambiguous acronym (it matches multiple distinct ontology concepts) and has been excluded from gene-list matching. Please type the full term name or an HPO/OMIM/MONDO ID.`;
            if (candidates && candidates.length) {
                const candidateItems = candidates.map((c) => `${c.accession} — ${c.name}`);
                messageContainer.append(phenotypeWarningWithDetails(msg, `${candidates.length} possible matches`, candidateItems));
            } else {
                messageContainer.append($("<li/>").addClass("warning").text(msg));
            }
        }
        const clearDiv = descriptionBox.siblings("div.clear");
        clearDiv.after(phenoMessages);
    }

}


function phenotypeMatchesToGridData(phenotypeMatches) {
    const gridData = [];
    const accessionSet = new Set(); // unique terms only
    for (let i=0 ; i<phenotypeMatches.length ; ++i) {
        const pm = phenotypeMatches[i];
        if (pm.ambiguous_alias) {
            continue;  // warning-only entry, not a real match
        }
        if (!accessionSet.has(pm.accession)) {
            const row = {
                'ontology_service': pm.ontology_service,
                'accession': pm.accession,
                'name': pm.name,
                'gene_symbols': pm.gene_symbols.join(', '),
            };
            gridData.push(row);
            accessionSet.add(pm.accession);
        }
    }
    return gridData;
}
