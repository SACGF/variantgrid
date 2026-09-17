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
            if (!(pm.ambiguous in ambiguous)) {
                ambiguous[pm.ambiguous] = new Set();
            }
            ambiguous[pm.ambiguous].add(pm.accession);
        }
        realMatches.push(pm);
    }

    // Longest first at the same start, so a match inside another opens inside it
    const matchGroups = groupPhenotypeMatchesByOntology(realMatches).sort(
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
            const serviceClass = group.ontology_service.toLowerCase();
            phenotypeHTML += `<span title="${escapeHtml(phenotypeMatchGroupTitle(group))}" class="ontology-service ${serviceClass}">`;
            openGroups.push(group);
        }

        const char = phenotypeText.charAt(i);
        if (char === '\n') {
            phenotypeHTML += "<br />";
        } else {
            phenotypeHTML += char;
        }
    }
    closeGroupsEndingBy(Infinity);

    if (excludeString) {
        const escapedRegex = excludeString.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
        const escapedText = excludeString.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
        const tooltip = "This string prevents the phenotype term from being officially matched and used. Remove it to signal human approval";
        phenotypeHTML = phenotypeHTML.replace(new RegExp(escapedRegex, 'g'),
            `<span class="phenotype-exclude-marker" title="${tooltip}">${escapedText}</span>`);
    }

    descriptionBox.html(phenotypeHTML);
    $(".term-match-ontology-service", descriptionBox);

    const ambiguousAcronyms = Object.keys(ambiguousAcronymCandidates);
    if (Object.keys(ambiguous).length > 0 || ambiguousAcronyms.length > 0) {
        const phenoMessages = $("<div/>").addClass("phenotype-messages");
        const messageContainer = $("<ul/>").addClass("messages");
        phenoMessages.append(messageContainer);

        for (const [text, termSet] of Object.entries(ambiguous)) {
            const terms = Array.from(termSet).join(', ');
            const msg = `Phenotype: ${text}' was ambiguous (matched >=2 times in the same ontology service): ${terms}. Please resolve by being more specific`;
            const listElement = $("<li/>").addClass("warning");
            listElement.text(msg);
            messageContainer.append(listElement);
        }
        for (const acronym of ambiguousAcronyms) {
            const candidates = ambiguousAcronymCandidates[acronym];
            const listElement = $("<li/>").addClass("warning");
            listElement.text(`'${acronym}' is an ambiguous acronym (it matches multiple distinct ontology concepts) and has been excluded from gene-list matching. Please type the full term name or an HPO/OMIM/MONDO ID.`);
            if (candidates && candidates.length) {
                const intro = $("<div/>").text("Possible matches:");
                const candList = $("<ul/>").addClass("ambiguous-candidates");
                for (const c of candidates) {
                    $("<li/>").text(`${c.accession} — ${c.name}`).appendTo(candList);
                }
                listElement.append(intro);
                listElement.append(candList);
            }
            messageContainer.append(listElement);
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
