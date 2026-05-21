var TERM_CLASSES = {"H" : "hpo", "O" : "omim", "G" : "gene"};  

function displayPhenotypeMatches(descriptionBox, phenotypeText, phenotypeMatches, excludeString) {
    function compareByStart(a,b) {
      if (a.offset_start < b.offset_start)
        return -1;
      if (a.offset_start > b.offset_start)
        return 1;
      return 0;
    }
    
    phenotypeMatches = phenotypeMatches.sort(compareByStart);
    let overlapping_matches = [];
    let ambiguousAcronymCandidates = {};  // acronym -> [{accession, name}, ...]
    const realMatches = [];
    for (let k = 0; k < phenotypeMatches.length; ++k) {
        const pm = phenotypeMatches[k];
        const acronym = pm.ambiguous_alias;
        if (acronym) {
            if (!(acronym in ambiguousAcronymCandidates)) {
                ambiguousAcronymCandidates[acronym] = pm.ambiguous_alias_candidates || [];
            }
            continue;  // warning-only; not a real match - skip highlighting/grid
        }
        realMatches.push(pm);
    }
    phenotypeMatches = realMatches;

    const phenoLen = phenotypeText.length;
    let phenotypeHTML = '';
    let phenoOffsetStart = 0;

    for(let i=0; i<phenoLen ; i++) {
        const char = phenotypeText.charAt(i);

        let sliceEnd = 0;
        for (let j=phenoOffsetStart ; j<phenotypeMatches.length ; ++j) {
            let pm = phenotypeMatches[j];
            if (pm.offset_start <= i) {
                sliceEnd = j + 1;
                const termClass = TERM_CLASSES[pm.match_type];
                phenotypeHTML += `<span title='${pm.match}' class='term-match ${termClass}'>`;
            } else {
                break;
            }
        }
        if (sliceEnd) {
            const sliced = phenotypeMatches.slice(phenoOffsetStart, sliceEnd);
            overlapping_matches = overlapping_matches.concat(sliced);
            phenoOffsetStart = sliceEnd;
        }
        
        if (char === '\n') {
            phenotypeHTML += "<br />";
        } else {
            phenotypeHTML += char;
        }
        
        for (let j=overlapping_matches.length - 1 ; j>= 0 ; --j) {
            let pm = overlapping_matches[j];
            if (pm.offset_end <= i) {
                phenotypeHTML += "</span>";
                overlapping_matches.splice(j, 1);
            }
        }
    }
    
    if (excludeString) {
        const escapedRegex = excludeString.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
        const escapedText = excludeString.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
        const tooltip = "This string prevents the phenotype term from being officially matched and used. Remove it to signal human approval";
        phenotypeHTML = phenotypeHTML.replace(new RegExp(escapedRegex, 'g'),
            `<span class="phenotype-exclude-marker" title="${tooltip}">${escapedText}</span>`);
    }

    descriptionBox.html(phenotypeHTML);
    $(".term-match", descriptionBox);

    const ambiguousAcronyms = Object.keys(ambiguousAcronymCandidates);
    if (ambiguousAcronyms.length > 0) {
        let phenoMessages = $("<div/>").addClass("phenotype-messages");
        let messageContainer = $("<ul/>").addClass("messages");
        phenoMessages.append(messageContainer);
        for (const acronym of ambiguousAcronyms) {
            const candidates = ambiguousAcronymCandidates[acronym];
            let listElement = $("<li/>").addClass("warning");
            listElement.text(`'${acronym}' is an ambiguous acronym (it matches multiple distinct ontology concepts) and has been excluded from gene-list matching. Please type the full term name or an HPO/OMIM/MONDO ID.`);
            if (candidates && candidates.length) {
                let intro = $("<div/>").text("Possible matches:");
                let candList = $("<ul/>").addClass("ambiguous-candidates");
                for (const c of candidates) {
                    $("<li/>").text(`${c.accession} — ${c.name}`).appendTo(candList);
                }
                listElement.append(intro);
                listElement.append(candList);
            }
            messageContainer.append(listElement);
        }
        let clearDiv = descriptionBox.siblings("div.clear");
        clearDiv.after(phenoMessages);
    }
}


function phenotypeMatchesToJqGridData(phenotypeMatches) {
    const gridData = [];
    const accessionSet = new Set(); // unique terms only
    for (let i=0 ; i<phenotypeMatches.length ; ++i) {
        const pm = phenotypeMatches[i];
        if (pm.ambiguous_alias) {
            continue;  // warning-only entry, not a real match
        }
        if (!accessionSet.has(pm.accession)) {
            const termType = TERM_CLASSES[pm.match_type];
            const row = {
                'term_type': termType,
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
