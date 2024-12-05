function displayPhenotypeMatches(descriptionBox, phenotypeText, phenotypeMatches) {
    function compareByStart(a,b) {
      if (a.offset_start < b.offset_start)
        return -1;
      if (a.offset_start > b.offset_start)
        return 1;
      return 0;
    }
    
    phenotypeMatches = phenotypeMatches.sort(compareByStart);
    let overlapping_matches = [];

    const phenoLen = phenotypeText.length;
    let phenotypeHTML = '';
    let phenoOffsetStart = 0;

    let ambiguous = {};
    function addToAmbiguous(key, value) {
        if (!(key in ambiguous)) {
            ambiguous[key] = new Set();
        }
        ambiguous[key].add(value);
    }


    for(let i=0; i<phenoLen ; i++) {
        const char = phenotypeText.charAt(i);

        let sliceEnd = 0;
        for (let j=phenoOffsetStart ; j<phenotypeMatches.length ; ++j) {
            let pm = phenotypeMatches[j];
            if (pm.ambiguous) {
                addToAmbiguous(pm.ambiguous, pm.accession)
            }
            if (pm.offset_start <= i) {
                sliceEnd = j + 1;
                phenotypeHTML += `<span title='${pm.match}' class='ontology-service ${pm.ontology_service}'>`;
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
    
    descriptionBox.html(phenotypeHTML);
    $(".term-match-ontology-service", descriptionBox);

    if (Object.keys(ambiguous).length > 0) {
        let phenoMessages = $("<div/>").addClass("phenotype-messages");
        let messageContainer = $("<ul/>").addClass("messages");
        phenoMessages.append(messageContainer);

        let msg;
        for (const [text, termSet] of Object.entries(ambiguous)) {
            const terms = Array.from(termSet).join(', ');
            msg = `Phenotype: ${text}' was ambiguous (matched >=2 times in the same ontology service): ${terms}. Please resolve by being more specific`;
            let listElement = $("<li/>").addClass("warning");
            listElement.text(msg);
            messageContainer.append(listElement);
        }
        let clearDiv = descriptionBox.siblings("div.clear")
        clearDiv.after(phenoMessages);
    }

}


function phenotypeMatchesToJqGridData(phenotypeMatches) {
    const gridData = [];
    const accessionSet = new Set(); // unique terms only
    for (let i=0 ; i<phenotypeMatches.length ; ++i) {
        const pm = phenotypeMatches[i];
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
