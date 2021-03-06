var TERM_CLASSES = {"H" : "hpo", "O" : "omim", "G" : "gene"};  

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
    
    descriptionBox.html(phenotypeHTML);
    $(".term-match", descriptionBox);
}


function phenotypeMatchesToJqGridData(phenotypeMatches) {
    const gridData = [];
    const accessionSet = new Set(); // unique terms only
    for (let i=0 ; i<phenotypeMatches.length ; ++i) {
        const pm = phenotypeMatches[i];
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
