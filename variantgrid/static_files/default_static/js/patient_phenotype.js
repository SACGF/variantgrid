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
    var overlapping_matches = [];

    var phenoLen = phenotypeText.length;
    var phenotypeHTML = '';
    var phenoOffsetStart = 0;
    
    for(var i=0; i<phenoLen ;i++) {
        var char = phenotypeText.charAt(i);
        var j, pm;
        
        var sliceEnd = 0;
        for (j=phenoOffsetStart ; j<phenotypeMatches.length ; ++j) {
            pm = phenotypeMatches[j];
            if (pm.offset_start <= i) {
                sliceEnd = j + 1; 
                var termClass = TERM_CLASSES[pm.match_type];
                var termTitle = termClass + ": " + pm.match;
                phenotypeHTML += "<span accession='" + pm.accession + "' title='" + termTitle + "' class='term-match " + termClass + "'>";
            } else {
                break;
            }
        }
        if (sliceEnd) {
              var sliced = phenotypeMatches.slice(phenoOffsetStart, sliceEnd);
              overlapping_matches = overlapping_matches.concat(sliced);
              phenoOffsetStart = sliceEnd;
        }
        
        if (char == '\n') {
            phenotypeHTML += "<br />";
        } else {
            phenotypeHTML += char;
        }
        
        for (j=overlapping_matches.length - 1 ; j>= 0 ; --j) {
            pm = overlapping_matches[j];
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
    var gridData = [];
    var accessionSet = new Set(); // unique terms only
    for (var i=0 ; i<phenotypeMatches.length ; ++i) {
        var pm = phenotypeMatches[i];
        if (!accessionSet.has(pm.accession)) {
            var termType = TERM_CLASSES[pm.match_type];
            var row = {
                'term_type' : termType,
                'accession' : pm.accession,
                'name' : pm.name,
                'gene_symbols' : pm.gene_symbols.join(', '), 
            };
            gridData.push(row);
            accessionSet.add(pm.accession);
        }
    }
    return gridData;
}
