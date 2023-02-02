function getEvidenceHover(geneSymbol, evidence, confidence) {
    let container = $("<div />");
    let geneTitle = $("<h3>" + geneSymbol + "</h3>").addClass(confidence);
    container.append(geneTitle);
    let modeOfInheritance = evidence["mode_of_inheritance"];
    if (modeOfInheritance) {
        mohDiv = $("<div />").append("<b>Mode of inheritance</b> " + modeOfInheritance);
        container.append(mohDiv);
    }
    let penetrance = evidence["penetrance"];
    if (penetrance) {
        penetranceDiv = $("<div />").append("<b>Penetrance</b> " + penetrance);
        container.append(penetranceDiv);
    }

    let evidences = evidence["evidence"];
    if (evidences) {
        container.append($("<h3 />").text("Sources"));
        evidencesUl = $("<ul />");
        for (let i=0 ; i<evidences.length ; ++i) {
            evidencesUl.append($("<li />").text(evidences[i]));
        }
        container.append(evidencesUl);
    }

    let phenotypes = evidence["phenotypes"];
    if (phenotypes) {
        container.append($("<h3 />").text("Phenotypes"));
        phenotypesUl = $("<ul />");
        for (let i=0 ; i<phenotypes.length ; ++i) {
            phenotypesUl.append($("<li />").text(phenotypes[i]));
        }
        container.append(phenotypesUl);
    }
    return container;
}

function getConfidenceCSS(confidence) {
    let cssClass = null;
    if (confidence) {
        const CONFIDENCE_CLASS = ["LowEvidence", "ModerateEvidence", "HighEvidence"];
        cssClass = CONFIDENCE_CLASS[parseInt(confidence) - 1];
    }
    return cssClass;
}


function addGeneEvidence(geneSymbol, evidence, geneContainer) {
    let evidenceDiv = $("<div />").addClass("evidence hover-detail");
    let confidence = evidence["confidence_level"];
    let evidenceHover = getEvidenceHover(geneSymbol, evidence, confidence);
    evidenceDiv.attr("data-content", evidenceHover.html());
    let cssClass = getConfidenceCSS(confidence);
    if (cssClass) {
        evidenceDiv.addClass(cssClass);
    }
    geneContainer.append(evidenceDiv);
}


function getDivFromPanelAppGeneEvidenceAPIResult(geneSymbol, panelAppEvidenceResultsList) {
    let newDiv = $("<div />");
    let summaryText = null;

    if (panelAppEvidenceResultsList.length) {
        let confidenceCount = {};
        for (let i=0 ; i<panelAppEvidenceResultsList.length ; ++i) {
            let evidence = panelAppEvidenceResultsList[i];
            addGeneEvidence(geneSymbol, evidence, newDiv);

            let panel = evidence["panel"];
            if (panel) {
                let diseaseName = panel["name"];
                let confidenceLabel = getConfidenceCSS(evidence["confidence_level"]);
                let existingCount = confidenceCount[confidenceLabel] || 0;
                confidenceCount[confidenceLabel] = existingCount + 1;
                let diseaseSpan = $("<span />").addClass("disease-name").text(diseaseName + ". ");
                newDiv.append(diseaseSpan);
            }
        }

        if (confidenceCount) {
            let summaryDiv = $("<div />").addClass("hidden gene-symbol-summary");
            // go high to low
            let evidenceList = [];
            for (let evidence of ["HighEvidence", "ModerateEvidence", "LowEvidence"]) {
                let count = confidenceCount[evidence];
                if (count) {
                    evidenceList.push(evidence + ": " + confidenceCount[evidence]);
                }
            }
            summaryText = evidenceList.join(", ");
            newDiv.append(summaryDiv);
        }
    } else {
        newDiv.append("No evidence.");
    }
    return {
        div: newDiv,
        summary: summaryText,
    };
}


function getPanelAppGeneEvidenceDiv(serverName, server_id, geneSymbol, summaryFunc) {
    // Returns a Div which is loading, then is filled in when API retrieved ok
    let panelAppDiv = $("<div />").addClass("loading icon16");
    $.ajax({
        type: "GET",
        url: Urls.api_panel_app_gene_evidence(server_id, geneSymbol),
        success: function(panelAppEvidenceResultsList) {
            const {div, summary} = getDivFromPanelAppGeneEvidenceAPIResult(geneSymbol, panelAppEvidenceResultsList);
            panelAppDiv.replaceWith(div);
            if (summaryFunc) {
                let summaryText = serverName + " - " + summary;
                summaryFunc(summaryText);
            }
        },
        error: function(qXHR, textStatus, errorThrown) {
            let errorText = "Failed...";
            if (qXHR.status == 404) {
                errorText = "Not found";
            }
            let newDiv = $("<div />").text(errorText);
            panelAppDiv.replaceWith(newDiv);
        }
    });
    return panelAppDiv;
}