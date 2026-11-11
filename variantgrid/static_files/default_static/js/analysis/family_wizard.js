// @ts-check
// analysis/templates/analysis/family_wizard.html
/* global displayPhenotypeMatches */ // patient_phenotype.js
/* global setupProbandSexMismatch */ // proband_sex_mismatch.js
/* global setupFamilyRoles */ // family_wizard_roles.js

/** The proband's HPO/OMIM matches, listed in the "Creating..." dialog while the form submits */
function setProbandPhenoText(patientDescriptionResults) {
    const pp = $("#proband-phenos");
    pp.empty();
    const displayedTerms = new Set();

    $("select", ".sample-select").each(function(i) {
        if ($(this).val() === 'P') {
            const results = patientDescriptionResults[i][1];
            if (results) {
                $("#proband-text").html("Using proband phenotypes:");
                for (const result of results) {
                    const ontologyService = result["ontology_service"];
                    const resultMatch = result["match"];
                    if (!displayedTerms.has(resultMatch) && (ontologyService === 'HPO' || ontologyService === 'OMIM')) {
                        pp.append($("<li>").text(resultMatch));
                        displayedTerms.add(resultMatch);
                    }
                }
            }
        }
    });
}

/**
 * @param {{patientDescriptionResults: Array, sampleSexes: Array, roleOptions: Object}} options
 */
function initFamilyWizard(options) {
    const patientDescriptionResults = options.patientDescriptionResults;

    patientDescriptionResults.forEach(function(pdr, i) {
        const phenotypeText = pdr[0];
        const phenotypeMatches = pdr[1];
        if (phenotypeText) {
            displayPhenotypeMatches($("#patient-results-" + (i + 1)), phenotypeText, phenotypeMatches);
        }
    });

    setupProbandSexMismatch(options.sampleSexes);

    $(".family-name-suggestion").click(function(e) {
        e.preventDefault();
        $("#id_name").val($(this).data("name"));
    });

    $("#create-family-button").click(function() {
        setProbandPhenoText(patientDescriptionResults);
        $("#form-submitted-dialog").modal("show");
    });

    setupFamilyRoles(options.sampleSexes, options.roleOptions);
}
