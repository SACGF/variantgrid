// @ts-check
// analysis/templates/analysis/family_wizard.html
const familyWizardData = readJsonData("family-wizard-data");
const patient_description_results = familyWizardData.patient_description_results;
const sample_sexes = familyWizardData.sample_sexes;

function setProbandPhenoText() {
    const pp = $("#proband-phenos");
    pp.empty();
    const displayedTerms = new Set();

    $("select", ".sample-select").each(function(i) {
        if ($(this).val() === 'P') {
            const results = patient_description_results[i][1];
            if (results) {
                $("#proband-text").html("Using proband phenotypes:");
                for (let m = 0; m < results.length; ++m) {
                    const result = results[m];
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
