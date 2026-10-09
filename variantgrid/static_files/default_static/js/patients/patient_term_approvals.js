// @ts-check
// patients/templates/patients/patient_term_approvals.html
/* global displayPhenotypeMatches */ // patient_phenotype.js
function initPatientTermApprovals(patientResults) {
    for (const p in patientResults) {
        const patientSelector = $(".unapproved-patient[patient_id=" + p + "]");
        const resultSelector = $(".results", patientSelector);
        const phenotypeMatches = patientResults[p];
        const phenotypeText = $("textarea.phenotype", patientSelector).val();
        displayPhenotypeMatches(resultSelector, phenotypeText, phenotypeMatches);
    }

    function greyOut(selector) {
        // Grey out the div as per http://stackoverflow.com/a/14461824/4846002
        selector.fadeTo('slow', .6);
        selector.append('<div style="position: absolute;top:0;left:0;width: 100%;height:100%;z-index:2;opacity:0.4;filter: alpha(opacity = 50)"></div>');
        const stopPropFn = function (e) {
            e.stopPropagation();
            e.preventDefault();
        };
        selector.bind("keydown", stopPropFn).bind("keypress", stopPropFn).bind("paste", stopPropFn);
    }

    $("button.accept-phenotype").click(function() {
        const parentContainer = $(this).parents(".unapproved-patient");
        greyOut(parentContainer);

        const patientId = $(this).attr("patient_id");
        const data = 'patient_id=' + patientId;
        $.ajax({
            type: "POST",
            data: data,
            url: Urls.approve_patient_term(),
            success: function(phenotypeMatches) {
                parentContainer.fadeOut();
            },
        });
    });
}
