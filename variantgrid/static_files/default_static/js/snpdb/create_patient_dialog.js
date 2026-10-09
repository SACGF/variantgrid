// @ts-check
// snpdb/templates/snpdb/create_patient_dialog.html
/* global setupPatientSamplenode */ // patient_samplenode.js
function setupCreatePatientDialog() {
    const form = $("#patient-form");
    const dialog = $("#new-patient-dialog");
    $("#create-patient-button", dialog).click(function() {
        form.submit();
    });

    function showFormErrors(errors) {
        const errorMessages = $(".error-messages", form).empty();
        for (const [field, messages] of Object.entries(errors)) {
            const label = $("label[for='id_" + field + "']", form).text().trim();
            const prefix = label ? label + ": " : "";
            $("<div/>").addClass("alert alert-danger")
                       .text(prefix + messages.join(", "))
                       .appendTo(errorMessages);
        }
    }

    setupPatientSamplenode({name: ''});

    form.ajaxForm({
        success: function(data) {
            if (data.error) {
                showFormErrors(data.error);
                return;
            }
            const patientSelect = $(window.activePatientSelect);
            setAutocompleteValue(patientSelect, data.patient_id, data.__str__);
            dialog.modal("hide");
        },
        error: function(data) {
            $(".error-messages", form).empty().append(
                $("<div/>").addClass("alert alert-danger").text("Error creating patient - please try again."));
        },
    });

    return dialog;
}

function addCreatePatientButton(dialog, patientSelect, sampleName, buttonContainer) {
    const form = dialog.find("form");
    const createPatientButton = $("<div/>")
        .addClass("click-to-add-button d-inline-block align-middle")
        .attr("title", "Create new patient...");

    createPatientButton.click(function() {
        form.resetForm().change();
        $(".error-messages", form).empty();
        $(".set-from-sample", form)
            .text("Set from sample ('" + sampleName + "')")
            .off("click")
            .click(function() {
                $("#id_" + $(this).data("field"), form).val(sampleName).change();
                return false;
            })
            .show();

        window.activePatientSelect = patientSelect;
        dialog.modal("show");
    });

    buttonContainer.append(createPatientButton);
}
