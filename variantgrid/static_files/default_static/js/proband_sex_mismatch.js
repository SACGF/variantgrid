/* Trio/Quad wizards - the proband is picked client side, so compare the chosen sample's patient
   record sex against the sex detected from its chrX genotypes here, and let the scientist decide.
   Wants a '#proband-sex-mismatch' container holding '#proband-sex-message' and '#id_proband_sex',
   and one select per sample inside '.sample-select', in the same order as sampleSexes. */
function setupProbandSexMismatch(sampleSexes) {
    function showProbandSexMismatch() {
        const probandSexSelect = $("#id_proband_sex");
        let sexes = null;
        $("select", ".sample-select").each(function(i) {
            if ($(this).val() == 'P') {
                sexes = sampleSexes[i];
            }
        });
        const differ = sexes && sexes["patient_sex"] != 'U' && sexes["detected_sex"] != 'U'
                       && sexes["patient_sex"] != sexes["detected_sex"];
        if (differ) {
            $("#proband-sex-message").html("The patient record says <b>" + sexes["patient_sex_label"]
                + "</b> but the sample's chrX genotypes detected <b>" + sexes["detected_sex_label"]
                + "</b>. Choose which to use:");
            if (!probandSexSelect.val()) {
                probandSexSelect.val(sexes["detected_sex"]);
            }
        } else {
            probandSexSelect.val('');
        }
        $("#proband-sex-mismatch").toggle(Boolean(differ));
    }

    showProbandSexMismatch();
    $("select", ".sample-select").change(showProbandSexMismatch);
}
