// @ts-check
// classification/templates/classification/discordance_report_action.html
/*
function notesVal() {
    let notes = $('#notes').val().trim();
    return notes.length > 0 ? notes : null;
}
function isConfirmed() {
    return $('#discordance-confirm').is(':checked');
}
 */

function isDiscordant(clinSigToBuckets) {
    const usedBuckets = {};
    let noBuckets = false;
    $('.clin-sig-change').each((index, elem) => {
        const bucket = clinSigToBuckets[$(elem).val()];
        if (bucket) {
            usedBuckets[bucket] = true;
        } else {
            noBuckets = true;
        }
    });
    $('#no-bucket').collapse(noBuckets ? 'show' : 'hide');
    const discordant = Object.keys(usedBuckets).length > 1;
    const resultingStatusDom = $('#resulting-status');
    resultingStatusDom.text(discordant ? "Continued Discordance" : "Pending Concordance");

    const resolveButton = $('#resolve-button');
    const resolution = $('#resolution');

    if (discordant) {
        resolution.val("discordant");
        resolveButton.text("Finish: Mark as Continued Discordance");
        resolveButton.removeClass("btn-primary");
        resolveButton.addClass("btn-danger");
        resultingStatusDom.addClass("text-danger");
        resultingStatusDom.removeClass("text-success");
    } else {
        resolution.val("concordant");
        resolveButton.text("Finish: Mark as Pending Concordance");
        resolveButton.addClass("btn-primary");
        resolveButton.removeClass("btn-danger");
        resultingStatusDom.removeClass("text-danger");
        resultingStatusDom.addClass("text-success");
    }
    // resolveButton.text(discordant ? "Mark as Continued Discordance" : "Mark as Pending Concordance")
    /*
    if (isConfirmed() && notesVal()) {
        resolveButton.removeClass('disabled');
    } else {
        resolveButton.addClass('disabled');
    }
     */
}
/*
function submitCheck() {
    let warnings = [];
    if (notesVal() == null) {
        warnings.push("Please provide text in the notes field.");
    }
    if (!isConfirmed()) {
        warnings.push("All changes must be marked as agreed upon.");
    }
    if (warnings.length > 0) {
        let warning = warnings.join("\n");
        window.alert(`Before submitting:\n${warning}`)
        return false;  // don't submit if checkbox isn't ticked
    }
}
 */
function showFormCheck() {
    const value = $("input[name='outcome']:checked").val();
    console.log(value);
    const allValues = {
        "#pending-changes": value == "agree",
        "#postponed-changes": value == "postpone"
    };
    Object.entries(allValues).forEach(entry => {
        const [key, value] = entry;
        if (value) {
            $(key).slideDown();
        } else {
            $(key).slideUp();
        }
    });
}

function initDiscordanceReportAction(clinSigToBuckets) {
    $('#resolve-form input, #resolve-form select').change(() => isDiscordant(clinSigToBuckets));
    isDiscordant(clinSigToBuckets);

    $("input[name='outcome']").change(showFormCheck);
    showFormCheck();
}
