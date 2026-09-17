/*
 * Patient phenotype terms (issue #1426)
 *
 * A page loads every patient's terms in one go as patientPhenotypes ({patient_id: {text, terms}}),
 * then draws cells from it. Changing a patient select fetches that one patient from the endpoint
 * and caches it, so the column follows the select before the page is saved.
 */

const PHENOTYPE_ONTOLOGY_CSS_CLASSES = {HPO: "hpo", OMIM: "omim", MONDO: "mondo"};

function patientPhenotypeHtml(phenotype) {
    if (!phenotype) {
        return "";
    }
    const chips = [];
    for (const [service, cssClass] of Object.entries(PHENOTYPE_ONTOLOGY_CSS_CLASSES)) {
        for (const term of (phenotype.terms[service] || [])) {
            chips.push(`<span class="${cssClass} phenotype-chip" title="${escapeHtml(term.id)}">${escapeHtml(term.name)}</span>`);
        }
    }
    if (!chips.length) {
        return "";
    }
    return `<span class="patient-phenotype" title="${escapeHtml(phenotype.text || "")}">${chips.join("")}</span>`;
}

/* Each patient's text once, skipping anything the existing text already holds */
function patientPhenotypeTextsToAdd(existingText, patientIds, patientPhenotypes) {
    const texts = [];
    for (const patientId of new Set(patientIds.filter(Boolean).map(String))) {
        const phenotype = patientPhenotypes[patientId];
        const patientText = phenotype && phenotype.text ? phenotype.text.trim() : "";
        if (patientText && ![existingText, ...texts].join("\n").includes(patientText)) {
            texts.push(patientText);
        }
    }
    return texts;
}

/*
 * The "Add all sample patient phenotypes" button seeds a phenotype textarea from the page's patients, one
 * text per line. Nothing is saved: the keyup runs the editor's live matching so the user reviews and trims
 * before submitting the form. getPatientIds is read on every refresh, as membership and patient selects
 * change before saving - returns the refresh for callers to run when they do.
 */
function setupPatientPhenotypeSeedButton(textarea, getPatientIds, patientPhenotypes) {
    const button = $("#add-sample-patient-phenotypes");
    const summary = $("#sample-patient-phenotypes-summary");

    function refresh() {
        const patientIds = getPatientIds();
        const toAdd = patientPhenotypeTextsToAdd(textarea.val() || "", patientIds, patientPhenotypes);
        let message;
        if (toAdd.length) {
            message = `${toAdd.length} phenotype${toAdd.length === 1 ? "" : "s"} will be copied in.`;
        } else if (patientPhenotypeTextsToAdd("", patientIds, patientPhenotypes).length) {
            message = "Every sample patient phenotype is already in the text.";
        } else {
            message = "No sample patients have phenotype text.";
        }
        button.prop("disabled", !toAdd.length);
        summary.text(message);
    }

    button.click(function () {
        const existingText = textarea.val() || "";
        const toAdd = patientPhenotypeTextsToAdd(existingText, getPatientIds(), patientPhenotypes);
        textarea.val([existingText, ...toAdd].filter(Boolean).join("\n"));
        textarea.trigger("keyup");
        refresh();
    });
    textarea.on("input change keyup", refresh);
    refresh();
    return refresh;
}

/* Terms for one patient, from the cache if we already have them */
function loadPatientPhenotype(patientPhenotypes, patientId, callback) {
    if (!patientId) {
        callback(null);
        return;
    }
    if (patientId in patientPhenotypes) {
        callback(patientPhenotypes[patientId]);
        return;
    }
    $.getJSON(Urls.patient_phenotype_terms(patientId), function (data) {
        patientPhenotypes[patientId] = data;
        callback(data);
    }).fail(function () {
        callback(null);
    });
}
