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

/*
 * Seed a phenotype textarea from the patients on the page - each patient's text once, on its own line,
 * skipping anything the textarea already holds. Nothing is saved: the keyup runs the editor's live
 * matching so the user reviews and trims before submitting the form.
 */
function addPatientPhenotypesToText(textarea, patientIds, patientPhenotypes) {
    let text = textarea.val() || "";
    const seen = new Set();
    for (const patientId of patientIds) {
        if (!patientId || seen.has(String(patientId))) {
            continue;
        }
        seen.add(String(patientId));
        const phenotype = patientPhenotypes[patientId];
        const patientText = phenotype && phenotype.text ? phenotype.text.trim() : "";
        if (!patientText || text.indexOf(patientText) !== -1) {
            continue;
        }
        text = text ? text + "\n" + patientText : patientText;
    }
    textarea.val(text);
    textarea.trigger("keyup");
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
