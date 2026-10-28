// The variant page's per-transcript UTRAnnotator and SpliceVault rows (#2426). transcriptFunc hands over the
// selected transcript's annotation (VariantTranscriptSelections._get_transcript_data): applyValue fills the scalar
// fields by id, this draws the '&'-joined records and hides a group the transcript has no values for.
// Function declarations only - the analysis editor loads the variant page's scripts into one document repeatedly.

function _renderTranscriptRecords(container, records, recordHtml) {
    container.empty();
    if (!records || !records.length) {
        container.append($("<span>", {class: "no-value", text: "-"}));
        return;
    }
    for (const record of records) {
        container.append(recordHtml(record));
    }
}

function _utrAnnotatorEffectHtml(effect) {
    const parts = Object.entries(effect).map(([key, value]) =>
        `${escapeHtml(key)} <span class="text-monospace">${escapeHtml(value)}</span>`);
    return $("<div>", {class: "text-muted small mt-1", html: parts.join(" &middot; ")});
}

function _spliceVaultEventHtml(event) {
    return $("<div>", {class: "text-monospace", text: event});
}

function applyTranscriptPluginRows(data) {
    const hasValue = name => data != null && data[name] != null;
    _renderTranscriptRecords($("#utr-annotator-effects"), data && data.utr_annotator_effects, _utrAnnotatorEffectHtml);
    $("#utr-annotator-orfs-values").toggle(hasValue("utr_annotator_existing_uorfs"));
    $("#utr-annotator-orfs-none").toggle(!hasValue("utr_annotator_existing_uorfs"));
    $("#splicevault-site-values").toggle(hasValue("splicevault_site_pos"));
    $("#splicevault-site-none").toggle(!hasValue("splicevault_site_pos"));
    _renderTranscriptRecords($("#splicevault-events"), data && data.splicevault_top_events_list, _spliceVaultEventHtml);
}
