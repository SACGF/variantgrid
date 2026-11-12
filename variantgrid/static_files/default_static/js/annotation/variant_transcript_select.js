// @ts-check
// annotation/templates/annotation/tags/variant_transcript_select.html - the transcript choice radios. The page's
// transcriptFunc(transcriptId, data) and geneFunc(geneId, geneSymbol) (either may be null) are called with the
// chosen transcript's annotation, once for the initial choice and again on each change.
// Function declarations only - the analysis editor loads the variant page's scripts into one document repeatedly.

function initVariantTranscriptSelect(options) {
    const {transcriptSelectJquery, variantTranscriptAnnotations, initialTranscriptId, transcriptFunc, geneFunc} = options;

    function setTranscriptInfo(transcriptId) {
        const data = variantTranscriptAnnotations[transcriptId];
        if (transcriptFunc) {
            transcriptFunc(transcriptId, data);
        }
        if (data && geneFunc) {
            geneFunc(data["gene_id"], data["gene_symbol"]);
        }
    }

    const transcriptRadio = $('input[type=radio]', '#variant-transcript-select-tag-container');
    transcriptRadio.change(function() {
        if ($(this).is(":checked")) {
            setTranscriptInfo($(this).val());
        }
    });
    transcriptRadio.appendTo($(transcriptSelectJquery));

    setTranscriptInfo(initialTranscriptId);
}
