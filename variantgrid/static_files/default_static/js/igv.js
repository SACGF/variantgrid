/* IGV desktop integration - links that drive the user's own IGV through its batch port. The port, and
   the path prefixes to rewrite to where the user has the BAMs/CRAMs mounted, come from the IGV Integration
   settings page (@see snpdb/models/models_user_settings.py:get_igv_data).

   A page shows the links by setting window.ANALYSIS_SETTINGS = {show_igv_links, igv_data} and defining
   getAlignmentFiles(), the BAM/CRAM paths to load alongside the locus. The settings are read off the analysis window,
   so variant details loaded inside an analysis follow the analysis'. igvPortUrls takes everything it
   needs as arguments. */

function getIgvSettings() {
    return getAnalysisWindow().ANALYSIS_SETTINGS || {};
}

// replaceDict is {serverPrefix: userPrefix} - {"/data/": "Z:\\"}
function replaceFilePrefix(replaceDict, alignmentFiles) {
    return alignmentFiles.filter(Boolean).map(alignmentFile => {
        for (const [fromValue, toValue] of Object.entries(replaceDict || {})) {
            if (alignmentFile.startsWith(fromValue)) {
                return alignmentFile.replace(fromValue, toValue);
            }
        }
        return alignmentFile;
    });
}

// Port requests, sent in order. IGV's /goto reads only the locus, so without alignment files the genome is
// switched by a /load of its own first (a no-op when IGV already has it). A /load with files takes all three.
function igvPortUrls(igvData, locus, alignmentFiles) {
    const genomeParam = "genome=" + igvData.genome;
    const files = replaceFilePrefix(igvData.replace_dict, alignmentFiles || []).join();
    if (files) {
        const params = [genomeParam, "file=" + files];
        if (locus) {
            params.push("locus=" + locus);
        }
        return [`${igvData.base_url}/load?${params.join("&")}`];
    }
    const urls = [`${igvData.base_url}/load?${genomeParam}`];
    if (locus) {
        urls.push(`${igvData.base_url}/goto?locus=${locus}`);
    }
    return urls;
}

let seenIgvError = false;

function openIgvLink(locus, alignmentFiles) {
    const igvData = getIgvSettings().igv_data;
    const [firstUrl, ...laterUrls] = igvPortUrls(igvData, locus, alignmentFiles);
    let request = $.ajax({url: firstUrl, suppressErrors: true});
    for (const url of laterUrls) {
        request = request.then(() => $.ajax({url: url, suppressErrors: true}));
    }
    request.fail(function(jqXHR, textStatus, errorThrown) {
        if (!seenIgvError) {
            console.log(jqXHR, textStatus, errorThrown);
            let message = `<p>Could not connect to IGV - is it running and accepting connections on ${igvData.base_url}?`;
            message += `<p>See also <a target='_blank' href='${Urls.igv_integration()}'>IGV Integration</a>`;
            createModal("igv-error-dialog", "IGV", message);
            seenIgvError = true;
        }
    });
}

// The locus goes into a javascript: URL, which is percent-decoded before it runs - only contig:position(-end)
// or gene symbols, space separated for a fusion's partners (@see GeneFusion.igv_locus), can go in
const IGV_LOCUS_REGEX = /^[\w.-]+(:\d+(-\d+)?)?( [\w.-]+(:\d+(-\d+)?)?)*$/;

// getAlignmentFilesFuncString names a page function returning the BAMs/CRAMs, called on click so it sees the page as it is then
function createIgvUrl(locus, getAlignmentFilesFuncString) {
    if (getIgvSettings().show_igv_links && IGV_LOCUS_REGEX.test(locus)) {
        const alignmentFiles = getAlignmentFilesFuncString ? `${getAlignmentFilesFuncString}()` : '[]';
        return `javascript:openIgvLink("${locus}", ${alignmentFiles})`;
    }
    return null;
}

function createIgvLink(locus, getAlignmentFilesFuncString) {
    const igvUrl = createIgvUrl(locus, getAlignmentFilesFuncString);
    if (igvUrl) {
        return createGridLink("Open " + locus + " in IGV", igvUrl, '', [], ['igv-link']);
    }
    return '';
}

// Fills the server rendered IGV placeholders in a fragment - '<span class="igv-locus"
// data-locus="chrX:66905968-66914514"></span>'. The link is built here rather than in the template
// so it follows ANALYSIS_SETTINGS.show_igv_links wherever the fragment was dropped in
function renderIgvLocusLinks(container) {
    $(".igv-locus", container || document).each(function() {
        $(this).html(createIgvLink($(this).data("locus"), 'getAlignmentFiles'));
    });
}
