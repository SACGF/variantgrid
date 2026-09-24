/* IGV desktop integration - links that drive the user's own IGV through its batch port. The port, and
   the path prefixes to rewrite to where the user has the BAMs mounted, come from the IGV Integration
   settings page (@see snpdb/models/models_user_settings.py:get_igv_data).

   A page shows the links by setting window.ANALYSIS_SETTINGS = {show_igv_links, igv_data} and defining
   getBams(), the BAM paths to load alongside the locus. The settings are read off the analysis window,
   so variant details loaded inside an analysis follow the analysis'. igvPortUrl takes everything it
   needs as arguments. */

function getIgvSettings() {
    return getAnalysisWindow().ANALYSIS_SETTINGS || {};
}

// replaceDict is {serverPrefix: userPrefix} - {"/data/": "Z:\\"}
function replaceFilePrefix(replaceDict, bamFiles) {
    return bamFiles.filter(Boolean).map(bamFile => {
        for (const [fromValue, toValue] of Object.entries(replaceDict || {})) {
            if (bamFile.startsWith(fromValue)) {
                return bamFile.replace(fromValue, toValue);
            }
        }
        return bamFile;
    });
}

// 'goto' the locus, or 'load' the BAMs at it when there are any
function igvPortUrl(igvData, locus, bamFiles) {
    const params = ["genome=" + igvData.genome];
    if (locus) {
        params.push("locus=" + locus);
    }
    let op = "goto";
    const files = replaceFilePrefix(igvData.replace_dict, bamFiles || []).join();
    if (files) {
        params.push("file=" + files);
        op = "load";
    }
    return `${igvData.base_url}/${op}?${params.join("&")}`;
}

let seenIgvError = false;

function openIgvLink(locus, bamFiles) {
    const igvData = getIgvSettings().igv_data;
    $.ajax({
        url: igvPortUrl(igvData, locus, bamFiles),
        error: function(jqXHR, textStatus, errorThrown) {
            if (!seenIgvError) {
                console.log(jqXHR, textStatus, errorThrown);
                let message = `<p>Could not connect to IGV - is it running and accepting connections on ${igvData.base_url}?`;
                message += `<p>See also <a target='_blank' href='${Urls.igv_integration()}'>IGV Integration</a>`;
                createModal("igv-error-dialog", "IGV", message);
                seenIgvError = true;
            }
        },
        suppressErrors: true,
    });
}

// getBamsFuncString names a page function returning the BAMs, called on click so it sees the page as it is then
function createIgvUrl(locus, getBamsFuncString) {
    if (getIgvSettings().show_igv_links) {
        const bams = getBamsFuncString ? `${getBamsFuncString}()` : '[]';
        return `javascript:openIgvLink("${locus}", ${bams})`;
    }
    return null;
}

function createIgvLink(locus, getBamsFuncString) {
    const igvUrl = createIgvUrl(locus, getBamsFuncString);
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
        $(this).html(createIgvLink($(this).data("locus"), 'getBams'));
    });
}
