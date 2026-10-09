// @ts-check
// variantopedia/templates/variantopedia/tag_stats.html
/* global loadTagStatsCard, renderTagStatsGenes, renderTagStatsReTagged, renderTagStatsTagGenesOverTime, renderTagStatsHeadline, renderTagStatsOverTime, renderTagStatsUser, renderTagStatsByLab, renderTagStatsCoOccurrence */ // tag_stats.js
function taggedVariantsUrl(tagIds) {
    const params = tagIds.map((t) => "tag=" + encodeURIComponent(t)).join("&");
    return Urls.genome_build_variant_tags(readJsonData("tag-stats-data").genome_build_name) + "?" + params;
}

function selectedValues(selector) {
    return $(selector).val() || [];
}

function alleleOriginParam() {
    return "allele_origin=" + encodeURIComponent($("input[name=allele_origin]:checked").val());
}

function loadGenesCard() {
    const params = $.param({
        gene_symbols: selectedValues("#id_genes-gene_symbols"),
        tags: selectedValues("#id_genes-tags"),
    }, true);
    loadTagStatsCard("tag-stats-genes",
                     Urls.tag_stats_genes() + "?" + alleleOriginParam() + "&" + params,
                     renderTagStatsGenes);
}

function loadReTaggedCard() {
    const tag = $("#id_re-tagged-tag").val() || "";
    loadTagStatsCard("tag-stats-re-tagged",
                     Urls.tag_stats_re_tagged() + "?" + alleleOriginParam() + "&tag=" + encodeURIComponent(tag),
                     (data, $content) => renderTagStatsReTagged(data, $content, Urls.view_allele));
}

function loadTagGenesOverTimeCard() {
    const tag = $("#id_gene-time-tag").val() || "";
    loadTagStatsCard("tag-stats-tag-genes-over-time",
                     Urls.tag_stats_tag_genes_over_time() + "?" + alleleOriginParam() + "&tag=" + encodeURIComponent(tag),
                     renderTagStatsTagGenesOverTime);
}

function showSelectedTagsVariants() {
    const tagIds = selectedValues("#id_co-occurrence-tags");
    if (tagIds.length) {
        window.location = taggedVariantsUrl(tagIds);
    }
}

function loadAllTagStatsCards() {
    const alleleOrigin = "?" + alleleOriginParam();
    loadTagStatsCard("tag-stats-headline", Urls.tag_stats_headline() + alleleOrigin,
                     renderTagStatsHeadline);
    loadTagStatsCard("tag-stats-over-time", Urls.tag_stats_over_time() + alleleOrigin,
                     renderTagStatsOverTime);
    loadTagStatsCard("tag-stats-your-tagging", Urls.tag_stats_for_user() + alleleOrigin,
                     renderTagStatsUser);
    loadTagStatsCard("tag-stats-by-lab", Urls.tag_stats_by_lab() + alleleOrigin,
                     renderTagStatsByLab);
    loadTagStatsCard("tag-stats-co-occurrence", Urls.tag_stats_co_occurrence() + alleleOrigin,
                     (data, $content) => renderTagStatsCoOccurrence(data, $content, taggedVariantsUrl));
    loadGenesCard();
    loadReTaggedCard();
    loadTagGenesOverTimeCard();
}

$(document).ready(() => {
    loadAllTagStatsCards();

    $("input[name=allele_origin]").change(() => {
        // Keep the selection in the URL so a refresh or a shared link stays on this origin
        history.replaceState(null, "", "?" + alleleOriginParam());
        loadAllTagStatsCards();
    });
    $("#genes-recalculate").click(loadGenesCard);
    $("#re-tagged-recalculate").click(loadReTaggedCard);
    $("#gene-time-recalculate").click(loadTagGenesOverTimeCard);
    $("#co-occurrence-show-variants").click(showSelectedTagsVariants);
});
