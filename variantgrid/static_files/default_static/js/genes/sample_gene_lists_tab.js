// @ts-check
// genes/templates/genes/sample_gene_lists_tab.html
// Note: This is a tab so can be re-loaded (ie don't use global let/const)
function populateSampleGeneListFromJSON(data) {
    const sampleGeneListContainer = $(".sample-gene-list-container[sample-gene-list-id=" + data.pk + "]");
    const details = $(".sample-gene-list-details", sampleGeneListContainer);
    const analysisTemplates = $(".sample-gene-list-analysis-templates", sampleGeneListContainer);
    const toggle = $("a.show-sample-gene-list-toggle", sampleGeneListContainer);

    if (data.visible) {
        toggle.hide();
        details.addClass("show");
        $("button.hide", details).show();
        $("button.un-hide", details).hide();
        analysisTemplates.show();
    } else {
        toggle.show();
        details.removeClass("show");
        $("button.hide", details).hide();
        $("button.un-hide", details).show();
        analysisTemplates.hide();
    }

    const myActiveIcon = $(".active-sample-gene-list", details);
    if (data.active) {
        $(".active-sample-gene-list").hide();  // hide all others
        myActiveIcon.show();
        $("button.make-active", details).hide();
    } else {
        myActiveIcon.hide();
        $("button.make-active", details).show();
    }
}

function drawInitialSampleGeneLists() {
    const sampleGeneListData = readJsonData("sample-gene-lists-tab-data").sample_gene_lists_data;
    for (let i=0 ; i<sampleGeneListData.length ; i++) {
        populateSampleGeneListFromJSON(sampleGeneListData[i]);
    }
    $("#sample-gene-lists").show();
}

function modifySampleGeneList(that, data) {
    const sampleGeneListContainer = $(that).parents(".sample-gene-list-container");
    const sampleGeneListId = sampleGeneListContainer.attr("sample-gene-list-id");

    $.ajax({
        type: "POST",
        data: data,
        url: Urls.api_sample_gene_list(sampleGeneListId),
        success: populateSampleGeneListFromJSON,
    });
}

$(document).ready(function() {
    drawInitialSampleGeneLists();

    $("button.make-active").click(function() {
        modifySampleGeneList(this, {active: true});
    });

    $("button.hide").click(function() {
        modifySampleGeneList(this, {visible: false});
    });

    $("button.un-hide").click(function() {
        modifySampleGeneList(this, {visible: true});
    });

    $("form#new-gene-list-form").ajaxForm({target: "#sample-gene-list-tab-container"});

});
