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
        // Only one list is active - the others lose the tick and can be made active again
        $(".active-sample-gene-list").hide();
        $("button.make-active").show();
        myActiveIcon.show();
        $("button.make-active", details).hide();
    } else {
        myActiveIcon.hide();
        $("button.make-active", details).show();
    }
}

// The DataTables ajax `data` hook of a genes table
function sampleGeneListGenesFilter(data, settings) {
    const genes = $("#" + settings.sTableId).closest(".sample-gene-list-genes");
    data["attention_only"] = $("input.attention-only", genes).is(":checked");
}

// A table is only created when its genes are first shown, so the tab doesn't wait on every list's
function initSampleGeneListGenes() {
    $(".sample-gene-list-genes").one("show.bs.collapse", function() {
        const genes = $(this);
        const sampleGeneListId = genes.closest(".sample-gene-list-container").attr("sample-gene-list-id");
        genes.append($("<table>", {
            id: "sample-gene-list-genes-datatable-" + sampleGeneListId,
            class: "sticky-header",
            "data-datatable-url": genes.attr("data-genes-datatable-url"),
            "data-datatable-data": "sampleGeneListGenesFilter",
        }));
    });

    $(".sample-gene-list-genes input.attention-only").change(function() {
        $(this).closest(".sample-gene-list-genes").find("table").DataTable().ajax.reload();
    });
}

function drawInitialSampleGeneLists(sampleGeneListData) {
    for (let i=0 ; i<sampleGeneListData.length ; i++) {
        const data = sampleGeneListData[i];
        populateSampleGeneListFromJSON(data);
        if (data.active && data.visible) {
            $(".sample-gene-list-container[sample-gene-list-id=" + data.pk + "] .sample-gene-list-genes").collapse("show");
        }
    }
    $("#sample-gene-list-items").show();
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

function initSampleGeneListsTab(sampleGeneListData) {
    initSampleGeneListGenes();
    drawInitialSampleGeneLists(sampleGeneListData);

    $("button.make-active").click(function() {
        modifySampleGeneList(this, {active: true});
    });

    $("button.hide").click(function() {
        modifySampleGeneList(this, {visible: false});
    });

    $("button.un-hide").click(function() {
        modifySampleGeneList(this, {visible: true});
    });

    $("form.sample-gene-list-form").ajaxForm({target: "#sample-gene-list-tab-container"});

    // The whole coverage table is on the sample page's QC tab
    $("a.sample-qc-tab-link").click(function(event) {
        event.preventDefault();
        $("a.nav-link[href='#sample-qc']").tab("show");
    });
}
