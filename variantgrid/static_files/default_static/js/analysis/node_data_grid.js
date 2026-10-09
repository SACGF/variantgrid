// @ts-check
// analysis/templates/analysis/node_data/node_data_grid.html
/* global nodeId:writable, nodeAlignmentsDict:writable */
/* global registerNodeGridDownloadButton, nodeGridHasData, export_grid, setupNodeGrid, gridLoadError, loadNodeGridData */ // grid.js
/* global isHorizontalMode, resizeGrid, bottomPaneGridHidden, showGridLoadingOverlay, registerDeferredGridLoad */ // analysis.js
/* global GRID */ // analysis/analysis_editor_and_grid.js
/* global revealSelectedTab */ // analysis/templates/analysis/node_editors/grid_editor.html
/* global load_node_editor, on_error_function */ // analysis/templates/analysis/node_data/base_node_data.html
function selectVariant(checkbox, gridData) {
    const variantId = $(checkbox).attr("variant_id");
    const checked = $(checkbox).is(":checked");

    const data = 'variant_id=' + variantId + '&checked=' + checked;
    $.ajax({
        type: "POST",
        data: data,
        url: Urls.set_variant_selected(gridData.analysisId, gridData.nodeId),
        success: function() {
            const aWin = getAnalysisWindow();
            const variants = aWin.selectedVariants[gridData.nodeId] || {};
            aWin.selectedVariants[gridData.nodeId] = variants;
            if (checked) {
                variants[variantId] = 1;
                revealSelectedTab(gridData.nodeId, true); // Could need to turn on...
            } else {
                delete variants[variantId];
            }
            checkAndMarkDirtyNodes(aWin);
        }
    });
}


function gridComplete(gridData) {
    const unique_code = gridData.nodeId + "_" + gridData.nodeVersion;
    if ($("#" + unique_code, "#node-data-container").length === 0) {
        return;  // user navigated away; ignore stale callback
    }

    const aWin = getAnalysisWindow();
    const variants = aWin.selectedVariants[gridData.nodeId];
    if (variants) {
        $("input.variant-select").each(function() {
            const variantId = $(this).attr("variant_id");
            if (variantId in variants) {
                $(this).prop('checked', true);
            }
        });
    }

    $("input.variant-select").click(function() { selectVariant(this, gridData); });
    registerComponent(unique_code, GRID);

    // The placeholder is the download route that matters most for big nodes - the user may never load
    // the grid at all - so keep its links in sync with any export already running for this node
    registerNodeGridDownloadButton("#placeholder-export-csv-" + gridData.nodeId, gridData.analysisId, gridData.nodeId,
                                   unique_code, 'csv', false, "CSV");
    registerNodeGridDownloadButton("#placeholder-export-vcf-" + gridData.nodeId, gridData.analysisId, gridData.nodeId,
                                   unique_code, 'vcf', false, "VCF");

    // A deferred grid draws once with no rows (so the editor's everythingLoaded can proceed). Only a
    // genuine row fetch counts as "loaded this session" - otherwise revisiting would auto-fire the
    // query the user never asked for.
    if (nodeGridHasData(gridData.nodeId, unique_code)) {
        // Record this node-version as loaded so a later revisit re-shows it automatically (page-cache hit).
        aWin.loadedGridVersions = aWin.loadedGridVersions || {};
        aWin.loadedGridVersions[gridData.nodeId] = gridData.nodeVersion;

        // Clear the grid-only overlay we put up when the user clicked "Show grid" (no-op otherwise).
        hideGridLoadingOverlay();

        if (typeof isHorizontalMode === "function" && isHorizontalMode()) {
            resizeGrid();  // fill the full width bottom panel the grid landed in
        }
    }
}


// Called by name from the IGV links - @see createIgvLink
function getAlignmentFiles() {
    const alignmentFiles = [];
    for (const k in nodeAlignmentsDict) {
        const sample_alignment_files = nodeAlignmentsDict[k];
        alignmentFiles.push(...sample_alignment_files);
    }
    return alignmentFiles;
}


/* CSV/VCF run as Celery jobs and hand back a file - the buttons show their progress, see
   analysis_downloads.js. They're built once the table is, so they can read its live ajax params. */
function addNodeGridExportButtons(unique_code, gridData) {
    const analysisId = gridData.analysisId;
    const gridNodeId = gridData.nodeId;
    const toolbar = $("#node-grid-toolbar-" + gridNodeId);
    const addButton = function(id, caption, title, exportType, useCanonicalTranscripts) {
        const link = $("<a>", {id: id, class: "btn btn-outline-secondary btn-sm", title: title,
                               href: "javascript:void(0)",
                               html: '<i class="fas fa-download"></i> ' + caption}).appendTo(toolbar);
        link.click(function() {
            export_grid(analysisId, gridNodeId, unique_code, exportType, useCanonicalTranscripts);
        });
        registerNodeGridDownloadButton("#" + id, analysisId, gridNodeId, unique_code,
                                       exportType, useCanonicalTranscripts, caption);
    };

    addButton("node-grid-export-csv-" + gridNodeId, "CSV", "Download as CSV", 'csv', false);

    const aWin = getAnalysisWindow();
    if (aWin.ANALYSIS_SETTINGS && aWin.ANALYSIS_SETTINGS.canonical_transcript_collection) {
        const ctc = aWin.ANALYSIS_SETTINGS.canonical_transcript_collection;
        addButton("node-grid-export-canonical-csv-" + gridNodeId, "Canonical transcript CSV",
                  "Download CSV using transcripts from " + ctc, 'csv', true);
    }

    addButton("node-grid-export-vcf-" + gridNodeId, "VCF", "Download as VCF", 'vcf', false);
}


/* Called by every node grid loaded. The grid callbacks are handed that load's values (gridData), so a callback
   that lands after the user has moved to another node still checks against its own node.
   gridData: analysisId, analysisVersion, nodeId, nodeVersion, extraFilters, nodeProbandSampleId,
             nodeProbandPatientId, alignmentsDict, gridAutoLoad */
function initNodeDataGrid(gridData) {
    nodeId = gridData.nodeId;
    // Tag pills are read against the node the grid is showing - a tagging for someone else, or for
    // nobody yet, is marked rather than looking like this proband's. @see VariantGridFormat.tags
    nodeProbandSampleId = gridData.nodeProbandSampleId;
    // A node above sample level is about a person without being about one of their VCFs
    nodeProbandPatientId = gridData.nodeProbandPatientId;
    nodeAlignmentsDict = gridData.alignmentsDict;
    $(document).ready(function() { nodeDataGridReady(gridData); });
}

function nodeDataGridReady(gridData) {
    const analysisId = gridData.analysisId;
    const gridNodeId = gridData.nodeId;
    const nodeVersion = gridData.nodeVersion;
    // Need to capture unique code and pass to functions as pages may be redefined by further DOM manipulation before load()s etc come back
    const unique_code = gridNodeId + "_" + nodeVersion;

    const node_view_url = Urls.node_view(analysisId, gridData.analysisVersion, gridNodeId, nodeVersion, gridData.extraFilters);
    load_node_editor(node_view_url, unique_code);

    const config_url = Urls.node_grid_config(analysisId, gridData.analysisVersion, gridNodeId, nodeVersion, gridData.extraFilters);
    const handler_url = Urls.node_grid_handler(analysisId);
    // Re-show instantly if this exact node-version was already loaded this session (page cache hit).
    const aWin = getAnalysisWindow();
    aWin.loadedGridVersions = aWin.loadedGridVersions || {};
    const alreadyLoaded = aWin.loadedGridVersions[gridNodeId] === nodeVersion;
    const wantsGrid = gridData.gridAutoLoad || alreadyLoaded;
    // Nothing to show while the Editor tab is up, so don't pay for the rows - the tab-show handler
    // in analysis.js fires the load we register below. A node over the threshold isn't deferred: its
    // placeholder shows with the tab, and the rows wait for "Show grid".
    const deferredForHiddenTab = wantsGrid && typeof bottomPaneGridHidden === "function" && bottomPaneGridHidden();
    const autoLoad = wantsGrid && !deferredForHiddenTab;

    const showGrid = function() {
        $("#grid-placeholder-" + gridNodeId).hide();
        $("#node-data-grid", "#" + unique_code).show();
    };

    if (autoLoad) {
        // We have retrieved it (or are about to) - show the grid, never the placeholder.
        showGrid();
    }
    const gridSetup = setupNodeGrid(config_url, handler_url, analysisId, gridNodeId, nodeVersion,
                                    unique_code, function() { gridComplete(gridData); }, gridLoadError, on_error_function, autoLoad);
    gridSetup.then(function(built) {
        if (built) {
            addNodeGridExportButtons(unique_code, gridData);
        }
    });

    const loadDeferredGrid = function() {
        showGrid();
        // Double-helix over just the grid section (editor stays usable) while phase 2 runs.
        showGridLoadingOverlay("#node-grid-container");
        // The Grid tab can come up before the table is built - the rows need the table to load into
        const loadOnceBuilt = function() {
            if ($("#" + unique_code, "#node-data-container").length === 0) {
                return;  // user moved to another node - the overlay is now that node's
            }
            if (!loadNodeGridData(gridNodeId, unique_code)) {  // fires deferred phase 2
                hideGridLoadingOverlay();  // no data callback is coming to clear it
            }
        };
        gridSetup.then(loadOnceBuilt, loadOnceBuilt);
    };

    if (deferredForHiddenTab) {
        registerDeferredGridLoad(loadDeferredGrid);
    }
    if (!autoLoad) {
        // Over the row-count threshold - the placeholder gives the user the choice
        $("#load-variants-" + gridNodeId).click(loadDeferredGrid);
    }
}
