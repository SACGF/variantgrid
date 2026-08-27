function startDualScreenMode() {
    const STAND_ALONE_URL = Urls.standalone_analysis_editor_and_grid(ANALYSIS_ID);
    if (!secondWindow) {
        saveSettingsOnResize = false;
        secondWindow = window.open(STAND_ALONE_URL, '_blank', 'toolbar=0,location=0,menubar=0');
        closeNodeEditorDrawer();  // the editor lives in the other window now
        $("#right-panel").empty();
        $('div#analysis-outer-container').removePane("east");
        $("#dual-screen-button").hide();
    } else {
        console.log("Already in dual screen mode!");
    }
}

function getGridAndEditorWindow() {
    return secondWindow || window;
}

function secondWindowClosing() {
    if (secondWindow) {
        $('div#analysis-outer-container').addPane("east");
        loadInitialGridEditor();
        secondWindow = null;
        $("#dual-screen-button").show();
    }
}

function loadNodeData(nodeId, extra_filters, fromSelectNode) {
    const win = getGridAndEditorWindow();
    win.loadGridAndEditorForNode(nodeId, extra_filters, fromSelectNode);
}

function isHorizontalMode() {
    return typeof ANALYSIS_HORIZONTAL_MODE !== "undefined" && ANALYSIS_HORIZONTAL_MODE;
}

// The height available to the grid inside the bottom panel - what's left under whatever sits above it
function gridPaneAvailableHeight() {
    const panel = $("#right-panel");
    const gridContainer = $("#node-grid-container");
    if (!panel.length || !gridContainer.length) {
        return 0;
    }
    const offsetInPanel = gridContainer.offset().top - panel.offset().top + panel.scrollTop();
    return panel.innerHeight() - offsetInPanel;
}

const MIN_GRID_HEIGHT = 100;
const GRID_BOTTOM_MARGIN = 10;

function resizeGrid() {
    const wrapper = $("#node-grid-container .dataTables_wrapper:visible");
    if (!wrapper.length) {
        return;
    }
    // The grid pane is the full window width in horizontal mode, so give the rows the pane height too -
    // otherwise the scroll body keeps them in a short box with the panel empty underneath
    const available = isHorizontalMode() ? gridPaneAvailableHeight() : null;

    wrapper.each(function() {
        const scrollBody = $(".dataTables_scrollBody", this);
        if (available && scrollBody.length) {
            // Everything but the rows - toolbar, column headers, pager
            const chrome = $(this).outerHeight(true) - scrollBody.height();
            scrollBody.css("max-height", Math.max(available - chrome - GRID_BOTTOM_MARGIN, MIN_GRID_HEIGHT) + "px");
        }
        const table = $("table.dataTable", this).first();
        if (table.length && $.fn.DataTable.isDataTable(table)) {
            // The pane changed width - re-sync the cloned scroll header with the body
            table.DataTable().columns.adjust();
        }
    });
}

function savePanelWidthSettings() {
    if (saveSettingsOnResize) {
        // analysis_panel_fraction is the fraction taken by the node canvas in both modes
        const canvasPanel = $("#left-panel");
        let analysis_panel_fraction;
        if (isHorizontalMode()) {
            analysis_panel_fraction = canvasPanel.outerHeight() / $("#analysis-and-toolbar-container").height();
        } else {
            analysis_panel_fraction = canvasPanel.outerWidth() / window.innerWidth;
        }
        const data = 'analysis_panel_fraction=' + analysis_panel_fraction;
        $.ajax({
            type: "POST",
            data: data,
            url: Urls.analysis_set_panel_size(ANALYSIS_ID),
        });
    }
}

function resizePanel() {
    clearTimeout(panelResizeTimeout);
    panelResizeTimeout = setTimeout(savePanelWidthSettings, panelResizeUpdateDelay);
}

/* Horizontal mode gives the node editor two homes: docked as a tab in the bottom panel beside the
   grid ("Swap"), or slid in over the canvas ("Drawer"). It's one editor either way - #node-editor-container
   is reparented between them, so everything that loads into it by id keeps working. The choice lives in
   localStorage while both are being compared. */
const BOTTOM_PANE_EDITOR = "editor";
const BOTTOM_PANE_GRID = "grid";
const BOTTOM_PANE_VARIANT = "variant";
const EDITOR_HOME_TAB = "tab";
const EDITOR_HOME_DRAWER = "drawer";
const EDITOR_HOME_STORAGE_KEY = "analysisEditorHome";
const DRAWER_MIN_WIDTH = 300;

let activeBottomPaneTab = BOTTOM_PANE_GRID;
// Which of the node editor's own tabs the strip is pointing at - 0 is the editor/grid pair
const NODE_EDITOR_TAB_EDITOR = 0;
let activeNodeEditorTab = NODE_EDITOR_TAB_EDITOR;
let activeVariantTab = null;  // variant id of the details tab that's up, if any
let paneBeforeVariantTab = BOTTOM_PANE_GRID;  // what to go back to when we leave the variant tabs

function getEditorHome() {
    if (!isHorizontalMode()) {
        return EDITOR_HOME_TAB;
    }
    return localStorage.getItem(EDITOR_HOME_STORAGE_KEY) === EDITOR_HOME_DRAWER ? EDITOR_HOME_DRAWER : EDITOR_HOME_TAB;
}

function bottomPaneTabsEnabled() {
    return isHorizontalMode() && getEditorHome() === EDITOR_HOME_TAB && $("#bottom-pane-tabs").length > 0;
}

// The grid tab defers its data fetch while hidden - editing a chain of nodes then costs no grid queries
function bottomPaneGridHidden() {
    return bottomPaneTabsEnabled() && activeBottomPaneTab !== BOTTOM_PANE_GRID;
}

// Set by node_data_grid.html when it built a grid whose rows we skipped because the Grid tab was hidden
let deferredGridLoad = null;
// A big grid still sits behind its placeholder until the user picks the Grid tab themselves
let deferredGridLoadNeedsUserRequest = false;

function registerDeferredGridLoad(loadFunc, needsUserRequest) {
    deferredGridLoad = loadFunc;
    deferredGridLoadNeedsUserRequest = Boolean(needsUserRequest);
}

// Editor tabs share a pane, so they're told apart by which node editor tab they point at
function updateMirroredTabHighlight() {
    const variantShowing = activeBottomPaneTab === BOTTOM_PANE_VARIANT;
    $(".bottom-pane-tab, .drawer-tab").each(function() {
        const tab = $(this);
        const pane = tab.attr("pane");
        let active;
        if (pane === BOTTOM_PANE_GRID) {
            active = activeBottomPaneTab === BOTTOM_PANE_GRID;
        } else if (pane === BOTTOM_PANE_VARIANT) {
            active = variantShowing && parseInt(tab.attr("variant_id")) === activeVariantTab;
        } else {
            // The drawer only ever shows the editor or a variant, so there's no pane to take into account
            const editorShowing = !variantShowing &&
                (!bottomPaneTabsEnabled() || activeBottomPaneTab === BOTTOM_PANE_EDITOR);
            active = editorShowing && parseInt(tab.attr("editor-tab")) === activeNodeEditorTab;
        }
        tab.toggleClass("active", active);
    });
}

/* The node editor's own tab strip is mirrored into whichever bar is always on screen, so horizontal
   mode only ever has one row of tabs: the bottom pane strip while docked (Editor | Grid | Summary...),
   or the drawer header while undocked. Rebuilt by each editor load - see base_editor.html */
function syncNodeEditorTabs() {
    const docked = getEditorHome() === EDITOR_HOME_TAB;
    $("#bottom-pane-node-tabs, #drawer-node-tabs").empty();
    const mirror = docked ? $("#bottom-pane-node-tabs") : $("#drawer-node-tabs");
    if (!mirror.length) {
        return;
    }
    activeNodeEditorTab = NODE_EDITOR_TAB_EDITOR;
    // Docked, the strip's own Editor and Grid tabs already cover the editor's first tab
    const firstTab = docked ? NODE_EDITOR_TAB_EDITOR + 1 : NODE_EDITOR_TAB_EDITOR;
    $("#node-editor-tabs > ul > li > a").slice(firstTab).each(function(i) {
        const editorTab = firstTab + i;
        $("<a>", {
            "class": docked ? "bottom-pane-tab" : "drawer-tab",
            "pane": BOTTOM_PANE_EDITOR,
            "editor-tab": editorTab,
            href: "javascript:void(0)",
            text: $(this).text(),
        }).appendTo(mirror).click(function() {
            showNodeEditorTab(editorTab);
        });
    });
    // The tabs name what's showing, so the drawer's own title gives up its space to them
    $("#node-editor-drawer").toggleClass("has-node-tabs", !docked && mirror.children().length > 0);
    updateMirroredTabHighlight();
}

/* Variant details from the grid open as their own closable tab rather than taking over the editor,
   so you can keep a few variants beside the node you're working on. They only ever come from the grid
   that's showing, and go when it does - see closeAllVariantDetailsTabs(). The details page uses page level
   ids, so only the tab being read is in the DOM - the others are remembered by url and re-loaded when
   picked, which is all a first open does anyway. Tabs are named by the page itself once it loads,
   see setVariantDetailsTabLabel() */
let openVariantTabs = [];  // {variantId, url, label}

function findVariantTab(variantId) {
    return openVariantTabs.find(function(tab) { return tab.variantId === variantId; });
}

function openVariantDetailsTab(variantId, url) {
    if (!isHorizontalMode()) {
        return false;
    }
    if (!findVariantTab(variantId)) {
        openVariantTabs.push({variantId: variantId, url: url, label: null});
    }
    showVariantDetailsTab(variantId);
    return true;
}

function loadVariantDetailsPane(tab) {
    const pane = $("<div>", {"class": "variant-details-pane"})
        .html('<div class="editor-loading"><i class="fa fa-spinner"></i> Loading variant details...</div>');
    $("#variant-details-container").empty().append(pane);
    pane.load(tab.url);
}

function showVariantDetailsTab(variantId) {
    const tab = findVariantTab(variantId);
    if (!tab) {
        return;
    }
    if (activeBottomPaneTab !== BOTTOM_PANE_VARIANT) {
        paneBeforeVariantTab = activeBottomPaneTab;
    }
    // The pane is also gone if the whole panel was reloaded under us (dual screen close)
    if (activeVariantTab !== variantId || !$(".variant-details-pane", "#variant-details-container").length) {
        activeVariantTab = variantId;
        loadVariantDetailsPane(tab);
    }
    renderVariantDetailsTabs();
    showBottomPaneTab(BOTTOM_PANE_VARIANT);
    openNodeEditorDrawer();
}

function closeVariantDetailsTab(variantId) {
    openVariantTabs = openVariantTabs.filter(function(tab) { return tab.variantId !== variantId; });
    if (activeVariantTab !== variantId) {
        renderVariantDetailsTabs();
        return;
    }

    activeVariantTab = null;
    $("#variant-details-container").empty();
    const nextTab = openVariantTabs[openVariantTabs.length - 1];
    if (nextTab) {
        showVariantDetailsTab(nextTab.variantId);
    } else {
        renderVariantDetailsTabs();
        showBottomPaneTab(paneBeforeVariantTab);
    }
}

/* The tabs belong to the grid they were opened from - anything that replaces it (selecting a node,
   a save, or toolbar content taking over the editor) takes its variant tabs with it */
function closeAllVariantDetailsTabs() {
    openVariantTabs = [];
    activeVariantTab = null;
    $("#variant-details-container").empty();
    renderVariantDetailsTabs();
    if (activeBottomPaneTab === BOTTOM_PANE_VARIANT) {
        showBottomPaneTab(paneBeforeVariantTab);
    }
}

// Called by the variant details page as it loads - until then the tab shows a placeholder
function setVariantDetailsTabLabel(variantId, label) {
    const tab = findVariantTab(variantId);
    if (tab) {
        tab.label = label;
        renderVariantDetailsTabs();
    }
}

function renderVariantDetailsTabs() {
    const docked = getEditorHome() === EDITOR_HOME_TAB;
    $("#bottom-pane-variant-tabs, #drawer-variant-tabs").empty();
    const strip = docked ? $("#bottom-pane-variant-tabs") : $("#drawer-variant-tabs");
    if (!strip.length) {
        return;
    }
    openVariantTabs.forEach(function(variantTab) {
        const tab = $("<a>", {
            "class": (docked ? "bottom-pane-tab" : "drawer-tab") + " variant-tab",
            "pane": BOTTOM_PANE_VARIANT,
            "variant_id": variantTab.variantId,
            href: "javascript:void(0)",
            text: variantTab.label || "Variant",
        }).appendTo(strip).click(function() {
            showVariantDetailsTab(variantTab.variantId);
        });
        $("<i>", {"class": "fa-solid fa-xmark close-variant-tab", title: "Close"}).appendTo(tab).click(function(event) {
            event.stopPropagation();
            closeVariantDetailsTab(variantTab.variantId);
        });
    });
    updateMirroredTabHighlight();
}

// Point the node editor (and so the data container, see switch_node_data) at one of its own tabs
function setNodeEditorTab(editorTab) {
    activeNodeEditorTab = editorTab;
    // ui-tabs-nav is added by the widget - the editor is only tabbed once it's finished loading
    const tabs = $("#node-editor-tabs:has(> ul.ui-tabs-nav)");
    if (tabs.length && tabs.tabs("option", "active") !== editorTab) {
        tabs.tabs("option", "active", editorTab);
    }
}

function showNodeEditorTab(editorTab) {
    setNodeEditorTab(editorTab);
    showBottomPaneTab(BOTTOM_PANE_EDITOR);
}

/* Which of the panes sharing the bottom panel is on screen. Docked, the panel shows exactly one of
   editor / variant details / grid; undocked the grid keeps the whole panel and the drawer body holds
   the editor or a variant details pane. */
function updatePaneVisibility() {
    const docked = bottomPaneTabsEnabled();
    const showVariant = activeBottomPaneTab === BOTTOM_PANE_VARIANT;
    const showEditor = docked ? activeBottomPaneTab === BOTTOM_PANE_EDITOR : !showVariant;

    $("#node-editor-container").toggle(showEditor);
    $("#variant-details-container").toggle(showVariant);
    $("#node-grid-container").toggle(!docked || activeBottomPaneTab === BOTTOM_PANE_GRID);
}

function gridPaneShown(userRequestedGrid) {
    if (deferredGridLoad && (userRequestedGrid || !deferredGridLoadNeedsUserRequest)) {
        const loadFunc = deferredGridLoad;
        deferredGridLoad = null;
        loadFunc();
    }
    resizeGrid();  // a table sized while hidden measures zero width
}

/* userRequestedGrid is the user asking for the grid itself (clicking the Grid tab) rather than the
   pane coming up for some other reason - that's the choice the big-grid placeholder asks for */
function showBottomPaneTab(name, userRequestedGrid) {
    if (!isHorizontalMode()) {
        return;
    }
    activeBottomPaneTab = name;
    updatePaneVisibility();
    updateMirroredTabHighlight();
    if ($("#node-grid-container").is(":visible")) {
        gridPaneShown(userRequestedGrid);
    }
}

function openNodeEditorDrawer() {
    if (getEditorHome() === EDITOR_HOME_DRAWER) {
        const drawer = $("#node-editor-drawer");
        drawer.css("top", $("#analysis-toolbar").outerHeight() + "px");  // start below the toolbar
        drawer.show();
    }
}

function closeNodeEditorDrawer() {
    $("#node-editor-drawer").hide();
}

/* Moves the editor into its current home and shows/hides the tab strip and drawer to match.
   Runs on every editor+grid load, so it also re-homes a container replaced by a panel reload. */
function applyEditorHome() {
    if (!isHorizontalMode()) {
        return;
    }
    const editorContainer = $("#node-editor-container");
    const variantContainer = $("#variant-details-container");
    const drawerBody = $("#node-editor-drawer-body");
    const tabs = $("#bottom-pane-tabs");

    if (getEditorHome() === EDITOR_HOME_DRAWER) {
        tabs.hide();
        drawerBody.children().not(editorContainer).not(variantContainer).remove();
        editorContainer.appendTo(drawerBody);
        variantContainer.appendTo(drawerBody);
    } else {
        closeNodeEditorDrawer();
        editorContainer.prependTo($("#grid-and-editor-container"));
        variantContainer.insertAfter(editorContainer);
        tabs.show();
    }
    updatePaneVisibility();
    // The editor's tabs and any open variant details mirror into whichever bar its new home has
    syncNodeEditorTabs();
    renderVariantDetailsTabs();
    if ($("#node-grid-container").is(":visible")) {
        gridPaneShown();
    }
}

function setEditorHome(home) {
    localStorage.setItem(EDITOR_HOME_STORAGE_KEY, home);
    applyEditorHome();
    if (home === EDITOR_HOME_DRAWER && $(".node-editor, #node-editor-wrapper", "#node-editor-container").length) {
        openNodeEditorDrawer();  // keep whatever the user was looking at on screen
    }
}

function setupNodeEditorDrawer() {
    const drawer = $("#node-editor-drawer");
    if (!drawer.length) {
        return;
    }

    $("#dock-editor-as-tab-button", drawer).click(function() {
        setEditorHome(EDITOR_HOME_TAB);
        showBottomPaneTab(BOTTOM_PANE_EDITOR);
    });
    $("#close-editor-drawer-button", drawer).click(closeNodeEditorDrawer);

    $(document).on("keydown", function(event) {
        const ESCAPE_KEY = 27;
        if (event.keyCode === ESCAPE_KEY && drawer.is(":visible")) {
            closeNodeEditorDrawer();
        }
    });

    // Drag the left edge to resize
    $(".drawer-resize-handle", drawer).on("mousedown", function(event) {
        event.preventDefault();
        const startX = event.pageX;
        const startWidth = drawer.outerWidth();

        const onMove = function(moveEvent) {
            const width = Math.max(startWidth + startX - moveEvent.pageX, DRAWER_MIN_WIDTH);
            drawer.css("width", width + "px");
        };
        const onUp = function() {
            $(document).off("mousemove", onMove).off("mouseup", onUp);
        };
        $(document).on("mousemove", onMove).on("mouseup", onUp);
    });
}

function setupBottomPaneTabs() {
    const tabs = $("#bottom-pane-tabs");
    if (!tabs.length) {
        return;
    }
    $(".bottom-pane-tab", tabs).click(function() {
        const pane = $(this).attr("pane");
        if (pane === BOTTOM_PANE_EDITOR) {
            showNodeEditorTab(NODE_EDITOR_TAB_EDITOR);
        } else {
            setNodeEditorTab(NODE_EDITOR_TAB_EDITOR);  // clicking Grid asks for the variant grid
            showBottomPaneTab(pane, pane === BOTTOM_PANE_GRID);
        }
    });
    $("#undock-editor-button", tabs).click(function() {
        setEditorHome(EDITOR_HOME_DRAWER);
        openNodeEditorDrawer();
    });
}

function viewTags() {
    unselectActive();
    loadGridAndEditorForNode(ANALYSIS_TAGS_NODE_ID);
}

function replaceEditorWindow(url) {
    const nodeEditorContainer = $("#node-editor-container");
    nodeEditorContainer.empty();
    $("#error-container").empty();
    const nodeDataContainer = $("#node-data-container");
    nodeDataContainer.empty();
    nodeDataContainer.removeAttr("node_url");
    nodeDataContainer.removeAttr("node_id");
    unselectActive();
    if (url) {
        // Toolbar content (analysis settings, input samples) is editor-only - put it where it can be seen
        revealNodeEditor();
        nodeEditorContainer.load(url);
    } else {
        closeNodeEditorDrawer();
    }
}

/* Content that only has an editor half (analysis settings, input samples, variant details) - bring
   the editor up in whichever home it's in */
function revealNodeEditor() {
    closeAllVariantDetailsTabs();
    syncNodeEditorTabs();  // what's loading has no tabs of its own
    showBottomPaneTab(BOTTOM_PANE_EDITOR);
    openNodeEditorDrawer();
}

function analysisSettings() {
    replaceEditorWindow(Urls.analysis_settings(ANALYSIS_ID));
}

function inputSamples() {
    replaceEditorWindow(Urls.analysis_input_samples(ANALYSIS_ID));
}


function layoutAnalysisPanels(showAnalysisVariables, initialAnalysisPanelFraction, nodeDataArray, nodeConnections, readOnly) {
    if (showAnalysisVariables) {
        $("#analysis-variables").show();
        Split(['#analysis-variables', '#analysis-and-toolbar-container'], {
            sizes: [15, 85],
            direction: 'vertical',
        });
    }

    const onDragEnd = function() {
        resizeGrid();
        if (!readOnly) {
            resizePanel();   // save panel widths
        }
    };
    const initialGridAndEditorFraction = 1.0 - initialAnalysisPanelFraction;
    const splitParams = {
        sizes: [initialAnalysisPanelFraction * 100, initialGridAndEditorFraction * 100],
        expandToMin: true,
        onDragEnd: onDragEnd,
    };
    if (isHorizontalMode()) {
        // Canvas on top, editor/grid along the bottom - both full window width
        splitParams.direction = 'vertical';
        // The canvas is full width now, so the counts legend goes in the toolbar beside the genome build
        $("#node-count-legend").insertAfter($("#analysis-genome-build"));
    }
    Split(['#left-panel', '#right-panel'], splitParams);

    setupNodeEditorDrawer();

    // Make clicking the background send a "click event" but not if you click on the .window
    // This is so we can make the editor switch out and then add nodes etc....
    $('#analysis-container').click(function() {

    });
    $('#analysis-container .window').click(function(event) {
        event.stopPropagation();
    });

    // Setup nodes
    addNodesToDOM('#analysis', nodeDataArray, readOnly);

    window.variantgridPipeline.init(readOnly);

    attatchAnalysisNodeConnections(nodeConnections, readOnly);

    messagePoller.update_loop();
}


/* Error handling */
function createJSEvent(details, severity, suppressErrors) {
    let data = 'app_name=analysis&event_name=javascript_error';
    data += '&severity=' + severity + '&details=' + encodeURIComponent(details);
    $.ajax({
        type: "POST",
        data: data,
        url: Urls.create_event(),
        suppressErrors: suppressErrors,
    });
}


function createJSErrorEvent(details, suppressErrors) {
    createJSEvent(details, 'E', suppressErrors);
}


function ajaxError(event, jqxhr, settings, thrownError) {
    //console.log("ajaxError:");
    //console.log(jqxhr);
    //console.log(settings);
    // #1123 - Chrome rendering error can be safely ignored
    if (event === 'ResizeObserver loop limit exceeded') {
        return;
    }

    if(settings.suppressErrors) {
        return;
    }

    const rj = jqxhr.responseJSON;
    if (rj) {
        if (rj.non_fatal) {
            console.log("Ignoring error");

            // Try and fix state up a bit
            if (rj.deleted_nodes) {
                deleteNodesFromDOM(rj.deleted_nodes, []);
            }

            return;
        }
    }

    // Ignore errors from navigating away from page
    // See https://stackoverflow.com/questions/9229005/how-to-handle-jquery-ajax-post-error-when-navigating-away-from-a-page
    if (jqxhr.readyState >= 4) {
        let message = "Error making request to server.";
        if (jqxhr.readyState == 4) {
            message += " status: " + jqxhr.status + " statusText: " + jqxhr.statusText;
        } else {
            message += "\n Unknown error, readyState: " + jqxhr.readyState;
        }
        createJSErrorEvent(message, true); // suppressErrors=True so we don't keep re-triggering errors
        showReloadPageErrorDialog($("#error-dialog"), message, true);
    }
}

function setupErrorHandlers() {
    $(document).ajaxError(ajaxError);

    window.onerror = function(msg, url, line) {
        // #1123 - Chrome rendering error can be safely ignored 
        if (msg === 'ResizeObserver loop limit exceeded') {
            return;
        }

        let details = "Error: " + msg + "\nurl: " + url + "\nline: " + line;
        details += "\nBrowser: " + navigator.userAgent;
        createJSErrorEvent(details);

        let userMessage = "<p>There was an problem running scripts on the page. Sorry about that!";
        userMessage += "<p>The error has been reported. In the mean time you can try reloading the page and avoiding whatever caused the error. ";
        userMessage += "<p>You could also try pressing <B>SHIFT</B> then click the <b>RELOAD</b> button (circular arrow) in your browswer to update your cache.";
        userMessage += "<p>Error was: <pre>" + details + "</pre>";
        showReloadPageErrorDialog($("#error-dialog"), userMessage, true);

        const suppressErrorAlert = false;
        return suppressErrorAlert;
    };
}

function setupNodeTypeSelect() {
    // Icon and source/filter colour come from NODE_TYPES - see node_types.get_node_display_data_by_class_name()
    function renderNodeTypeItem(className, label) {
        // Class name on the row picks up the node's accent colour - see analysis_nodes.css
        const wrapper = $("<div>", {"class": "node-type-item " + className});
        const nodeType = NODE_TYPES[className];
        if (nodeType) {
            wrapper.attr("node_classification", nodeType.classification);
        }
        renderNodeIcon(nodeType && nodeType.icon).appendTo(wrapper);
        $("<span>", {text: label}).appendTo(wrapper);
        return wrapper;
    }

    $.widget( "custom.iconselectmenu", $.ui.selectmenu, {
      _renderItem: function( ul, item ) {
        const li = $( "<li>" );
        if ( item.disabled ) {
          li.addClass( "ui-state-disabled" );
        }
        return li.append(renderNodeTypeItem(item.element.attr("value"), item.label)).appendTo( ul );
      },

      _renderButtonItem: function( item ) {
        return renderNodeTypeItem(item.element.attr("value"), item.label)
            .addClass("ui-selectmenu-text");
      }
    });

    $("#id_node_types").iconselectmenu({width: false})
        .iconselectmenu( "menuWidget" )
        .addClass( "node-type-menu" );

}

function addVariantTag(variantId, nodeId, tagId, successFunc) {
    setVariantTag(variantId, nodeId, tagId, successFunc, 'add');
}

function removeVariantTag(variantId, tagId, successFunc) {
    setVariantTag(variantId, null, tagId, successFunc, 'del');
}

function setNumVariantTags(pulseLabel) {
    pulseLabel = getValue(pulseLabel, true);
    const numberOfTags = $("#number-of-tags");
    const numTags = Object.keys(variantTags).length;
    let label = "";
    if (numTags) {
        label = numTags;
        if (pulseLabel) {
            numberOfTags.addClass("red-pulse");
            setTimeout(function() {
                numberOfTags.removeClass("red-pulse");
            }, 2000);
        }
    }
    numberOfTags.text(label);
}

function setVariantTag(variantId, nodeId, tagId, successFunc, op) {
    let data = 'variant_id=' + variantId;
    data += '&tag_id=' + tagId;
    data += '&op=' + op;
    data += '&analysis_id=' + ANALYSIS_ID;
    if (nodeId) {
        data += '&node_id=' + nodeId;
    }

    const success = function () {
        const aWin = getAnalysisWindow();
        const tagList = aWin.variantTags[variantId] || [];
        if (op == 'add') {
            tagList.push(tagId);
        } else if (op == 'del') {
            if (tagList) {
                removeItemFromArray(tagId, tagList);
            }
        }
        if (Object.keys(tagList).length > 0) {
            aWin.variantTags[variantId] = tagList;
        } else {
            delete aWin.variantTags[variantId];
        }
        setNumVariantTags();
        checkAndMarkDirtyNodes(aWin);

        if (successFunc) {
            successFunc();
        }
    };

    $.ajax({
        type: "POST",
        data: data,
        url: Urls.set_variant_tag('A'),
        success : success,
    });
}

function analysisVariable(nodeId, field, op, successCallback) {
    const url = Urls.analysis_template_variable(ANALYSIS_ID, nodeId);
    const data = 'field=' + field + '&op=' + op;
    $.ajax({
        type: "POST",
        data: data,
        url: url,
        success: successCallback,
    });
}

function getNodeFieldWrapper(nodeId, field) {
    const nodeEditor = $("#node-editor-wrapper[node_id=" + nodeId + "]");
    return $(".analysis-variable-node-field-wrapper[field=" + field + "]", nodeEditor);
}

function lockNodeField(nodeId, field, lock) {
    const NODE_FIELD_LOCKED_CSS_CLASS = "node-field-locked";
    const nodeFieldWrapper = getNodeFieldWrapper(nodeId, field);
    const addVariableButton = $(".add-analysis-variable-button", nodeFieldWrapper);

    if (lock) {
        nodeFieldWrapper.addClass(NODE_FIELD_LOCKED_CSS_CLASS);
        addVariableButton.hide();
    } else {
        nodeFieldWrapper.removeClass(NODE_FIELD_LOCKED_CSS_CLASS);
        addVariableButton.show();
    }
    $("#id_" + field, nodeFieldWrapper).attr("disabled", lock);
}

function addAnalysisVariableButton(nodeId, field, readOnly) {
    const ANALYSIS_VARIABLES_HELP_DIV = $("#analysis-variables-help");
    const NODE_CONTAINER_CLASS = "analysis-variable-node";
    const ANALYSIS_VARIABLE_CLASS = "analysis-variable";

    const nodeFieldSet = analysisNodeVariables[nodeId] || new Set();
    analysisNodeVariables[nodeId] = nodeFieldSet;
    nodeFieldSet.add(field);

    $("#av-example-add-button", ANALYSIS_VARIABLES_HELP_DIV).button({icon: "ui-icon-arrowthick-1-n"});
    $("#av-example-button", ANALYSIS_VARIABLES_HELP_DIV).button();

    const node = getNode(nodeId);
    const avContainer = $("#analysis-variables");
    let existingNodeContainer = $("." + NODE_CONTAINER_CLASS + "[node-id=" + nodeId + "]", avContainer);
    if (existingNodeContainer.length === 0) {
        existingNodeContainer = $("<div />").addClass("left").addClass(NODE_CONTAINER_CLASS).attr("node-id", nodeId);
        $(".node-overlay", node).addClass("variable-node");
        let nodeName = $(".node-name", node).text();
        if (!nodeName) {
            nodeName = node.attr("node_class") + ": " + nodeId;
        }
        const nodeTitle = $("<div />").append($("<b/>").html(nodeName));
        existingNodeContainer.append(nodeTitle);
        $("#insert-analysis-variables-before-here", avContainer).before(existingNodeContainer);
    }

    let fieldAnalysisVariable = $("." + ANALYSIS_VARIABLE_CLASS + "[field=" + field + "]", existingNodeContainer);
    if (fieldAnalysisVariable.length === 0) {
        const attrDict = {"field": field};
        fieldAnalysisVariable = $("<button />").addClass("btn btn-primary-outline").addClass(ANALYSIS_VARIABLE_CLASS).attr(attrDict).html(field);
        if (!readOnly) {
            fieldAnalysisVariable.click(function () {
                function removeVariableButton() {
                    nodeFieldSet.delete(field);
                    // TODO: Also node
                    fieldAnalysisVariable.remove();
                    // If empty, remove node container
                    if ($("." + ANALYSIS_VARIABLE_CLASS, existingNodeContainer).length === 0) {
                        existingNodeContainer.remove();
                        $(".node-overlay", node).removeClass("variable-node");
                    }
                    lockNodeField(nodeId, field, false);
                    checkTemplateSave();
                }

                analysisVariable(nodeId, field, 'del', removeVariableButton);
            });
        }
        existingNodeContainer.append(fieldAnalysisVariable);
    }
    checkTemplateSave();
}

function checkTemplateSave() {
    const analysisVariables = $(".analysis-variable-node", "#analysis-variables");
    const templateSave = $("button#analysis-template-save-version", "#analysis-template-version");
    if (analysisVariables.length) {
        templateSave.prop("disabled", false);
        templateSave.prop("title", "");
    } else {
        templateSave.prop("disabled", true);
        templateSave.prop("title", "Cannot save - no Analysis variables!");
    }
}

function setupAnalysisTemplateTopBar(analysisTemplateId) {
    const templateInfo = $("#analysis-template-info");
    const atVersion = $("#analysis-template-version");
    $("button#analysis-template-save-version", atVersion).button().click(function() {
        const analysisNameTemplate = $("#id_analysis_name_template").val();

        atVersion.hide();
        const savingMessage = $("<div />").html("saving...");
        templateInfo.append(savingMessage);

        $.ajax({
            type: "POST",
            data: {analysis_name_template: analysisNameTemplate},
            url: Urls.analysis_template_save(analysisTemplateId),
            success: function(data) {
                savingMessage.remove();
                let message = "Save completed";
                let severity = "info";
                let messageTime = 1000;

                if (data.version) {
                    const savedDate = data.created ? new Date(data.created).toLocaleString() : "";
                    const versionText = "v." + data.version + (savedDate ? " saved " + savedDate : "");
                    $("#latest-template-version").html(versionText);
                }
                if (data.error) {
                    message = data.error;
                    severity = "error";
                    messageTime = 5000;
                }
                const sm = $("<li />", {class: severity}).html(message);
                const saveMessage = $("<ul />", {class: 'messages'}).append(sm);
                templateInfo.append(saveMessage);
                saveMessage.fadeOut(messageTime, function () {
                    atVersion.fadeIn();
                });
            },
        });
    });
    checkTemplateSave();
}

function addInitialAnalysisVariables(analysisVariablesArray, readOnly) {
    for (let i=0 ; i<analysisVariablesArray.length; ++i) {
        const av = analysisVariablesArray[i];
        addAnalysisVariableButton(av[0], av[1], readOnly);
    }
}


function _getAnalysisWindow() {
   if (opener && opener.secondWindow == window) {
       return opener;
   }
   return window;
}

function loadGridAndEditorForNode(nodeId, extra_filters, fromSelectNode) {
    /* fromSelectNode - means we came from clicking on analysis (ie not reload etc)
        so won't reload if already selected */

    const gridAndEditorContainer = $("#grid-and-editor-container");
    const dataContainer = $("#node-data-container", gridAndEditorContainer);
    if (nodeId) {
        let load_node_url = Urls.node_load(ANALYSIS_ID, nodeId);
        if (extra_filters) {
            load_node_url += "?extra_filters=" + extra_filters;
        }

        if (fromSelectNode) {
            if (dataContainer.attr("node_url") === load_node_url) {
                return;
            }
        }

        dataContainer.attr("node_url", load_node_url);
        removeGridLoadingOverlay();  // tear down any in-progress grid overlay from the previous node
        registerDeferredGridLoad(null);  // the pending load belonged to the node we're leaving
        $("#node-editor-container").empty();
        closeAllVariantDetailsTabs();  // they belonged to the grid we're replacing
        showNodeEditorTab(NODE_EDITOR_TAB_EDITOR);  // a node comes up on its editor, so its grid stays deferred
        openNodeEditorDrawer();  // selecting a node is what opens the drawer
        showLoadingOverlay();
        dataContainer.load(load_node_url, function() {
            $(this).attr('node_id', nodeId);
        });
    } else {
        removeGridLoadingOverlay();
        closeNodeEditorDrawer();
        closeAllVariantDetailsTabs();
        $("#node-editor-container").html("Please select a node");
        dataContainer.empty();
    }
}

function layout_analysis_editor_and_grid() {
    //console.log("editor_and_grid layout panels...");

    if (typeof(resizePanel) == 'undefined') { // Defined in analysis, missing in stand alone
        //console.log("resizePanel = null");
        resizePanel = null;
    }

    setupBottomPaneTabs();
    applyEditorHome();
}


function cancelNodeLoad(analysisId, nodeId) {
	$.ajax({
	    type: "POST",
	    url: Urls.node_cancel_load(analysisId, nodeId),
	    success: function(data) {
	    	console.log("Stopped loading");
		}
	});
}

function getLoadedNodeId() {
    const loadedNodeId = $('#node-data-container').attr('node_id');
    return parseInt(loadedNodeId);
}


function showLoadingOverlay() {
    const oc = $("#overlay-container");
    if (!oc.is(":visible")) {
        // Move to right-panel (with top z-order), then things can load underneath.
        oc.appendTo("#right-panel");
        oc.show();

        $("#loading-message").remove();

        const canvasAttributes = {class: 'node-load-animation'};
        const container = $("#animation-container");
        const canvas = $("<canvas />", canvasAttributes);
        canvas.attr({ width: 50, height: 180 });
        canvas.css('opacity', 0.35);
        canvas.DoubleHelix({fps: 20, spinSpeed: 4});
        container.append(canvas);
    }
}

// Hide overlay by fading out and removing canvas (don't want to animate it anymore)
function hideLoadingOverlay() {
    $("#overlay-container").fadeOut(function() {
        $('canvas.node-load-animation').each(function() {
             this.active = false;
        });
        $("#animation-container", this).empty();
        window.showingLoadOverlay = false;
    });
}


// Loading overlay scoped to a single container (e.g. the grid section) so the rest of the node
// (the editor) stays visible and usable. Builds its own double-helix inside the container rather
// than reparenting the shared #overlay-container - which flashed over the whole panel before
// snapping into place.
function showGridLoadingOverlay(container) {
    if ($("#grid-loading-overlay").length) {
        return;
    }
    const $container = $(container);
    // Guarantee the container is the positioning context so the absolute overlay fills only it
    // (not a wider ancestor like #right-panel), and that it has room for the helix even when the
    // grid is still empty mid-load. Don't rely on the CSS rule alone - the panel HTML can be
    // served from browser cache without it.
    if ($container.css("position") === "static") {
        $container.css("position", "relative");
    }
    if ($container.height() < 200) {
        $container.css("min-height", "200px");
    }
    const overlay = $("<div>", {id: "grid-loading-overlay"});
    // Pin the overlay to the visible grid area (panel height minus the editor above it). A fixed
    // height keeps the centred helix from jumping down as rows load and the container grows.
    const available = gridPaneAvailableHeight();
    if (available > 200) {
        overlay.css("height", available + "px");
    }
    $container.append(overlay);

    // Show one of the user's chosen full-size DNA "reading" effects (dark-mode glyphs on a clear
    // background, over the overlay's white), picked at random. ANALYSIS_LOADING_ANIMATIONS holds the
    // user's menu ids; map them to VGLoaders ids here. Fall back to the double-helix if the loaders
    // aren't available.
    const LOADER_IDS = {flowcell: "flowcell", ripple: "base-fx-ripple", pileup: "pileup", matrix: "base-fx-matrix"};
    let menu = [];
    if (typeof VGLoaders !== "undefined") {
        const prefs = (typeof ANALYSIS_LOADING_ANIMATIONS !== "undefined" && ANALYSIS_LOADING_ANIMATIONS) || [];
        menu = prefs.filter(a => LOADER_IDS[a]);
        if (!menu.length) {
            menu = ["flowcell", "ripple", "pileup"];  // safety default if prefs are empty/unknown
        }
    }

    if (menu.length) {
        const choice = menu[Math.floor(Math.random() * menu.length)];
        const stage = $("<div>", {class: "grid-loading-stage"});
        stage.css({position: "absolute", top: 0, left: 0, right: 0, bottom: 0});
        overlay.append(stage);
        overlay.data("stopLoader", VGLoaders.start(LOADER_IDS[choice], stage[0], {theme: "dark", clearBackground: true}));
    } else {
        // Fallback: classic double-helix.
        const canvas = $("<canvas />", {class: 'node-load-animation'});
        canvas.attr({width: 50, height: 180});
        canvas.css('opacity', 0.35);
        canvas.DoubleHelix({fps: 20, spinSpeed: 4});
        overlay.append(canvas);
    }
}

function hideGridLoadingOverlay() {
    // Fade out (like the original hideLoadingOverlay) so the loaded grid eases in rather than jumping.
    const overlay = $("#grid-loading-overlay");
    const stopLoader = overlay.data("stopLoader");
    overlay.fadeOut(function() {
        if (typeof stopLoader === "function") {
            stopLoader();  // stop the ripple/flowcell animation
        }
        $(this).find('canvas.node-load-animation').each(function() {
            this.active = false;  // stop the helix animation before removing (fallback)
        });
        $(this).remove();
    });
}

// Synchronous teardown (no fade) for when we abandon a load - e.g. switching to another node
// before the grid finished. Stops the animation loop so it doesn't keep running detached.
function removeGridLoadingOverlay() {
    const overlay = $("#grid-loading-overlay");
    if (!overlay.length) {
        return;
    }
    const stopLoader = overlay.data("stopLoader");
    if (typeof stopLoader === "function") {
        stopLoader();  // cancel the ripple/flowcell requestAnimationFrame
    }
    overlay.find('canvas.node-load-animation').each(function() {
        this.active = false;  // stop the helix animation (fallback)
    });
    overlay.remove();
}


function finishedLoadingEditor(node_id, version_id) {
    const everythingLoaded = function () {
        hideLoadingOverlay();
    };

    const unique_code = node_id + "_" + version_id; // make sure only attach editor to grid that requested
    registerComponent(unique_code, EDITOR, everythingLoaded);
}


function setVisibleSliderValue(inputSelector, sliderSelector, value) {
    const container = sliderSelector.parents(".slider-container");
    const sliderValue = $(".slider-value", container);
    let decimalPlaces = inputSelector.attr("decimal_places");
    if (typeof decimalPlaces == 'undefined') {
        decimalPlaces = 2;
    }
    const floatVal = parseFloat(value);
    if (!isNaN(floatVal)) {
        value = floatVal.toFixed(decimalPlaces);
    }
    sliderValue.html(value);
}

function setupSlider(inputSelector, sliderSelector) {
    // Returns sliderValue
    const sliderMinVal = Number(inputSelector.attr("min"));
    const sliderMaxVal = Number(inputSelector.attr("max"));
    const sliderVal = inputSelector.val();
    const container = sliderSelector.parents(".slider-container");

    setVisibleSliderValue(inputSelector, sliderSelector, sliderVal); // set initial

    sliderSelector.slider({
        min: sliderMinVal,
        max: sliderMaxVal,
        step: Number(inputSelector.attr("step") || 1),
        value: sliderVal,
        change: function (event, ui) {
            inputSelector.val(ui.value);
        },
        slide: function (event, ui) {
            inputSelector.val(ui.value);
            setVisibleSliderValue(inputSelector, sliderSelector, ui.value);
        },
    });

    $(".min-value", container).html(sliderMinVal);
    $(".max-value", container).html(sliderMaxVal);
}
