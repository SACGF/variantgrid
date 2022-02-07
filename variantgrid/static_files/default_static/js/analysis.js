function startDualScreenMode() {
    const STAND_ALONE_URL = Urls.standalone_analysis_editor_and_grid(ANALYSIS_ID);
    if (!secondWindow) {
        saveSettingsOnResize = false;
        secondWindow = window.open(STAND_ALONE_URL, '_blank', 'toolbar=0,location=0,menubar=0');
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
    let win = getGridAndEditorWindow();
    win.loadGridAndEditorForNode(nodeId, extra_filters, fromSelectNode);
}

function resizeGrid(pane, $Pane, paneState) {
    let grid = $('.ui-jqgrid-btable:visible');
    if(grid.length) {
        grid.each(function(index) {
            let gridId = $(this).attr('id');
            $('#' + gridId).setGridWidth(paneState.innerWidth - 2);
        });
    }
}

function savePanelWidthSettings() {
    if (saveSettingsOnResize) {
        let analysis_panel_fraction = $("#analysis-and-toolbar-container").outerWidth() / window.innerWidth;
        let data = 'analysis_panel_fraction=' + analysis_panel_fraction;
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

function viewTags() {
    unselectActive();
    loadGridAndEditorForNode(ANALYSIS_TAGS_NODE_ID);
}

function replaceEditorWindow(url) {
    let nodeEditorContainer = $("#node-editor-container");
    nodeEditorContainer.empty();
    $("#error-container").empty();
    let nodeDataContainer = $("#node-data-container");
    nodeDataContainer.empty();
    nodeDataContainer.removeAttr("node_url");
    nodeDataContainer.removeAttr("node_id");
    unselectActive();
    if (url) {
        nodeEditorContainer.load(url);
    }
}

function analysisSettings() {
    replaceEditorWindow(Urls.analysis_settings(ANALYSIS_ID));
}

function inputSamples() {
    replaceEditorWindow(Urls.analysis_input_samples(ANALYSIS_ID));
}


function layoutAnalysisPanels(showAnalysisVariables, initialGridAndEditorWidth, nodeDataArray, nodeConnections, readOnly) {
    if (showAnalysisVariables) {
        $("#analysis-variables").show();
    }

    const centerLayoutParams = {
        minWidth: 200,
    };
    if (!readOnly) {
        centerLayoutParams.onresize = resizePanel;  // save panel widths
    }

    $('div#content').layout();
    $('div#analysis-outer-container').layout({
        north: {
            size: "15%",
            spacing_closed: 0, // HIDE resizer & toggler when 'closed'
            slidable: false,
            initClosed: !showAnalysisVariables,
        },
        center: centerLayoutParams,
        east: { onresize: resizeGrid,
                triggerEventsOnLoad: true,
                size: initialGridAndEditorWidth,
                // Setting minSize in this pane cases 'InternalError: too much recursion' with sizeMidPanes
                // So we'll just reset min size upon loading each time
        }
    });
    $('div#analysis-and-toolbar-container').layout({
        north: {
            minSize: 32,
        },
    });

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

    let rj = jqxhr.responseJSON;
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

        let suppressErrorAlert = false;
        return suppressErrorAlert;
    };
}

function setupNodeTypeSelect() {
    $.widget( "custom.iconselectmenu", $.ui.selectmenu, {
      _renderItem: function( ul, item ) {
        let li = $( "<li>" ),
          wrapper = $( "<div>", { text: item.label } );

        if ( item.disabled ) {
          li.addClass( "ui-state-disabled" );
        }

        $( "<span>", {
          style: item.element.attr( "data-style" ),
          "class": "ui-icon " + item.element.attr( "value" )
        })
          .appendTo( wrapper );

        return li.append( wrapper ).appendTo( ul );
      }
    });

    $("#id_node_types").iconselectmenu()
        .iconselectmenu( "menuWidget" )
        .addClass( "ui-menu-icons customicons" );

}

function addVariantTag(variantId, nodeId, tagId, successFunc) {
    setVariantTag(variantId, nodeId, tagId, successFunc, 'add');
}

function removeVariantTag(variantId, tagId, successFunc) {
    setVariantTag(variantId, null, tagId, successFunc, 'del');
}

function setNumVariantTags(pulseLabel) {
    pulseLabel = getValue(pulseLabel, true);
    let numberOfTags = $("#number-of-tags");
    let numTags = Object.keys(variantTags).length;
    let label = "";
    if (numTags) {
        label = "x" + numTags;
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
        let tagList = aWin.variantTags[variantId] || [];
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
    let url = Urls.analysis_template_variable(nodeId);
    let data = 'field=' + field + '&op=' + op;
    $.ajax({
        type: "POST",
        data: data,
        url: url,
        success: successCallback,
    });
}

function getNodeFieldWrapper(nodeId, field) {
    let nodeEditor = $("#node-editor-wrapper[node_id=" + nodeId + "]");
    return $(".analysis-variable-node-field-wrapper[field=" + field + "]", nodeEditor);
}

function lockNodeField(nodeId, field, lock) {
    const NODE_FIELD_LOCKED_CSS_CLASS = "node-field-locked";
    let nodeFieldWrapper = getNodeFieldWrapper(nodeId, field);
    let addVariableButton = $(".add-analysis-variable-button", nodeFieldWrapper);

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

    let nodeFieldSet = analysisNodeVariables[nodeId] || new Set();
    analysisNodeVariables[nodeId] = nodeFieldSet;
    nodeFieldSet.add(field);

    $("#av-example-add-button", ANALYSIS_VARIABLES_HELP_DIV).button({icon: "ui-icon-arrowthick-1-n"});
    $("#av-example-button", ANALYSIS_VARIABLES_HELP_DIV).button();

    let node = getNode(nodeId);
    let avContainer = $("#analysis-variables");
    let existingNodeContainer = $("." + NODE_CONTAINER_CLASS + "[node-id=" + nodeId + "]", avContainer);
    if (existingNodeContainer.length === 0) {
        existingNodeContainer = $("<div />").addClass("left").addClass(NODE_CONTAINER_CLASS).attr("node-id", nodeId);
        $(".node-overlay", node).addClass("variable-node");
        let nodeName = $(".node-name", node).text();
        if (!nodeName) {
            nodeName = node.attr("node_class") + ": " + nodeId;
        }
        let nodeTitle = $("<div />").append($("<b/>").html(nodeName));
        existingNodeContainer.append(nodeTitle);
        $("#insert-analysis-variables-before-here", avContainer).before(existingNodeContainer)
    }

    let fieldAnalysisVariable = $("." + ANALYSIS_VARIABLE_CLASS + "[field=" + field + "]", existingNodeContainer);
    if (fieldAnalysisVariable.length === 0) {
        let attrDict = {"field": field};
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
    let analysisVariables = $(".analysis-variable-node", "#analysis-variables");
    let templateSave = $("button#analysis-template-save-version", "#analysis-template-version");
    if (analysisVariables.length) {
        templateSave.prop("disabled", false);
        templateSave.prop("title", "");
    } else {
        templateSave.prop("disabled", true);
        templateSave.prop("title", "Cannot save - no Analysis variables!");
    }
}

function setupAnalysisTemplateTopBar(analysisTemplateId) {
    let templateInfo = $("#analysis-template-info");
    let atVersion = $("#analysis-template-version");
    $("button#analysis-template-save-version", atVersion).button().click(function() {
        let analysisNameTemplate = $("#id_analysis_name_template").val();

        atVersion.hide();
        let savingMessage = $("<div />").html("saving...");
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
                    $("#latest-template-version").html(data.version);
                }
                if (data.error) {
                    message = data.error;
                    severity = "error";
                    messageTime = 5000;
                }
                let sm = $("<li />", {class: severity}).html(message)
                let saveMessage = $("<ul />", {class: 'messages'}).append(sm);
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
        let av = analysisVariablesArray[i];
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

    let gridAndEditorContainer = $("#grid-and-editor-container");
    let dataContainer = $("#node-data-container", gridAndEditorContainer);
    if (nodeId) {
        let load_node_url = Urls.node_load(nodeId);
        if (extra_filters) {
            load_node_url += "?extra_filters=" + extra_filters;
        }

        if (fromSelectNode) {
            if (dataContainer.attr("node_url") === load_node_url) {
                return;
            }
        }

        dataContainer.attr("node_url", load_node_url);
        $("#node-editor-container", gridAndEditorContainer).empty();
        showLoadingOverlay();
        dataContainer.load(load_node_url, function() {
            $(this).attr('node_id', nodeId);
        });
    } else {
        $("#node-editor-container", gridAndEditorContainer).html("Please select a node");
        dataContainer.empty();
    }
}

function layout_analysis_editor_and_grid() {
    //console.log("editor_and_grid layout panels...");

    if (typeof(resizePanel) == 'undefined') { // Defined in anaylsis, missing in stand alone
        //console.log("resizePanel = null");
        resizePanel = null;
    }
}

function getLoadedNodeId() {
    const loadedNodeId = $('#node-data-container').attr('node_id');
    return parseInt(loadedNodeId);
}


function showLoadingOverlay() {
    const oc = $("#overlay-container");
    if (!oc.is(":visible")) {
        // Move to right-panel (with top z-order), then things can load underneath.
        oc.show();
        oc.appendTo("#right-panel");
    
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


function finishedLoadingEditor(node_id, version_id) {
    const everythingLoaded = function () {
        hideLoadingOverlay();
    };

    const unique_code = node_id + "_" + version_id; // make sure only attach editor to grid that requested
    registerComponent(unique_code, EDITOR, everythingLoaded);
}


function setupSlider(inputSelector, sliderSelector, enableInput) {
    let sliderMinVal = Number(inputSelector.attr("min"));
    let sliderMaxVal = Number(inputSelector.attr("max"));
    let sliderVal = inputSelector.val();
    let container = sliderSelector.parents(".slider-container");
    let sliderValue = $(".slider-value", container);
    sliderValue.html(sliderVal); // set initial

    sliderSelector.slider({
        min: sliderMinVal,
        max: sliderMaxVal,
        step: Number(inputSelector.attr("step") || 1),
        value: sliderVal,
        change: function( event, ui ) {
            inputSelector.val(ui.value);
        },
        slide: function( event, ui ) {
            sliderValue.html(ui.value);
        },
    });

    $(".min-value", container).html(sliderMinVal);
    $(".max-value", container).html(sliderMaxVal);

    if (enableInput) {
        let row = enableInput.parents("tr");
        let filterRequired = $(".filter-required", row);
        let frInputs = $("input", filterRequired);

        enableInput.change(function () {
            if ($(this).is(":checked")) {
                sliderSelector.slider("enable");
                frInputs.prop("disabled", false);
            } else {
                sliderValue.html("");
                sliderSelector.slider({value: 0})
                sliderSelector.slider("disable");
                inputSelector.val("");
                frInputs.prop("disabled", true);
            }
        });
        if (sliderVal === '') {
            enableInput.prop("checked", false);
            frInputs.prop("disabled", true);
        } else {
            enableInput.prop("checked", true);
            frInputs.prop("disabled", false);
        }
    }

    if (sliderVal === '') {
        sliderSelector.slider("disable");
    } else {
        sliderSelector.slider("enable");
    }
}
