/* The node grid's <table>. One per node-version pane, so it's looked up inside the pane rather
   than by id alone - a stale pane may still be in the DOM while the next one loads. */
function getGrid(nodeId, unique_code) {
	return $("#grid-" + nodeId, "#" + unique_code);
}

// The DataTableDefinition driving a node's grid, by node id - set by setupNodeGrid
const nodeGridDefinitions = {};

function getNodeGridDefinition(nodeId) {
	return nodeGridDefinitions[nodeId];
}

function getNodeDataTable(nodeId, unique_code) {
	const table = getGrid(nodeId, unique_code);
	if (table.length && $.fn.DataTable.isDataTable(table)) {
		return table.DataTable();
	}
	return null;
}

function getAnalysisDownloadTracker() {
	return getAnalysisWindow().analysisDownloadTracker;
}

/* Everything a node export is identified by, plus the URL that launches it. Built off the grid's live
   ajax params so the export sees whatever the user has filtered the grid down to.

   Works before the rows have ever been fetched: the placeholder on a big node offers CSV/VCF without
   the user loading the grid, so the params come from the definition's postData when the table has
   yet to make a request. Paging and sorting are dropped - the export orders by genome position. */
function nodeGridExportInfo(analysisId, nodeId, unique_code, export_type, use_canonical_transcripts, caption) {
	const definition = getNodeGridDefinition(nodeId);
	const dataTable = getNodeDataTable(nodeId, unique_code);
	let gridParams = {};
	if (dataTable) {
		gridParams = $.extend({}, dataTable.ajax.params());
	} else if (definition) {
		gridParams = $.extend({}, definition.serverParams.postData);
	}
	delete gridParams['start'];
	delete gridParams['order[0][column]'];
	delete gridParams['order[0][dir]'];
	gridParams['rows'] = 0; // no pagination
	gridParams['length'] = 0;
	gridParams['export_type'] = export_type;
	if (use_canonical_transcripts) {
		gridParams['use_canonical_transcripts'] = true;
	}

	const gridCaption = (definition && definition.serverParams.gridName) || ("Node " + nodeId);
	return {
		nodeId: nodeId,
		nodeVersion: gridParams['version_id'],
		exportType: export_type,
		useCanonicalTranscripts: Boolean(use_canonical_transcripts),
		gridParams: gridParams,
		caption: caption || export_type.toUpperCase(),
		label: `${gridCaption} (${export_type.toUpperCase()})`,
		url: Urls.node_grid_export(analysisId) + "?" + EncodeQueryData(gridParams),
	};
}

function export_grid(analysisId, nodeId, unique_code, export_type, use_canonical_transcripts) {
	const info = nodeGridExportInfo(analysisId, nodeId, unique_code, export_type, use_canonical_transcripts);
	const tracker = getAnalysisDownloadTracker();
	if (tracker) {
		tracker.download(info);
	} else {
		window.location = info.url;
	}
}

/* Keep a download button in sync with its export - a spinner and percentage while it runs, then a
   download icon. Safe to call for a node whose export was never started. */
function registerNodeGridDownloadButton(element, analysisId, nodeId, unique_code, export_type,
                                        use_canonical_transcripts, caption) {
	const tracker = getAnalysisDownloadTracker();
	if (tracker) {
		// Resolve here - in dual screen mode the tracker lives in the other window's document
		tracker.registerButton($(element), function() {
			return nodeGridExportInfo(analysisId, nodeId, unique_code, export_type,
			                          use_canonical_transcripts, caption);
		});
	}
}

function load_variant_details(variant_id) {
    const aWin = getAnalysisWindow();
    let variant_details_url = Urls.view_variant(variant_id);
    if (aWin.ANALYSIS_SETTINGS) {
        if (aWin.ANALYSIS_SETTINGS.open_variant_details_in_new_window) {
            const newWin = aWin.open(variant_details_url, '_blank');
            newWin.focus();
            return;
        } else {
            const annotation_version_id = aWin.ANALYSIS_SETTINGS["annotation_version"];
            if (annotation_version_id) {
                variant_details_url = Urls.variant_details_annotation_version(variant_id, annotation_version_id);
            }
        }
    }
    // Horizontal mode gives each variant its own closable tab, leaving the node editor where it was
    if (typeof openVariantDetailsTab === "function" && openVariantDetailsTab(variant_id, variant_details_url)) {
        return;
    }
    const editorContainer = $("#node-editor-container");
    editorContainer.html('<div class="editor-loading"><i class="fa fa-spinner"></i> Loading variant details...</div>');
    editorContainer.load(variant_details_url);
}

function getAnalysisWindow() {
    let aWin = window;
    if (typeof _getAnalysisWindow !== "undefined") {
        try {
            // Not sure what's going on here (leaking between tabs?) but I very rarely get
            // Error: Uncaught SecurityError: Blocked a frame with origin "X" from accessing a cross-origin frame.
            aWin = _getAnalysisWindow();
        } catch(e) {
            if (typeof RAISED_GET_ANALYSIS_WINDOW_JS_ERROR == "undefined") {
                RAISED_GET_ANALYSIS_WINDOW_JS_ERROR = true;
                console.log(e);
                const exception_string = e.message + '\n' + e.stack;
                createJSEvent(exception_string, 'W', true); // log to server
            }
        }
    } 
    return aWin;
}


function get_igv_data() {
    const aWin = getAnalysisWindow();
    return aWin.ANALYSIS_SETTINGS["igv_data"];
}

function replaceFilePrefix(replaceDict, bamFiles) {
    let replacedBamFiles = [];
    if (replaceDict) {
        for (let i=0 ; i<bamFiles.length ; ++i) {
            let bamFile = bamFiles[i];
            if (bamFile) {
                for (const fromValue in replaceDict) {
                    const toValue = replaceDict[fromValue];
                    if (bamFile.startsWith(fromValue)) {
                        bamFile = bamFile.replace(fromValue, toValue);
                        break;
                    }
                }
                replacedBamFiles.push(bamFile);
            }
        }
    } else {
        replacedBamFiles = bamFiles;
    }
    return replacedBamFiles;
}


function create_igv_url(locus, inputBams) {
    const IGV_DATA = get_igv_data();
    let url = IGV_DATA['base_url'];
    const params = ["genome=" + IGV_DATA['genome']];
    if (locus) {
    	params.push("locus=" + locus);
    }
    let bamFiles = [];
    const manual_zygosity_cohort = IGV_DATA["manual_zygosity_cohort"];
    if (manual_zygosity_cohort && manual_zygosity_cohort.length) {
    	for (let i=0 ; i<manual_zygosity_cohort.length ; ++i) {
    		bamFiles.push(manual_zygosity_cohort[i]);
    	}
    } else {
        bamFiles = inputBams;
    }

    let op = 'goto';
    if (bamFiles.length > 0) {
        const replaceDict = IGV_DATA["replace_dict"];
        const replacedBamFiles = replaceFilePrefix(replaceDict, bamFiles);
        const joinedFiles = replacedBamFiles.join();
        if (joinedFiles) {
    		params.push("file=" + joinedFiles);
    		op = "load";
    	}
	}
	url += '/' + op + '?';
	url += params.join("&");
	return url;
}


seen_igv_error = false;

function open_igv_link(locus, getBamsFunc) {
    const url = create_igv_url(locus, getBamsFunc);

    $.ajax({
        url: url,
        error: function(jqXHR, textStatus, errorThrown) {
            if (!seen_igv_error) {
                console.log(jqXHR);
                console.log(textStatus);
                console.log(errorThrown);

                const IGV_DATA = get_igv_data();
                const base_url = IGV_DATA['base_url'];
                const igvIntegrationUrl = Urls.igv_integration();

                let message = "<p>Could not connect to IGV - is it running and accepting connections on " + base_url + "?";
                message += "<p>See also <a target='_blank' href='" + igvIntegrationUrl + "'>IGV Integration</a>";
                
                $("#error-dialog").html(message).dialog({
                    minWidth: 500,
                    buttons: [
                        {   text: "OK",
                            click: function() {
                                $(this).dialog("close");
                            },
                        },
                    ],                
                });
                seen_igv_error = true;
            }
        },
        suppressErrors: true,
    });
}


function noBamsHere() {
    return [];    
}


function getViewVariantUrl(variantLink) {
    const variantId = $(variantLink).parent(".variant_id-container").attr("variant_id");
    return Urls.view_variant(variantId);
}

function setFullscreenVariantLink() {
    const url = getViewVariantUrl(this);
    $(this).attr('href', url);
}

function restoreVariantLink() {
    const orig_href = $(this).attr('orig_href');
    $(this).attr('href', orig_href);
}


function createGridLink(title, url, contents, extraLinkClasses, extraIconClasses) {
    let linkCSS = ['grid-link'];
    if (extraLinkClasses) {
        linkCSS = linkCSS.concat(extraLinkClasses);
    }
    let iconCSS = ['grid-link-icon', 'user-tag-colored'];
    if (extraIconClasses) {
        iconCSS = iconCSS.concat(extraIconClasses);
    }

    const gridBox = "<div class='" + iconCSS.join(' ') + "'>" + contents + "</div>";
    const link = $("<a/>").attr({
        "class": linkCSS.join(' '),
        "title": title,
        "orig_href": url,
        "href": url
    });
    link.append(gridBox);
    return link.prop("outerHTML");
}

function createIgvUrl(locus, getBamsFuncString) {
    const aWin = getAnalysisWindow();
    if (aWin.ANALYSIS_SETTINGS && aWin.ANALYSIS_SETTINGS['show_igv_links']) {
        if (!getBamsFuncString) {
            getBamsFuncString = 'noBamsHere';
        }
        return 'javascript:open_igv_link("' + locus + '", ' + getBamsFuncString + '())';
    }
    return null;
}

function create_igv_link(locus, getBamsFuncString) {
    const igvUrl = createIgvUrl(locus, getBamsFuncString);
    if (igvUrl) {
        return createGridLink("Open " + locus + " in IGV", igvUrl, '', [], ['igv-link']);
    }
    return '';
}

function showGridCell(gridColumn) {
    // Every cell carries a dt-<column> class (@see RichColumn.css_classes)
    const selector = $("td.dt-" + gridColumn);
    if (selector.length) {
        selector[0].scrollIntoView();
    }
}

function inAnalysis() {
    // Variant selection and filter child links are wired up by the analysis page - other pages showing
    // the same grids (sample variants tab, variantopedia) render them read only
    return typeof ANALYSIS_ID !== "undefined";
}


function showTagAutocomplete(variantId) {
    const container = $("#tag-entry-container-" + variantId);
    const addTagButton = $(".show-tag-autocomplete", container.parent());
    const nodeId = container.parents("#node-data-container").attr("node_id");
    
    addTagButton.hide();
    container.load(Urls.tag_autocomplete_form(), function() {
        const tagSelect = $("select#id_tag", container);
        tagSelect.change(function() {
            const tag = $(this).val();
            if (tag) {
                const successFunc = function () {
                    const vtHtml = getVariantTagHtml(variantId, tag);
                    const newTag = $(vtHtml);
                    newTag.click(tagClickHandler);
                    container.parent().append(newTag);
                    container.empty();
                    addTagButton.show();
                };
                addVariantTag(variantId, nodeId, tag, successFunc);
            }
        });

        // Can't call open() until element is fully initialised
        $(document).on("dal-element-initialized", function (e) {
            if (e.detail.element == tagSelect[0]) {
                tagSelect.select2("open").trigger("focus");
            }
        });
    });

}


function getVariantTagHtml(variantId, tag, readOnly, tagLabel, extraClasses, title) {
    if (typeof(tagLabel) === 'undefined') {
        tagLabel = tag;
    }
    if (typeof(title) === 'undefined') {
        title = `Tagged as ${tag}`;
    }
    const outerClasses = ["grid-tag", "tagged-" + tag];
    if (!readOnly) {
        outerClasses.push("grid-tag-deletable");
    }
    if (extraClasses) {
        outerClasses.push(...extraClasses);
    }
    return `<span class='${outerClasses.join(' ')}' title='${title}' variant_id='${variantId}' tag_id='${tag}'><span class='user-tag-colored'>${tagLabel}</span></span>`;
}


// Tags with no entry in variantTagOrder sort as 0, ties broken alphabetically.
// Customised per-collection on the tag colors page - see issue #343
function sortVariantTags(aWin, tagList) {
    const tagOrder = aWin.variantTagOrder || {};
    return tagList.slice().sort(function(a, b) {
        const diff = (tagOrder[a] || 0) - (tagOrder[b] || 0);
        if (diff) {
            return diff;
        }
        return a.localeCompare(b);
    });
}


function tagClickHandler() {
    const gridTag = $(this);
    const innerSpan = $(".user-tag-colored", gridTag);

    function removeClickHandler() {
        const tagId = gridTag.attr('tag_id');
        const variantId = gridTag.attr('variant_id');
        const removeTagCallback = function () {
            gridTag.remove();
        };
        removeVariantTag(variantId, tagId, removeTagCallback);
    }
    deleteItemClickHandler(gridTag, innerSpan, removeClickHandler);
}


// Kicked off after every grid draw, alongside the page's own gridComplete
function gridCompleteExtra() {
    const aWin = getAnalysisWindow();
    if (!aWin.variantTagsReadOnly) {
        $(".grid-tag-deletable").click(tagClickHandler);
    }

    // We want to be able to right click to open full screen link new tab
    // but normal click does JS / load() call to open in the editor.
    const variantLink = $('a.variant-link');
    variantLink.on("contextmenu", setFullscreenVariantLink);
    variantLink.on("mouseup", setFullscreenVariantLink);
    variantLink.on('click', restoreVariantLink);
    variantLink.on('mouseout', restoreVariantLink);

}


/* The analysis node grid. Two endpoints: config_url answers the table definition (columns, widths,
   renderers - it varies by node version, so it caches with the node) and the handler serves the rows
   for the whole analysis. The node's own state - node_id, version_id, ccc_id, ccc_version_id,
   extra_filters, zygosity_samples_hash and a FilterNode's filters - comes back on the definition as
   'postData' and goes up as the ajax params of every row request.

   autoLoad false builds the table but holds the row query back until loadNodeGridData() asks for it,
   which is how a node over the auto-load row count waits behind its placeholder.

   Resolves once the table is built (or the definition reported node errors), so the caller can wire
   up its row interactions. */
function setupNodeGrid(config_url, handler_url, analysisId, nodeId, versionId, unique_code,
                       gridComplete, gridLoadError, on_error_function, autoLoad) {
    if (typeof autoLoad === "undefined") { autoLoad = true; }

    const definition = new DataTableDefinition({
        dom: getGrid(nodeId, unique_code),
        definitionUrl: config_url,
        url: handler_url,
        // The widths come from the definition and the CSS lays the table out table-layout: fixed, so
        // there is nothing to measure - and on a deferred grid the adjust draw is the fetch we're holding
        adjustColumns: false,
        // FilterNode rules are saved against the node, so they ride along with its other state
        data: function(data) {
            return $.extend({}, definition.serverParams.postData);
        },
        onDefinition: function(defn) {
            if (defn.errors) {
                on_error_function(defn.errors);
                return false;
            }
            return true;
        },
        onData: function(json) {
            if (json.non_fatal && json.deleted_nodes) {
                deleteNodesFromDOM(json.deleted_nodes, []);
                if (json.message) {
                    $("#node-editor-container").text(json.message);
                }
            }
        },
        onBeforeSend: function(jqXHR) {
            window.activeGridRequestXHR = jqXHR;
        },
        onLoadError: gridLoadError,
    });
    nodeGridDefinitions[nodeId] = definition;

    // Only one node grid query at a time - a user clicking through nodes would otherwise leave
    // several multi-minute queries running against each other
    if (window.activeGridRequestXHR) {
        window.activeGridRequestXHR.abort();
        window.activeGridRequestXHR = null;
    }

    return definition.setup().then(function(built) {
        if (!built) {
            return null;  // node errors - on_error_function has already put them on the page
        }
        const dataTable = built.dataTable;
        dataTable.on('draw.dt', function() {
            gridComplete();
            gridCompleteExtra();
        });
        if (autoLoad) {
            loadNodeGridData(nodeId, unique_code);
        } else {
            // The rows are being held back, but the table itself is up - say so, or the editor's
            // everythingLoaded (which waits on both halves) sits behind its overlay until the user
            // clicks "Show grid". gridComplete tells the two apart with nodeGridHasData()
            gridComplete();
            gridCompleteExtra();
        }
        return built;
    });
}

// Fire the row query for a node whose table was built but held back (@see setupNodeGrid autoLoad)
function loadNodeGridData(nodeId, unique_code) {
    const dataTable = getNodeDataTable(nodeId, unique_code);
    if (dataTable) {
        dataTable.ajax.reload();
    }
}

// True once the grid has actually fetched rows - a built-but-deferred table has made no request
function nodeGridHasData(nodeId, unique_code) {
    const dataTable = getNodeDataTable(nodeId, unique_code);
    return Boolean(dataTable && dataTable.ajax.json());
}

function gridLoadError(jqXHR, textStatus, errorThrown) {
    if (errorThrown == "abort") {
        console.log("task aborted...");
        return; // manually aborted
    }
    let errorMessage = errorThrown;
    const rj = jqXHR.responseJSON;
    if (rj) {
        if (rj.message) {
            errorMessage = rj.message;
        }
    }

    const ec = $("#error-container");
    ec.empty();
    ec.html("<ul class='messages'><li class='error'>Grid failed to load due to: " + errorMessage + "</li></ul>");
    $("#node-data-container").empty();
    hideLoadingOverlay();
    if (typeof hideGridLoadingOverlay === 'function') {
        hideGridLoadingOverlay();  // clear the grid-only overlay from a deferred "Show grid" load
    }
}

