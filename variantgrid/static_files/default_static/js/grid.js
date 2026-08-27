// LEGACY-JQGRID: node grid element lookup - callers move to $(...).DataTable() in plan Stage 4
function getGrid(nodeId, unique_code) {
	return $("#grid-" + nodeId, "#" + unique_code);
}

function getAnalysisDownloadTracker() {
	return getAnalysisWindow().analysisDownloadTracker;
}

/* Everything a node export is identified by, plus the URL that launches it. Built off the grid's live
   postData so the export sees whatever the user has filtered the grid down to. */
function nodeGridExportInfo(analysisId, nodeId, unique_code, export_type, use_canonical_transcripts, caption) {
	const grid = getGrid(nodeId, unique_code);
	const gridParams = $.extend({}, grid.jqGrid('getGridParam', 'postData'));
	gridParams['rows'] = 0; // no pagination
	gridParams['export_type'] = export_type;
	if (use_canonical_transcripts) {
		gridParams['use_canonical_transcripts'] = true;
	}

	const gridCaption = grid.jqGrid('getGridParam', 'caption') || ("Node " + nodeId);
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
    const selector = $("td[aria-describedby*='" + gridColumn + "']");
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


// LEGACY-JQGRID: the pages still on {% jqgrid %} register the renderers under their old formatter
// names. jqGrid calls (cellvalue, options, rowObject), so map options.colModel onto the context
// object the DataTables renderers take. Delete with the last {% jqgrid %} variant grid.
function jqGridFormatterContext(options) {
    const colModel = (options && options.colModel) || {};
    return {extra: {analysisNode: colModel.analysisNode}, kwargs: colModel.formatter_kwargs};
}

jQuery.extend($.fn.fmatter , {
    'detailsLink': (v, o, row) => VariantGridFormat.detailsLink(v, null, row, jqGridFormatterContext(o)),
    'geneSymbolLink': (v, o, row) => VariantGridFormat.geneSymbolLink(v, null, row, jqGridFormatterContext(o)),
    'gnomadFilteredFormatter': (v, o, row) => VariantGridFormat.gnomadFiltered(v, null, row, jqGridFormatterContext(o)),
    'tagsFormatter': (v, o, row) => VariantGridFormat.tags(v, null, row),
    'tagsGlobalFormatter': (v, o, row) => VariantGridFormat.tagsGlobal(v, null, row),
    'clinvarLink': (v) => VariantGridFormat.clinvarLink(v),
    'cosmicLink': (v) => VariantGridFormat.cosmicLink(v),
    'omimLink': (v) => VariantGridFormat.omimLink(v),
    'formatClinGenAlleleId': (v) => VariantGridFormat.clinGenAlleleId(v),
    'formatDBSNP': (v) => VariantGridFormat.dbsnp(v),
    'formatOntologyTerms': (v) => VariantGridFormat.ontologyTerms(v),
    'formatPubMed': (v) => VariantGridFormat.pubMed(v),
    'geneSymbolNewWindowLink': (v) => VariantGridFormat.geneSymbolNewWindowLink(v),
    'unitAsPercentFormatter': (v) => VariantGridFormat.unitAsPercent(v),
    'formatMasterMindMMID3': (v) => VariantGridFormat.masterMind(v),
    'formatMavedbUrnLinks': (v) => VariantGridFormat.mavedbUrn(v),
});


// LEGACY-JQGRID: used by setupGrid below.
// We need to do this, so that we don't send up a changing timestamp and thus never get cached
function deleteNdParam(postData) {
    const myPostData = $.extend({}, postData); // make a copy of the input parameter
    myPostData._filters = myPostData.filters;
    delete myPostData.nd;
    return myPostData;
}

// LEGACY-JQGRID: DataTables persists rows per page through the definition's gridName
// (@see DataTableDefinition.setupDom). FIXME: Duplicated in jqgrid.html
function setRowChangeCallbacks(grid, gridName) {
	$(".ui-pg-selbox").change(function() {
        const gridRows = $(this).val();
        const data = 'grid_name=' + gridName + '&grid_rows=' + gridRows;
        $.ajax({
		    type: "POST",
		    data: data,
		    url: Urls.set_user_row_config(),
		});
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


// This is always kicked off after grid is loaded (after passed in function gridComplete)
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



// LEGACY-JQGRID: the analysis node grid's jqGrid setup - replaced by a DataTableDefinition build in
// node_data_grid.html when the node grid converts (plan Stage 4)
function setupGrid(config_url, analysisId, nodeId, versionId, unique_code, gridComplete, gridLoadError, on_error_function, autoLoad) {
    if (typeof autoLoad === "undefined") { autoLoad = true; }
	$(function () {
    	$.getJSON(config_url, function(data) {
            const errors = data["errors"];
            if (errors) {
				on_error_function(errors);
    		} else {
                const postData = data["postData"] || {};
                // TODO: From issue #1041 6/6/2018 - remove this when nodes config cache expires in 1 week.
                if (typeof postData["node_id"] == "undefined") {
                    postData["node_id"] = nodeId;
                }
				// end obsolete code... 
				
				data["postData"] = postData;
				data["serializeGridData"] = deleteNdParam;
				data["shrinkToFit"] = false;

                const pagerId = '#pager-' + nodeId;
                data["pager"] = pagerId;
				data["gridComplete"] = function() {
				    gridComplete();
				    gridCompleteExtra();
				};
				data["loadError"] = gridLoadError;
                data["loadComplete"] = function(data) {
                    if (data.non_fatal) {
                        // console.log("jQGrid: Non fatal node error...");
                        // console.log(data);
                        if (data.deleted_nodes) {
                            deleteNodesFromDOM(data.deleted_nodes, []);
                            if (data.message) {
                                $("#node-editor-container").text(data.message);
                            }
                        }
                    }
                };
                // height: auto screws up on firefox
                if (typeof(data["height"]) === "undefined" || data["height"] === "auto") {
                    data["height"] = null;
                }

                // You can only have 1 active grid request
                data["loadBeforeSend"] = function(xhr) {
                    window.activeGridRequestXHR = xhr;
                };
                if (window.activeGridRequestXHR) {
                    window.activeGridRequestXHR.abort();
                    window.activeGridRequestXHR = null;
                }

                // Remember the server datatype so a deferred load can flip back to it.
                window.nodeGridServerDatatype = window.nodeGridServerDatatype || {};
                window.nodeGridServerDatatype[nodeId] = data["datatype"] || "json";
                if (!autoLoad) {
                    data["datatype"] = "local";  // build colModel/pager/nav, fetch no rows
                }

                const grid = getGrid(nodeId, unique_code);
                grid.jqGrid(data).navGrid(pagerId,
	                	{add: false, edit: false, del: false, view: false, search:false},
			       		{}, // edit options
			        	{}, // add options
			       	 	{}, // del options 
			        	{ multipleSearch:true, closeOnEscape:true }, // search options 
			        	{} // view options 
		        	);

				setRowChangeCallbacks(grid, data["caption"]);

                const csvButtonId = `node-grid-export-csv-${nodeId}`;
                grid.jqGrid(
		            'navButtonAdd', pagerId, {
		            id : csvButtonId,
		            caption : "CSV",
		            buttonicon : "ui-icon-arrowthickstop-1-s",
		            onClickButton : function() {
		            	export_grid(analysisId, nodeId, unique_code, 'csv');
		            },
		            title : "Download as CSV",
		            cursor : "pointer"
		        });
                registerNodeGridDownloadButton(`#${csvButtonId}`, analysisId, nodeId, unique_code, 'csv',
                                               false, "CSV");

                const aWin = getAnalysisWindow();
                if (aWin.ANALYSIS_SETTINGS && aWin.ANALYSIS_SETTINGS.canonical_transcript_collection) {
                    const ctc = aWin.ANALYSIS_SETTINGS.canonical_transcript_collection;
                    const ctcButtonId = `node-grid-export-canonical-csv-${nodeId}`;
                    grid.jqGrid(
                        'navButtonAdd', pagerId, {
                        id : ctcButtonId,
                        caption : "Canonical transcript CSV",
                        buttonicon : "ui-icon-arrowthickstop-1-s",
                        onClickButton : function() {
                            export_grid(analysisId, nodeId, unique_code, 'csv', true);
                        },
                        title : "Download CSV using transcripts from " + ctc,
                        cursor : "pointer"
                    });
                    registerNodeGridDownloadButton(`#${ctcButtonId}`, analysisId, nodeId, unique_code, 'csv',
                                                   true, "Canonical transcript CSV");
                }

                const vcfButtonId = `node-grid-export-vcf-${nodeId}`;
		        grid.jqGrid(
		            'navButtonAdd', pagerId, {
		            id : vcfButtonId,
		            caption : "VCF",
		            buttonicon : "ui-icon-arrowthickstop-1-s",
		            onClickButton : function() {
		            	export_grid(analysisId, nodeId, unique_code, 'vcf');
		            },
		            title : "Download as VCF",
		            cursor : "pointer"
	        	});
                registerNodeGridDownloadButton(`#${vcfButtonId}`, analysisId, nodeId, unique_code, 'vcf',
                                               false, "VCF");
			}
	    });
	});
}

// Fire the deferred row query (phase 2) for a node whose grid was config-loaded only.
function loadNodeGridData(nodeId, unique_code) {
    const grid = getGrid(nodeId, unique_code);
    const datatype = (window.nodeGridServerDatatype || {})[nodeId] || "json";
    grid.jqGrid('setGridParam', {datatype: datatype}).trigger('reloadGrid');
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

