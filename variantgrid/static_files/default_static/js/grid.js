function getGrid(nodeId, unique_code) {
	return $("#grid-" + nodeId, "#" + unique_code);
}

function export_grid(analysisId, nodeId, unique_code, export_type, use_canonical_transcripts) {
	let grid = getGrid(nodeId, unique_code);
	let gridParam = grid.jqGrid('getGridParam', 'postData');
	gridParam['rows'] = 0; // no pagination
	gridParam['export_type'] = export_type;
	gridParam['use_canonical_transcripts'] = use_canonical_transcripts;

	let querystring = EncodeQueryData(gridParam);
	let url = Urls.node_grid_export(analysisId) + "?" + querystring;
	window.location = url;
}

function load_variant_details(variant_id) {
    const aWin = getAnalysisWindow();
    let variant_details_url = Urls.view_variant(variant_id);
    if (aWin.ANALYSIS_SETTINGS) {
        if (aWin.ANALYSIS_SETTINGS.open_variant_details_in_new_window) {
            let newWin = aWin.open(variant_details_url, '_blank');
            newWin.focus();
            return;
        } else {
            let annotation_version_id = aWin.ANALYSIS_SETTINGS["annotation_version"];
            if (annotation_version_id) {
                variant_details_url = Urls.variant_details_annotation_version(variant_id, annotation_version_id);
            }
        }
    }
    let editorContainer = $("#node-editor-container");
    editorContainer.html('<div class="editor-loading"><img src="/static/images/spinner.gif" /> Loading variant details...</div>');
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
                for (let fromValue in replaceDict) {
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
    let igvUrl = createIgvUrl(locus, getBamsFuncString);
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

function isNodeVisible(options) {
    // Default to True so that any cached grid data won't be missing new field
    let nodeVisible = true;
    let analysisNode = options.colModel.analysisNode;
    if (analysisNode) {
        nodeVisible = analysisNode.visible;
    }
    return nodeVisible;
}


function detailsLink(variantId, options, rowData) {
    let nodeVisible = isNodeVisible(options);
    let kwargs = options.colModel.formatter_kwargs;
    if (kwargs) {
        nodeVisible = kwargs.node_visible;
    }

    let variantBoxes = [];
    if (nodeVisible) {
        let variant_selector = "<input type='checkbox' class='variant-select' variant_id=" + variantId + ">";
        variantBoxes.push(variant_selector);
    }

    let detailsUrl = "javascript:load_variant_details(" + variantId + ");";
    let detailsLink = createGridLink('View details', detailsUrl, '', ['variant-link'], ['view-details-link']);
    variantBoxes.push(detailsLink);

    // ClinVar
    let cvHighestPath = rowData["clinvar__highest_pathogenicity"];
    let cvClinSig = rowData["clinvar__clinical_significance"];

    let linkUrl = null;
    let extraLinkClasses = ['node-count-legend-C'];
    let extraIconClasses = [];
    let cvTitle = "ClinVar: ";
    if (cvHighestPath !== null) {
        cvTitle += cvClinSig;
        linkUrl = 'javascript:showGridCell("clinvar__clinical_significance")';
    } else {
        cvTitle += "not classified";
        extraIconClasses.push("no-entry");
        cvHighestPath = '';
    }
    let cvLink = createGridLink(cvTitle, linkUrl, cvHighestPath, extraLinkClasses, extraIconClasses);
    variantBoxes.push(cvLink);

    // Internally Classified
    let intMaxClass = rowData["max_internal_classification"];
    let intClassified = rowData["internally_classified"];

    linkUrl = null;
    let icTitle = "Internally Classified: ";
    extraLinkClasses = ['node-count-legend-G'];
    extraIconClasses = [];

    if (intMaxClass !== null) {
        icTitle += intClassified;
        linkUrl = 'javascript:showGridCell("max_internal_classification")';
    } else {
        icTitle += "not classified";
        extraIconClasses.push("no-entry");
        intMaxClass = '';
    }
    let icLink = createGridLink(icTitle, linkUrl, intMaxClass, extraLinkClasses, extraIconClasses);
    variantBoxes.push(icLink);

    let locus = rowData["locus__contig__name"] + ":" + rowData["locus__position"];
    let igvLink = create_igv_link(locus, 'getBams');
    variantBoxes.push(igvLink);
    return "<span class='variant_id-container' variant_id=" + variantId + ">" + variantBoxes.join('') + "</span>";
}

function clinvarLink(clinvar_variation_id) {
    let clinvar_string = '';
    if (clinvar_variation_id) {
        clinvar_string = "<a title='View ClinVar entry in new window' target='_blank' href='http://www.ncbi.nlm.nih.gov/clinvar/variation/" + clinvar_variation_id + "'>" + clinvar_variation_id + "</a>";
    }
    return clinvar_string;
}

function cosmicLink(cosmic_ids) {
    const COSMIC_PREFIX = "COSV";
    const COSMIC_LEGACY_PREFIX = "COSM";

    let cosmic_string = '';
    if (cosmic_ids) {
        const cosmic_ids_list = cosmic_ids.split("&");
        const cosmic_links = [];
        for(let i=0 ; i<cosmic_ids_list.length ; ++i) {
            let cosmic_id = cosmic_ids_list[i];
            if (cosmic_id.startsWith(COSMIC_PREFIX)) {
                // #2637 - COSMIC switched to using COSV in 2019, I can't find a direct link but search works then you select transcript
                cosmic_id = "<a title='View COSMIC entry in new window' target='_blank' href=' https://cancer.sanger.ac.uk/cosmic/search?q=" + cosmic_id + "'>" + cosmic_id + "</a>";
            } else if (cosmic_id.startsWith(COSMIC_LEGACY_PREFIX)) {
                const cosmic_id_int = cosmic_id.replace("COSM", "");
                cosmic_id = "<a title='View COSMIC entry in new window' target='_blank' href='http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=" + cosmic_id_int + "'>" + cosmic_id + "</a>";
            }
            cosmic_links.push(cosmic_id);
        }
        
        cosmic_string = cosmic_links.join();
    }
    return cosmic_string;
}


function omimLink(omim_id) {
    let omim_string = '';
    if (omim_id) {
        omim_string = "<a title='View OMIM entry in new window' target='_blank' href='https://www.omim.org/entry/" + omim_id + "'>" + omim_id + "</a>";
    }
    return omim_string;
}


function _geneSymbolLink(geneSymbolColumn, filterChildLink) {
    let columnString = '';
    if (geneSymbolColumn) {
        const geneSymbolList = geneSymbolColumn.split(",");
        const geneSymbolLinks = [];
        for(let i=0 ; i<geneSymbolList.length ; ++i) {
            let geneSymbol = geneSymbolList[i];
            let geneLinkString = '';
            if (filterChildLink) {
                let filterGeneLink = "javascript:createFilterChild(\"gene_symbol\", \"" + geneSymbol + "\");";
                geneLinkString = "<a class='grid-link' title='Filter to " + geneSymbol + "' href='" + filterGeneLink + "'><div class='grid-link-icon GeneListNode'></div></a>";
                geneLinkString += " <a class='left' target='_blank' title='View gene in new window' href='" + Urls.view_gene_symbol(geneSymbol) + "'>" + geneSymbol + "</a> ";
            } else {
                // not left
                geneLinkString += " <a target='_blank' title='View gene in new window' href='" + Urls.view_gene_symbol(geneSymbol) + "'>" + geneSymbol + "</a> ";
            }
            geneSymbolLinks.push(geneLinkString);
        }
        columnString = geneSymbolLinks.join();
    }
    return columnString;
}

function geneSymbolLink(geneSymbol, options) {
    let filterChildLink = isNodeVisible(options); // don't create kids for analysis wide tag nodes
    return _geneSymbolLink(geneSymbol, filterChildLink);
}

function geneSymbolNewWindowLink(geneSymbol) {
    return _geneSymbolLink(geneSymbol, false);
}


function showTagAutocomplete(variantId) {
    let container = $("#tag-entry-container-" + variantId);
    let addTagButton = $(".show-tag-autocomplete", container.parent());
    let nodeId = container.parents("#node-data-container").attr("node_id");
    
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


function getVariantTagHtml(variantId, tag, readOnly, tagLabel) {
    if (typeof(tagLabel) === 'undefined') {
        tagLabel = tag;
    }
    let outerClasses = ["grid-tag", "tagged-" + tag]
    if (!readOnly) {
        outerClasses.push("grid-tag-deletable");
    }
    return `<span class='${outerClasses.join(' ')}' title='Tagged as ${tag}' variant_id='${variantId}' tag_id='${tag}'><span class='user-tag-colored'>${tagLabel}</span></span>`;
}


// This is driven entirely off variantTags (not passed through SQL->JQGrid)
// This is so we can add/remove tags without wrecking cache  
function tagsFormatter(tagsCellValue, a, rowData) {
    let variantId = rowData['id'];
    let tagHtml = "";
    let aWin = getAnalysisWindow();

    if (!aWin.variantTagsReadOnly) {
        tagHtml += "<a class='show-tag-autocomplete' href='javascript:showTagAutocomplete(" + variantId + ")'><span class='add-variant-tag' title='Tag variant..'></span></a>";
        tagHtml += "<span id='tag-entry-container-" + variantId + "'></span>";
    }

    let tagList = aWin.variantTags[variantId];
    if (tagList) {
        for (let i=0 ; i<tagList.length ; ++i) {
            let tag = tagList[i];
            tagHtml += getVariantTagHtml(variantId, tag, aWin.variantTagsReadOnly);
        }        
    }
    return tagHtml;
}


function tagsGlobalFormatter(value, a, rowData) {
    if (!value) {
        return "";
    }
    let variantId = rowData['id'];
    let tags = value.split("|");
    let tagCounts = {};
    for (let i=0 ; i<tags.length ; ++i) {
        let tag = tags[i];
        let count = tagCounts[tag] || 0;
        tagCounts[tag] = count + 1;
    }

    let tagGlobalHtml = "";
    let sortedKeys = Object.keys(tagCounts).sort();
    for (let i=0 ; i<sortedKeys.length ; ++i) {
        let tag = sortedKeys[i];
        let tagCount = tagCounts[tag];
        let tagLabel = tag;
        if (tagCount > 1) {
            tagLabel = `${tag} x ${tagCount}`;
        }
        tagGlobalHtml += getVariantTagHtml(variantId, tag, true, tagLabel);
    }
    return tagGlobalHtml;
}


function classifyAndCloseButton(tagId) {
    let classifyUrl = Urls.create_classification_for_variant_tag(tagId);
    window.open(classifyUrl, "_blank");
    $(`#tag-button-${tagId}`).remove();
}

function formatVariantTagFirstColumn(variantString, options, rowObject) {
    const REQUIRES_CLASSIFICATION = "RequiresClassification";
    const variantURL = Urls.view_variant(rowObject["variant__id"]);
    let cellValue = "<a href='" + variantURL + "' target='_blank'>" + variantString + "</a>";
    if (rowObject.tag__id == REQUIRES_CLASSIFICATION) {
        let tagId = rowObject.id;
        cellValue += `<a id="tag-button-${tagId}" class="btn btn-primary new-classification-button" href="javascript:classifyAndCloseButton(${tagId})"><i class="fas fa-plus-circle"></i> New Classification</a>`;
    }
    return cellValue;
}


function gnomADVariant(rowData) {
    let chrom = rowData["locus__contig__name"];
    if (chrom.startsWith("chr")) {
        chrom = chrom.substr(3);
    }
    return [chrom, rowData["locus__position"], rowData["locus__ref__seq"], rowData["alt__seq"]].join("-");
}


function gnomadFilteredFormatter(gnomadFilteredCellValue, a, rowData) {
    let gnomadFilteredString = '';
    if (gnomadFilteredCellValue !== null) {
        const filterDiv = $("<div/>").addClass("gnomad-flag-label");
        if (gnomadFilteredCellValue) {
            filterDiv.addClass("gnomad-flagged");
            filterDiv.text("Fail");
        } else {
            filterDiv.text("Pass");
        }
        const gv = gnomADVariant(rowData);
        const dataset = ANALYSIS_SETTINGS["genome_build"] === 'GRCh38'? 'gnomad_r3' : 'gnomad_r2_1';
        const url = `http://gnomad.broadinstitute.org/variant/${gv}?dataset=${dataset}`;
        const gnomADLink = $("<a />").addClass("gnomad-link").attr({
            "href": url,
            "target": "_blank",
            "title": "View in gnomAD"
        });
        gnomADLink.append(filterDiv);
        gnomadFilteredString = gnomADLink.get(0).outerHTML;
    }
    return gnomadFilteredString;
}


function formatClinGenAlleleId(cellValue) {
    // warning: doesn't use settings.CLINGEN_ALLELE_REGISTRY_DOMAIN as static JS
    if (cellValue) {
        let ca_id = cellValue;
        let url = `http://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_caid?caid=${ca_id}`;
        cellValue = `<a href="${url}" target="_blank">${ca_id}</a>`;
    } else {
        cellValue = "";
    }
    return cellValue;
}


function splitAndLink(rawValue, split, buildLinkFunc) {
    let formattedValue = '';
    if (rawValue) {
        const raw_value_list = rawValue.split(split);
        const links = [];
        for(let i=0 ; i<raw_value_list.length ; ++i) {
            let value = raw_value_list[i];
            links.push(buildLinkFunc(value));
        }
        formattedValue = links.join();
    }
    return formattedValue;
}

function formatDBSNP(dbsnp_rs_ids) {
    function buildDBSNPLink(dbsnp_id) {
        return "<a title='View dbSNP in new window' target='_blank' href='https://www.ncbi.nlm.nih.gov/snp/" + dbsnp_id + "'>" + dbsnp_id + "</a>"
    }
    return splitAndLink(dbsnp_rs_ids, "&", buildDBSNPLink);
}


function formatPubMed(pubmed) {
    function buildPubMedLink(pubmed_id) {
        return "<a title='View PubMed article in new window' target='_blank' href='https://pubmed.ncbi.nlm.nih.gov/" + pubmed_id + "'>" + pubmed_id + "</a>"
    }
    return splitAndLink(pubmed, "&", buildPubMedLink);
}


function formatOntologyTerms(ontology_terms) {
    function buildOntologyTermLink(ontology_term) {
        let termSlug = ontology_term.split(" ")[0].replace(":", "_");
        let url = Urls.ontology_term(termSlug);
        return "<a title='View Ontology Term in new window' target='_blank' href='" + url + "'>" + ontology_term + "</a>"
    }
    return splitAndLink(ontology_terms, " | ", buildOntologyTermLink);
}


function unitAsPercentFormatter(unitValue) {
    let percentValue = "";
    // Allele Frequency missing data passed as "." to match VCF
    // Shows falsey values (eg 0.0) or '.' as blank
    if (unitValue && unitValue !== ".") {
        let percent = (100.0 * unitValue).toPrecision(3);
        percentValue = percent + "%";
    }
    return percentValue;
}


function formatMasterMindMMID3(value) {
    function buildMasterMindLink(mmid3) {
        return "<a title='View MasterMind in new window' target='_blank' href='https://mastermind.genomenon.com/detail?mutation=" + mmid3 + "'>" + mmid3 + "</a>"
    }
    return splitAndLink(value, "&", buildMasterMindLink);
}


jQuery.extend($.fn.fmatter , {
    'detailsLink' : detailsLink,
    'tagsFormatter' : tagsFormatter,
    'tagsGlobalFormatter' : tagsGlobalFormatter,
    'clinvarLink' : clinvarLink,
    'cosmicLink' : cosmicLink,
    'omimLink' : omimLink,
    'formatClinGenAlleleId': formatClinGenAlleleId,
    'formatDBSNP': formatDBSNP,
    'formatOntologyTerms': formatOntologyTerms,
    'formatPubMed': formatPubMed,
    'geneSymbolLink' : geneSymbolLink,
    'geneSymbolNewWindowLink' : geneSymbolNewWindowLink,
    'gnomadFilteredFormatter' : gnomadFilteredFormatter,
    'unitAsPercentFormatter' : unitAsPercentFormatter,
    'formatMasterMindMMID3': formatMasterMindMMID3,
});


// We need to do this, so that we don't send up a changing timestamp and thus never get cached
function deleteNdParam(postData) {
    const myPostData = $.extend({}, postData); // make a copy of the input parameter
    myPostData._filters = myPostData.filters;
    delete myPostData.nd;
    return myPostData;
}

// FIXME: Duplicated in jqgrid.html
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



function setupGrid(config_url, analysisId, nodeId, versionId, unique_code, gridComplete, gridLoadError, on_error_function) {
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
                                $("#node-editor-container").text(data.message)
                            }
                        }
                    }
                }
                // height: auto screws up on firefox
                if (typeof(data["height"]) === "undefined" || data["height"] === "auto") {
                    data["height"] = null;
                }

                // You can only have 1 active grid request
                data["loadBeforeSend"] = function(xhr) {
                    window.activeGridRequestXHR = xhr;
                }
                if (window.activeGridRequestXHR) {
                    window.activeGridRequestXHR.abort();
                    window.activeGridRequestXHR = null;
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

				setRowChangeCallbacks(grid, data["caption"])

                grid.jqGrid(
		            'navButtonAdd', pagerId, {
		            caption : "CSV",
		            buttonicon : "ui-icon-arrowthickstop-1-s",
		            onClickButton : function() {
		            	export_grid(analysisId, nodeId, unique_code, 'csv');
		            },
		            title : "Download as CSV",
		            cursor : "pointer"
		        });

                const aWin = getAnalysisWindow();
                if (aWin.ANALYSIS_SETTINGS && aWin.ANALYSIS_SETTINGS.canonical_transcript_collection) {
                    const ctc = aWin.ANALYSIS_SETTINGS.canonical_transcript_collection
                    grid.jqGrid(
                        'navButtonAdd', pagerId, {
                        caption : "Canonical transcript CSV",
                        buttonicon : "ui-icon-arrowthickstop-1-s",
                        onClickButton : function() {
                            export_grid(analysisId, nodeId, unique_code, 'csv', true);
                        },
                        title : "Download CSV using transcripts from " + ctc,
                        cursor : "pointer"
                    });
                }

		        grid.jqGrid(
		            'navButtonAdd', pagerId, {
		            caption : "VCF",
		            buttonicon : "ui-icon-arrowthickstop-1-s",
		            onClickButton : function() {
		            	export_grid(analysisId, nodeId, unique_code, 'vcf');
		            },
		            title : "Download as VCF",
		            cursor : "pointer"
	        	});
			}
	    });
	});
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
}

