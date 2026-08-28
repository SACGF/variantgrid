/* Client renderers for the variant grids (the AbstractVariantGrid family).

   DataTables renderer signature is (data, type, row, ctx), where ctx is added by DataTableDefinition:
     ctx.extra    grid wide metadata from the definition JSON (JqGrid.get_datatable_extra)
     ctx.kwargs   this column's renderKwargs (the colmodel's formatter_kwargs)

   Shared page helpers (createGridLink, IGV, tags, load_variant_details) live in grid.js.
   @see library/django_utils/jqgrid_datatable_adapter.py for the formatter -> renderer mapping */

const VariantGridFormat = (function() {
    "use strict";
    const VariantGridFormat = function() {};
    VariantGridFormat.prototype = {};
    return VariantGridFormat;
})();


// Filter children and variant selection only exist inside the analysis editor, and analysis-wide
// tag nodes (node.visible false) don't get them either
function _isNodeVisible(ctx) {
    if (!inAnalysis()) {
        return false;
    }
    // Default to true so that any cached grid data won't be missing new field
    const analysisNode = ctx && ctx.extra && ctx.extra.analysisNode;
    return analysisNode ? analysisNode.visible : true;
}


function _genomeBuildName(ctx) {
    if (ctx && ctx.extra && ctx.extra.genomeBuild) {
        return ctx.extra.genomeBuild;
    }
    return (typeof ANALYSIS_SETTINGS !== "undefined") ? ANALYSIS_SETTINGS["genome_build"] : null;
}


function splitAndLink(rawValue, split, buildLinkFunc) {
    let formattedValue = '';
    if (rawValue) {
        const raw_value_list = rawValue.split(split);
        const links = [];
        for(let i=0 ; i<raw_value_list.length ; ++i) {
            const value = raw_value_list[i];
            links.push(buildLinkFunc(value));
        }
        formattedValue = links.join();
    }
    return formattedValue;
}


// Keys are Classification.clinical_significance codes / somatic summary tiers - values give the box
// contents and tooltip label
const GERMLINE_CLASSIFICATION_BOXES = {
    '0': {display: 'O', label: 'Other'},
    '1': {display: '1', label: 'Benign'},
    '2': {display: '2', label: 'Likely Benign'},
    '3': {display: '3', label: 'VUS'},
    '4': {display: '4', label: 'Likely Pathogenic'},
    '5': {display: '5', label: 'Pathogenic'},
    'U': {display: '', label: 'Unclassified'},
};

const SOMATIC_CLASSIFICATION_BOXES = {
    'tier_1': {display: '1', label: 'Tier I'},
    'tier_1_or_2': {display: '1/2', label: 'Tier I/II'},
    'tier_2': {display: '2', label: 'Tier II'},
    'tier_3': {display: '3', label: 'Tier III'},
    'tier_4': {display: '4', label: 'Tier IV'},
    'U': {display: '', label: 'No tier'},
};

function internalClassificationBox(title, maxClassification, classifiedSummary, gridColumn, boxLookup, originCssClass) {
    let linkUrl = null;
    let contents = '';
    const extraIconClasses = [originCssClass];  // Standard germline/somatic colour, see global.scss
    if (maxClassification != null) {  // != also catches undefined (cached grid data w/o new fields)
        const box = boxLookup[maxClassification] || {display: maxClassification, label: maxClassification};
        // classifiedSummary is every record's value joined by '|' - show the pretty labels
        const summaryLabels = (classifiedSummary || '').split('|').map(function(cs) {
            const summaryBox = boxLookup[cs];
            return summaryBox ? summaryBox.label : cs;
        });
        title += ": " + summaryLabels.join('|');
        linkUrl = 'javascript:showGridCell("' + gridColumn + '")';
        contents = box.display;
        if (contents.length > 1) {  // e.g. Tier "1/2" - shrink to fit the 16px box
            extraIconClasses.push("grid-link-small-text");
        }
    } else {
        title += ": not classified";
        extraIconClasses.push("no-entry");
    }
    return createGridLink(title, linkUrl, contents, [], extraIconClasses);
}


VariantGridFormat.detailsLink = (variantId, type, rowData, ctx) => {
    let nodeVisible = _isNodeVisible(ctx);
    const kwargs = ctx && ctx.kwargs;
    if (kwargs) {
        nodeVisible = kwargs.node_visible;
    }

    const variantBoxes = [];
    if (nodeVisible) {
        const variant_selector = "<input type='checkbox' class='variant-select' variant_id=" + variantId + ">";
        variantBoxes.push(variant_selector);
    }

    const detailsUrl = "javascript:load_variant_details(" + variantId + ");";
    const detailsLink = createGridLink('View details', detailsUrl, '', ['variant-link'], ['view-details-link']);
    variantBoxes.push(detailsLink);

    // ClinVar
    let cvHighestPath = rowData["clinvar__highest_pathogenicity"];
    const cvClinSig = rowData["clinvar__clinical_significance"];

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
    const cvLink = createGridLink(cvTitle, linkUrl, cvHighestPath, extraLinkClasses, extraIconClasses);
    variantBoxes.push(cvLink);

    // Internally Classified - one box each for germline and somatic, coloured with the standard
    // clinical significance pill colours (.cs-* / .scs-* in global.scss)
    variantBoxes.push(internalClassificationBox(
        "Internally Classified (Germline)", rowData["max_internal_classification"],
        rowData["internally_classified"], "max_internal_classification",
        GERMLINE_CLASSIFICATION_BOXES, "allele-origin-G"));
    variantBoxes.push(internalClassificationBox(
        "Internally Classified (Somatic)", rowData["max_internal_somatic_classification"],
        rowData["internally_classified_somatic"], "max_internal_somatic_classification",
        SOMATIC_CLASSIFICATION_BOXES, "allele-origin-S"));

    const locus = rowData["locus__contig__name"] + ":" + rowData["locus__position"];
    const igvLink = create_igv_link(locus, 'getBams');
    variantBoxes.push(igvLink);
    return "<span class='variant_id-container' variant_id=" + variantId + ">" + variantBoxes.join('') + "</span>";
};


// This is driven entirely off variantTags (not passed through SQL->the grid)
// This is so we can add/remove tags without wrecking cache
VariantGridFormat.tags = (tagsCellValue, type, rowData) => {
    const variantId = rowData['id'];
    let tagHtml = "";
    const aWin = getAnalysisWindow();
    const readOnly = aWin.variantTagsReadOnly || !inAnalysis();

    if (!readOnly) {
        tagHtml += "<a class='show-tag-autocomplete' href='javascript:showTagAutocomplete(" + variantId + ")'><span class='add-variant-tag' title='Tag variant..'></span></a>";
        tagHtml += "<span id='tag-entry-container-" + variantId + "'></span>";
    }

    const tagList = (aWin.variantTags || {})[variantId];
    if (tagList) {
        const sortedTags = sortVariantTags(aWin, tagList);
        for (let i=0 ; i<sortedTags.length ; ++i) {
            const tag = sortedTags[i];
            tagHtml += getVariantTagHtml(variantId, tag, readOnly);
        }
    }
    return tagHtml;
};


VariantGridFormat.tagsGlobal = (value, type, rowData) => {
    if (!value) {
        return "";
    }
    const variantId = rowData['id'];
    const aWin = getAnalysisWindow();
    // In an analysis this is the analysis variant_tag_stale_days setting (null/undefined = staleness off)
    const staleDays = aWin.variantTagStaleDays;
    let staleCutoff = null;  // ISO date - payload dates compare lexically
    if (staleDays) {
        staleCutoff = new Date(Date.now() - staleDays * 86400 * 1000).toISOString().slice(0, 10);
    }

    // Entries are "tag:date" - see get_variantgrid_extra_annotate
    const tagStats = {};
    const entries = value.split("|");
    for (let i=0 ; i<entries.length ; ++i) {
        const entry = entries[i];
        const sep = entry.lastIndexOf(":");
        const tag = sep >= 0 ? entry.slice(0, sep) : entry;
        const date = sep >= 0 ? entry.slice(sep + 1) : null;
        let stats = tagStats[tag];
        if (!stats) {
            stats = {total: 0, fresh: 0, mostRecent: null};
            tagStats[tag] = stats;
        }
        stats.total += 1;
        if (date && (!staleCutoff || date >= staleCutoff)) {
            stats.fresh += 1;
        }
        if (date && (!stats.mostRecent || date > stats.mostRecent)) {
            stats.mostRecent = date;
        }
    }

    let tagGlobalHtml = "";
    const sortedKeys = sortVariantTags(aWin, Object.keys(tagStats));
    for (let i=0 ; i<sortedKeys.length ; ++i) {
        const tag = sortedKeys[i];
        const stats = tagStats[tag];
        let tagLabel = tag;
        if (stats.total > 1) {
            tagLabel = `${tag} x ${stats.total}`;
        }
        let title;
        const extraClasses = [];
        if (staleCutoff) {
            if (stats.fresh === 0) {  // Most recent event is older than the cutoff
                extraClasses.push("grid-tag-stale");
                tagLabel += " <i class='fas fa-clock'></i>";
                title = `Tagged as ${tag} - no events within the last ${staleDays} days, most recent ${stats.mostRecent}`;
            } else if (stats.total > 1) {
                tagLabel = `${tag} x ${stats.total} (${stats.fresh} fresh)`;
                title = `${stats.fresh} of ${stats.total} tag events within the last ${staleDays} days, most recent ${stats.mostRecent}`;
            }
        }
        tagGlobalHtml += getVariantTagHtml(variantId, tag, true, tagLabel, extraClasses, title);
    }
    return tagGlobalHtml;
};


VariantGridFormat.clinvarLink = (clinvar_variation_id) => {
    let clinvar_string = '';
    if (clinvar_variation_id) {
        clinvar_string = "<a title='View ClinVar entry in new window' target='_blank' href='http://www.ncbi.nlm.nih.gov/clinvar/variation/" + clinvar_variation_id + "'>" + clinvar_variation_id + "</a>";
    }
    return clinvar_string;
};


VariantGridFormat.cosmicLink = (cosmic_ids) => {
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
};


VariantGridFormat.omimLink = (omim_id) => {
    let omim_string = '';
    if (omim_id) {
        omim_string = "<a title='View OMIM entry in new window' target='_blank' href='https://www.omim.org/entry/" + omim_id + "'>" + omim_id + "</a>";
    }
    return omim_string;
};


function _geneSymbolLink(geneSymbolColumn, filterChildLink) {
    let columnString = '';
    if (geneSymbolColumn) {
        const geneSymbolList = geneSymbolColumn.split(",");
        const geneSymbolLinks = [];
        for(let i=0 ; i<geneSymbolList.length ; ++i) {
            const geneSymbol = geneSymbolList[i];
            let geneLinkString = '';
            if (filterChildLink) {
                const filterGeneLink = "javascript:createFilterChild(\"gene_symbol\", \"" + geneSymbol + "\");";
                geneLinkString = "<a class='grid-link' title='Filter to " + geneSymbol + "' href='" + filterGeneLink + "'><i class='fa-solid fa-list-check gene-list-node-icon'></i></a>";
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


VariantGridFormat.geneSymbolLink = (geneSymbol, type, row, ctx) => {
    return _geneSymbolLink(geneSymbol, _isNodeVisible(ctx));
};


VariantGridFormat.geneSymbolNewWindowLink = (geneSymbol) => {
    return _geneSymbolLink(geneSymbol, false);
};


function mavedbUrnToUrls(mavedbUrn) {
  const experimentSets = new Set();
  const EXPERIMENT_URL_PREFIX = "https://www.mavedb.org/#/experiment-sets/";
  if (mavedbUrn) {
    mavedbUrn.split("&").forEach(urn => {
      const es = urn.rsplit("-", 2)[0];
      experimentSets.add(es);
    });
  }
  const urls = {};
  [...experimentSets].sort().forEach(urn => {
    urls[urn] = EXPERIMENT_URL_PREFIX + urn;
  });
  return urls;
}

// JavaScript does not have a built-in rsplit function, so we need to implement it.
String.prototype.rsplit = function(sep, maxsplit) {
  const split = this.split(sep);
  return maxsplit ? [split.slice(0, -maxsplit).join(sep), ...split.slice(-maxsplit)] : split;
};


VariantGridFormat.mavedbUrn = (mavedbUrn) => {
    const urls = mavedbUrnToUrls(mavedbUrn);
    return Object.entries(urls).map(([key, value]) => `<a target="_blank" href="${value}">${key}</a>`).join(" ");
};


function gnomADVariant(rowData) {
    let chrom = rowData["locus__contig__name"];
    if (chrom.startsWith("chr")) {
        chrom = chrom.substr(3);
    }
    return [chrom, rowData["locus__position"], rowData["locus__ref__seq"], rowData["alt__seq"]].join("-");
}


VariantGridFormat.gnomadFiltered = (gnomadFilteredCellValue, type, rowData, ctx) => {
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
        const dataset = _genomeBuildName(ctx) === 'GRCh38' ? 'gnomad_r3' : 'gnomad_r2_1';
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
};


VariantGridFormat.clinGenAlleleId = (cellValue) => {
    // warning: doesn't use settings.CLINGEN_ALLELE_REGISTRY_DOMAIN as static JS
    if (!cellValue) {
        return "";
    }
    const url = `http://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_caid?caid=${cellValue}`;
    return `<a href="${url}" target="_blank">${cellValue}</a>`;
};


VariantGridFormat.dbsnp = (dbsnp_rs_ids) => {
    function buildDBSNPLink(dbsnp_id) {
        return "<a title='View dbSNP in new window' target='_blank' href='https://www.ncbi.nlm.nih.gov/snp/" + dbsnp_id + "'>" + dbsnp_id + "</a>";
    }
    return splitAndLink(dbsnp_rs_ids, "&", buildDBSNPLink);
};


VariantGridFormat.pubMed = (pubmed) => {
    function buildPubMedLink(pubmed_id) {
        return "<a title='View PubMed article in new window' target='_blank' href='https://pubmed.ncbi.nlm.nih.gov/" + pubmed_id + "'>" + pubmed_id + "</a>";
    }
    return splitAndLink(pubmed, "&", buildPubMedLink);
};


VariantGridFormat.ontologyTerms = (ontology_terms) => {
    function buildOntologyTermLink(ontology_term) {
        const termSlug = ontology_term.split(" ")[0].replace(":", "_");
        const url = Urls.ontology_term(termSlug);
        return "<a title='View Ontology Term in new window' target='_blank' href='" + url + "'>" + ontology_term + "</a>";
    }
    return splitAndLink(ontology_terms, " | ", buildOntologyTermLink);
};


VariantGridFormat.masterMind = (value) => {
    function buildMasterMindLink(mmid3) {
        return "<a title='View MasterMind in new window' target='_blank' href='https://mastermind.genomenon.com/detail?mutation=" + mmid3 + "'>" + mmid3 + "</a>";
    }
    return splitAndLink(value, "&", buildMasterMindLink);
};


VariantGridFormat.unitAsPercent = (unitValue) => {
    // Allele Frequency missing data passed as "." to match VCF
    // Shows falsey values (eg 0.0) or '.' as blank
    if (!unitValue || unitValue === ".") {
        return "";
    }
    return (100.0 * unitValue).toPrecision(3) + "%";
};


// renderKwargs: url_name, url_object_column, icon_css_class
VariantGridFormat.link = (cellValue, type, row, ctx) => {
    const kwargs = (ctx && ctx.kwargs) || {};
    const cssClasses = ["icon24", "left", "margin-r-5"];
    const iconList = kwargs.icon_css_class ? [kwargs.icon_css_class] : ["view-details-link"];
    const icons = iconList.map(icon => `<div class='${cssClasses.concat(icon).join(" ")}'></div>`).join("");
    const urlObject = kwargs.url_object_column ? row[kwargs.url_object_column] : cellValue;
    const url = Urls[kwargs.url_name](urlObject);
    return `<a class='grid-link' href='${url}'>${icons}<div class='display-text'>${cellValue}</div></a>`;
};
