// @ts-check
// VCTable: DataTables client renderers for classification values (hgvs, condition, classification, identifier...),
// named as client_renderer strings by grids shown site-wide, so base.html loads this on every page.

const VCTable = (function() {
    const VCTable = function() {};
    VCTable.prototype = {};
    return VCTable;
})();

VCTable.format_condition = function(condition_json) {
    if (!condition_json) {
        return $('<span>', {text: "-", class:'no-value'});
    }
    const MAX_TERM_NAME_LENGTH = 50;
    const termDoms = [];
    if (condition_json.resolved_terms) {
        for (const term of condition_json.resolved_terms) {
            termDoms.push($('<div>', {
                class: 'ontology-term semicolon-sep',
                html: [
                    $('<a>', {
                        class: 'hover-link',
                        html: $('<span>', {class: 'term-id', text: term.term_id}),
                        href: Urls.ontology_term(term.term_id.replace(':', '_'))
                    }),
                    " ",
                    limitLengthSpan(term.name, MAX_TERM_NAME_LENGTH).addClass('term-name')
                ]
            }));
        }
    }
    if (condition_json.plain_text_terms) {
        for (const term of condition_json.plain_text_terms) {
            termDoms.push($('<div>', {class:'ontology-term free-text semicolon-sep', html: limitLengthSpan(term, MAX_TERM_NAME_LENGTH)}));
        }
    }

    if (!termDoms.length) {
        return $('<div>', {class: 'ontology-term free-text', html: limitLengthSpan(condition_json.display_text, MAX_TERM_NAME_LENGTH)});
    }

    const dom = $('<div>');
    const MAX_VISIBLE_CONDITIONS = 3;
    if (termDoms.length > MAX_VISIBLE_CONDITIONS) {
        dom.append(termDoms.slice(0, MAX_VISIBLE_CONDITIONS));
        $('<a>', {
            class: 'show-more-conditions hover-link',
            href: '#',
            text: `+ ${termDoms.length - MAX_VISIBLE_CONDITIONS} more conditions`,
            title: 'click to show all conditions'
        }).appendTo(dom);
        $('<div>', {class: 'hidden-conditions', style: 'display:none'}).append(termDoms.slice(MAX_VISIBLE_CONDITIONS)).appendTo(dom);
    } else {
        dom.append(termDoms);
    }
    if (condition_json.resolved_terms && condition_json.resolved_terms.length > 1 && condition_json.resolved_join) {
        $('<div>', {class: 'font-italic', text:condition_json.resolved_join === 'C' ? ' Co-occurring' : ' Uncertain'}).appendTo(dom);
    }
    return dom;
};

// grid cells are rendered as static HTML, so the show-more toggle needs a delegated handler
$(document).on('click', '.show-more-conditions', function() {
    $(this).hide().siblings('.hidden-conditions').show();
    return false;
});

VCTable.format_hgvs = (parts) => {
    if (typeof(parts) == 'string') {
        parts = parts.trim();
        if (!parts.length) {
            return $('<span>', {text:'-', class:'no-value'});
        }
        return $('<span>', {text:limitLength(parts, 100)});
    } else if (!parts) {
        return $('<span>', {text:'-', class:'no-value'});
    }
    const genomeBuild = parts.genome_build;
    const transcript = parts.transcript;
    const geneSymbol = parts.gene_symbol;
    const cNomen = parts.c_nomen || parts.variant; // older code called cNomen 'variant'
    const allele = parts.allele;
    const variantId = parts.variant_id;
    const alleleId = parts.allele_id;
    const alleleInfoId = parts.allele_info_id;
    // let validationInclude = parts.validation_include;
    // let alleleInfoStatus = parts.allele_info_status;
    const icon = parts.icon;
    const tooltip = parts.tooltip;
    const pHgvs = parts.p_hgvs;
    let url = null;
    const error = parts.error;

    if (error) {
        return $('<span>', {html: [
            $('<i>', {class: 'fa-solid fa-circle-exclamation text-danger'}),
            document.createTextNode(` ${error}`)
        ]});
    }

    if (variantId) {
        url = Urls.view_allele_from_variant(variantId);
    } else if (alleleId) {
        url = Urls.view_allele(alleleId);
    }
    // also turn into a link

    const outterDom = $('<div>').css("display", "inline-block");
    const dom = $('<div>').appendTo(outterDom);

    if (allele) {
        outterDom.prepend($('<dom>', {class: 'font-weight-bold', text: allele}));
    }

    if (genomeBuild && (parts.desired === false || parts.normalized === false || parts.resolved_without_hgvs || parts.always_show_genome_build)) {
        const genomeBuildWrapper = $('<div>');
        let addNewLine = false;
        if (parts.normalized === false) {
            $('<span>', {html: 'not resolved<br/>showing imported ', style:'color:#888'}).appendTo(genomeBuildWrapper);
        } else if (parts.resolved_without_hgvs) {
            $('<span>', {html: 'resolved, no c.HGVS<br/>showing imported ', style:'color:#888'}).appendTo(genomeBuildWrapper);
        } else if (parts.desired === false) {
            $('<span>', {html: 'not lifted-over ', style:'color:#888'}).appendTo(genomeBuildWrapper);
        } else {
            addNewLine = true;
        }

        dom.append(genomeBuildWrapper);
        dom.append($('<span>', {text: genomeBuild, style:'font-weight:500; color:#888'}));
        if (addNewLine) {
            dom.append("<br/>");
        } else {
            dom.append(' ');
        }
        // <span style="white-space: nowrap"><span>{{ c_hgvs.transcript }}</span>{% if c_hgvs.gene_symbol %}(<span class="text-secondary" style="letter-spacing: 0.5px">{{ c_hgvs.gene_symbol }}</span>){% endif %}:</span><span style="display:inline-block;word-break: break-all">{{ c_hgvs.raw_c }}</span>
    }

    let cDom = $('<span>', {class:'c-hgvs-body'});
    if (transcript && cNomen) {
        cDom.append($('<span>', {text: transcript, class:'c-hgvs-transcript'}));
        if (geneSymbol) {
            cDom.append($('<span>', {class: 'c-hgvs-gene-symbol-b', html:[
                $('<span>', {text: "(", class: 'bracket'}),
                $('<span>', {class: 'c-hgvs-gene-symbol', text: geneSymbol}),
                $('<span>', {text: ")", class: 'bracket bracket-close'})
            ]}));
        }
        cDom.append($('<span>', {text: ":", class: 'colon'}));
        // used to be display:inline-block; but that doesn't underline
        cDom.append($('<span>', {class:'c-hgvs-nomen', text: limitLength(cNomen, 100)}));
    } else {
        cDom.append(document.createTextNode(limitLength(parts.full || '', 100)));
    }
    if (url) {
        cDom = $('<a>', {href: url, html: cDom, class: 'hover-link'});
    }
    dom.append(cDom);
    if (pHgvs) {
        dom.append($('<span>', {class: 'p-hgvs', text: limitLength(pHgvs), title: parts.p_hgvs_full || pHgvs}));
    }
    if (icon) {
        const iconData = {"class": icon};
        if (tooltip) {
            iconData.title = tooltip;
        }
        const iconDom = $('<i>', iconData);
        iconDom.appendTo(dom);
    }

    return outterDom;
};

VCTable.hgvs = (data, type, row) => {
    return VCTable.format_hgvs(data);
};

VCTable.condition = (data, type, row) => {
    return VCTable.format_condition(data);
};

VCTable.latest_curation_and_link = (data, type, row) => {
    const lastCurated = data["curation_date"];
    if (lastCurated) {
        const dom = $("<div>");
        let content = $("<data>", {"class": "convert-timestamp time-ago", "data-date": lastCurated});

        const classificationId = data["classification_id"];
        if (classificationId) {
             content = $('<a>', {
                href: Urls.view_classification(classificationId),
                class: 'hover-link',
                html: content
            });
            dom.append(content);
        }
        const dateType = data["date_type"];
        if (dateType) {
            dom.append($("<div>", {text: dateType}));
        }
        return dom;
    } else {
        return $("<span>", {text:"-", class: "no-value"});
    }
};

VCTable.somatic_clinical_significance = (data, type, row) => {
    let value = null;
    if (data === null) {
        return "";
    }
    if (typeof(data) === "string") {
        const parts = data.split("|");
        value = {};
        value[SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE] = parts[0];
        if (parts.length > 1) {
            value["amp_level"] = parts[1];
        }
    } else {
        value = data;
    }
    const pendingVal = value["pending"];
    const scs = value["clinical_significance"];

    if (scs) {
        const scsKey = EKeys.cachedKeys.key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE);
        const scsLabel = scsKey.prettyValue(scs);
        const pendingClass = pendingVal ? 'strike' : '';
        const dom = ($('<span>', {text: scsLabel.val, 'class': `c-pill cs cs-${scs} ${pendingClass}`}));

        const highest_level = value["amp_level"];
        if (highest_level) {
            dom.append(`<span class="amp-level">${highest_level}</span>`);
        }
        if (value["diff"]) {
            dom.append('<i class="fa-solid fa-asterisk ml-1" title="Multiple values have been recorded - showing latest"></i>');
        }
        if (pendingVal) {
            const domParts = [dom];
            if (pendingVal == 'in-review') {
                domParts.push($('<div>', {
                    class: `cs c-pill cs-in-review`,
                    html: "In-Review"
                }));
            } else {
                domParts.push($('<div>', {
                    class: `cs c-pill cs-${pendingVal.toLowerCase()}`,
                    html: scsKey.prettyValue(pendingVal).val
                }));
            }
            return $('<div>', {html: domParts});
        }

        return dom;
    } else {
        // FIXME can still have different values with no-value
        return $('<div>', {class: 'c-pill cs-none no-value', html: 'No Data'});
    }
};

VCTable.classification = (data, type, row) => {
    // as in render the "classification" known as clinical_significance value
    if (data === null) {
        return ""; // support for dirty groups still processing
    }
    let value = data;
    if (typeof(data) == "string") {
        value = {"classification": data};
    }

    const csKey = EKeys.cachedKeys.key(SpecialEKeys.CLINICAL_SIGNIFICANCE);
    const csVal = value["classification"] || value[SpecialEKeys.CLINICAL_SIGNIFICANCE];
    const pendingVal = value["pending"];

    const domParts = [];
    let dom;

    if (!csVal || !csVal.length) {
        dom = $('<div>', {class: 'c-pill cs-none no-value', html: 'No Data'});
    } else {
        const pendingClass = pendingVal ? 'strike' : '';
        dom = $('<div>', {class: `cs c-pill cs-${csVal.toLowerCase()} ${pendingClass}`, text: csKey.prettyValue(csVal).val});
    }
    domParts.push(dom);
    if (value["diff"]) {
        dom.append(' <i class="fa-solid fa-asterisk ml-1" title="Multiple values have been recorded - showing latest"></i>');
    }
    if (pendingVal) {
        if (pendingVal == 'in-review') {
            domParts.push($('<div>', {
                class: `cs c-pill cs-in-review`,
                html: "In-Review"
            }));
        } else {
            domParts.push($('<div>', {
                class: `cs c-pill cs-${pendingVal.toLowerCase()}`,
                html: csKey.prettyValue(pendingVal).val
            }));
        }
    }
    return $('<div>', {html: domParts});
};

VCTable.evidence_key = (key_name, data, type, row) => {
    const csKey = EKeys.cachedKeys.key(key_name);
    let span;
    if (data == null) {
        span = $('<span>', {class: 'no-value', text: '-'});
    } else {
        span = $('<span>');
        csKey.formatValue(data, span);
    }
    return span;
};

VCTable.allele_origin_bucket_label = (allele_origin_bucket, override_text = "", alignment="vertical") => {
    let allele_origin_label = override_text;
    if (!allele_origin_label) {
        if (allele_origin_bucket == "S") {
            allele_origin_label = "SOMATIC";
        } else if (allele_origin_bucket == "G") {
            allele_origin_label = "GERMLINE";
        } else if (allele_origin_bucket == "U") {
            allele_origin_label = "UNKNOWN";
        } else {
            allele_origin_label = "???";
        }
    }
    return $('<div>', {
        class: `allele-origin-box ${alignment} allele-origin-${allele_origin_bucket}`,
        html: [
            $('<div>', {
                class: 'allele-origin-text',
                text: allele_origin_label
            })
        ]
    });
};

VCTable.alleleGroupingIdentifier = (data, type, row) => {
    const dom = $("div");
    for (const allele of data.alleles) {
        dom.append(VCTable.format_hgvs(allele));
    }
    return dom;
};

// The grouping grid's expand URL (TableFormat.expandAjax looks it up on window): while the grid is filtered to a
// sample, the expanded records lead with that sample's
function classificationGroupingDetailUrl(groupingId) {
    const url = Urls.classification_grouping_detail(groupingId);
    const sampleId = $('#id_sample').val();
    return sampleId ? `${url}?sample=${sampleId}` : url;
}

VCTable.groupIdentifier = (data, type, row) => {
    const id = data.id;
    const dirty = data.dirty;
    const org_name = data.org_name;
    const lab_name = data.lab_name;
    const research = data.research;
    const research_icon = data.research_icon;
    const shareLevel = data.share_level;
    const allele_origin_bucket = data.allele_origin_bucket;

    const alleleOriginDiv = VCTable.allele_origin_bucket_label(allele_origin_bucket);

    const shareInfo = EKeys.shareLevelInfo(shareLevel);
    const icon = $('<img>', {src: shareInfo.icon, class:'share-icon'});

    const classification_count = data.classification_count;

    const labLine = [
        icon,
        $('<span>', {text: `${org_name} / ${lab_name}`})
    ];
    if (research) {
        labLine.push($('<span>', {class: 'badge badge-info ml-1', title: 'Research lab', text: `${research_icon} Research`}));
    }
    const dom = $('<div>', {html: labLine});

    if (allele_origin_bucket === "S") {
        dom.append($('<div>', {class:'testing-context', text: data.testing_context_bucket_label}));
    }

    if (dirty) {
        dom.append($("<div class='mt-2'><i class=\"fa-solid fa-clock\"></i> Data is currently being updated</div>"));
    } else if (classification_count === 0) {
        dom.append("-Invalid Record - no Classifications");
    } else if (data.sample_record_count) {
        // A grouping holds many samples' records - say how many of them made the sample filter match
        dom.append($('<div>', {class:'text-muted text-small', text: `${data.sample_record_count} of ${classification_count} records from this sample`}));
    } else if (classification_count > 1) {
        dom.append($('<div>', {class:'text-muted text-small', text: `${classification_count} records`}));
    }

    // matches
    if (data.matches) {
        const searchTerm = data.search;
        const match_doms = [];
        const ekeys = EKeys.cachedKeys;
        for (const [key, value] of Object.entries(data.matches)) {
            match_doms.push($('<div>', {class: 'search-result', html:[
                $('<span>', { text: ekeys.key(key).label + ' : '}),
                    highlightTextAsDom(searchTerm, value)
                ]
            }));
        }
        dom.append($('<div>', {html: match_doms}));
    }

    const indicatorClassName = `allele-origin-indicator allele-origin-horizontal allele-origin-${allele_origin_bucket}`;
    const fullDom = $('<div>', {style: 'margin-left: 5px; margin-top: -10px; display:flex; flex-direction:row', html:[
        alleleOriginDiv,
        dom
    ]});

    // if (data.share_level == "lab" || data.share_level == "institution") {
    //     fullDom.css({"opacity": 0.5});
    // }

    return fullDom;
};

VCTable.identifier = (data, type, row) => {
    const id = data.id;
    const org_name = data.org_name;
    const lab_name = data.lab_name;
    const lab_record_id = data.lab_record_id;
    const shareLevel = data.share_level;
    const allele_origin_bucket = data.allele_origin_bucket;

    const alleleOriginDiv = VCTable.allele_origin_bucket_label(allele_origin_bucket);

    const shareInfo = EKeys.shareLevelInfo(shareLevel);
    const icon = $('<img>', {src: shareInfo.icon, class:'share-icon'});

    const classification_count = data.classification_count;

    const content = [
        icon,
        $('<span>', {text: `${org_name} / ${lab_name}`}),
        '<br/>',
        limitLengthSpan(lab_record_id, 30)
    ];
    const link = $('<a>', {
        href: Urls.view_classification(id),
        class: 'hover-link',
        html: content
    });

    let dom;
    if (data.matches) {
        const searchTerm = data.search;
        const match_doms = [link];
        const ekeys = EKeys.cachedKeys;
        for (const [key, value] of Object.entries(data.matches)) {
            match_doms.push($('<div>', {class: 'search-result', html:[
                $('<span>', { text: ekeys.key(key).label + ' : '}),
                    highlightTextAsDom(searchTerm, value)
                ]
            }));
        }
        dom = $('<div>', {html: match_doms});
    } else if (link) {
        dom = link;
    } else {
        dom = $('<div>', {html: content});
    }

    const indicatorClassName = `allele-origin-indicator allele-origin-horizontal allele-origin-${allele_origin_bucket}`;
    const fullDom = $('<div>', {style: 'display:flex; flex-direction:row', html:[
        alleleOriginDiv,
        dom
    ]});
    return fullDom;
};

VCTable.sample = (sample_name, type, row) => {
    let dom = $('<span/>');
    if (sample_name && row.sample_id) {
        const link = $('<a>', {
            href: Urls.view_sample(row.sample_id),
            class: 'hover-link',
            html: [
                $('<span>', {text: sample_name}),
            ]
        });
        dom = link;
    }
    return dom;
};
