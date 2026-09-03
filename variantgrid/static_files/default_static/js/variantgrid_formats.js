/* Client renderers for the variant grids (the AbstractVariantGrid family).

   DataTables renderer signature is (data, type, row, ctx), where ctx is added by DataTableDefinition:
     ctx.extra    grid wide metadata from the definition JSON (DatatableConfig.get_extra)
     ctx.kwargs   this column's own settings, where it was declared with client_renderer_kwargs

   Shared page helpers (createGridLink, IGV, tags, load_variant_details) live in grid.js.
   A column names its renderer in RichColumn(client_renderer=...) - @see snpdb/grids.py */

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


// Multi-value columns arrive joined by the separator their VEP column declares (@see
// VEPColumnDef.separator) - the exports carry the raw value, the cell reads as a list
const VEP_SEPARATOR = "&";  // what a column falls back to, @see vep_field_formatters.VEP_SEPARATOR
const VALUE_JOIN = ", ";

function _separator(ctx) {
    return (ctx && ctx.kwargs && ctx.kwargs.separator) || VEP_SEPARATOR;
}

function splitValues(rawValue, separator) {
    return String(rawValue).split(separator || VEP_SEPARATOR).filter(v => v !== '');
}

function splitAndLink(rawValue, separator, buildLinkFunc) {
    if (!rawValue) {
        return '';
    }
    return splitValues(rawValue, separator).map(buildLinkFunc).join(VALUE_JOIN);
}

// Columns with nothing to link - the value is just easier to read split up
VariantGridFormat.separated = (value, type, rowData, ctx) => {
    if (value == null || value === '') {
        return '';
    }
    return escapeHtml(splitValues(value, _separator(ctx)).join(VALUE_JOIN));
};

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

// Chip text and .cs-*/.scs-* colour class (global.scss) per classification code
const GERMLINE_CLASSIFICATION_CHIPS = {
    '0': {text: 'O', css: 'cs-none'},
    '1': {text: 'B', css: 'cs-b'},
    '2': {text: 'LB', css: 'cs-lb'},
    '3': {text: 'VUS', css: 'cs-vus'},
    '4': {text: 'LP', css: 'cs-lp'},
    '5': {text: 'P', css: 'cs-p'},
};

const SOMATIC_CLASSIFICATION_CHIPS = {
    'tier_1': {text: 'Tier I', css: 'scs-tier_1'},
    'tier_1_or_2': {text: 'Tier I/II', css: 'scs-tier_1_or_2'},
    'tier_2': {text: 'Tier II', css: 'scs-tier_2'},
    'tier_3': {text: 'Tier III', css: 'scs-tier_3'},
    'tier_4': {text: 'Tier IV', css: 'scs-tier_4'},
};

const CLINVAR_PATHOGENICITY_CHIPS = {  // ClinVar.highest_pathogenicity (ClinVarPathogenicity)
    0: {text: 'CV', css: 'cs-none'},
    1: {text: 'B', css: 'cs-b'},
    2: {text: 'LB', css: 'cs-lb'},
    3: {text: 'VUS', css: 'cs-vus'},
    4: {text: 'LP', css: 'cs-lp'},
    5: {text: 'P', css: 'cs-p'},
};

const CLINVAR_ONCOGENICITY_CHIPS = {  // ClinVar.highest_oncogenicity (ClinVarOncogenicity)
    1: {text: 'B', css: 'cs-b'},
    2: {text: 'LB', css: 'cs-lb'},
    3: {text: 'VUS', css: 'cs-vus'},
    4: {text: 'LO', css: 'cs-lp'},
    5: {text: 'O', css: 'cs-p'},
};

// ClinVar.somatic_tier (SomaticClinicalSignificance) - the AMP tier, in the short form the chip has
// room for. It's a choice field, so unlike the internal tiers above the row carries the label
const CLINVAR_SOMATIC_TIER_CHIPS = {
    'Tier I': {text: 'I', css: 'scs-tier_1'},
    'Tier I/II': {text: 'I/II', css: 'scs-tier_1_or_2'},
    'Tier II': {text: 'II', css: 'scs-tier_2'},
    'Tier III': {text: 'III', css: 'scs-tier_3'},
    'Tier IV': {text: 'IV', css: 'scs-tier_4'},
};

// Each chip sits in a fixed slot (global.scss) so the four origins line up down the column - a
// reader scans one origin without their eye tracking sideways
const CLINVAR_CHIP_SLOT = 'cs-chip-clinvar';
const GERMLINE_CHIP_SLOT = 'cs-chip-germline';
const CLINVAR_SOMATIC_CHIP_SLOT = 'cs-chip-clinvar-somatic';
const SOMATIC_CHIP_SLOT = 'cs-chip-somatic';

// Chips are spans rather than links, so a click on one falls through to the row click handler and
// expands the row - where the full ClinVar and classification detail is (@see variantGridRowDetail)
function _classificationChip(cssClasses, innerHtml, title) {
    return `<span class='cs-chip ${cssClasses.join(' ')}' title='${escapeHtml(title)}'>${innerHtml}</span>`;
}

function _emptyChip(slotCssClass, title) {
    return _classificationChip([slotCssClass, 'cs-chip-empty'], '&mdash;', title);
}

function _internalClassificationChip(originLabel, originCssClass, slotCssClass, maxClassification,
                                     classifiedSummary, chipLookup, boxLookup) {
    // null = not classified; undefined = a row cached before these fields existed - treat the same
    if (maxClassification == null) {
        return _emptyChip(slotCssClass, originLabel + ": not classified");
    }
    const chip = chipLookup[maxClassification] || {text: maxClassification, css: 'cs-none'};
    const records = (classifiedSummary || '').split('|');
    const summaryLabels = records.map(cs => (boxLookup[cs] || {label: cs}).label);
    let inner = escapeHtml(chip.text);
    if (records.length > 1) {
        inner += ` <span class='cs-chip-count'>&times;${records.length}</span>`;
    }
    return _classificationChip([slotCssClass, chip.css, originCssClass], inner,
                               `${originLabel}: ${summaryLabels.join(' | ')}`);
}

function _clinvarChip(rowData, ctx) {
    const highestPath = rowData["clinvar__highest_pathogenicity"];
    if (highestPath == null) {
        return _emptyChip(CLINVAR_CHIP_SLOT, "ClinVar: not classified");
    }
    const chip = CLINVAR_PATHOGENICITY_CHIPS[highestPath] || {text: String(highestPath), css: 'cs-none'};
    const inner = `<span class='cs-chip-src'>CV</span>${escapeHtml(chip.text)}`
                + _clinvarStarsHtml(rowData["clinvar__review_status"], ctx);
    const title = "ClinVar: " + (rowData["clinvar__clinical_significance"] || chip.text);
    return _classificationChip([CLINVAR_CHIP_SLOT, chip.css], inner, title);
}

// clinvarStars is keyed by the review status label, which is what the choice field's column carries
// @see variant_grid_client_extra
function _clinvarStarsHtml(reviewStatus, ctx) {
    const starsLookup = (ctx && ctx.extra && ctx.extra.clinvarStars) || {};
    const stars = starsLookup[reviewStatus];
    if (stars === undefined) {
        return '';
    }
    let starHtml = '';
    for (let i = 0; i < 4; ++i) {
        starHtml += `<span class='${i < stars ? '' : 'off'}'>&#9733;</span>`;
    }
    return ` <span class='cs-chip-stars'>${starHtml}</span>`;
}

// ClinVar's somatic view of the variant. Oncogenicity is the call most records carry, so it leads;
// a record with only an AMP tier shows the tier instead. Both go in the tooltip either way.
function _clinvarSomaticChip(rowData, ctx) {
    const oncogenicity = rowData["clinvar__highest_oncogenicity"];
    const tier = rowData["clinvar__somatic_tier"];
    const oncChip = CLINVAR_ONCOGENICITY_CHIPS[oncogenicity];
    const tierChip = CLINVAR_SOMATIC_TIER_CHIPS[tier];
    if (!oncChip && !tierChip) {
        return _emptyChip(CLINVAR_SOMATIC_CHIP_SLOT, "ClinVar somatic: no call");
    }
    const titleParts = [];
    if (oncChip) {
        titleParts.push("oncogenicity " + (rowData["clinvar__oncogenic_classification"] || oncChip.text));
    }
    if (tierChip) {
        titleParts.push("AMP tier " + tierChip.text);
    }
    const chip = oncChip || tierChip;
    const reviewStatus = oncChip ? rowData["clinvar__oncogenic_review_status"]
                                 : rowData["clinvar__somatic_review_status"];
    const inner = `<span class='cs-chip-src'>CV</span>${escapeHtml(chip.text)}`
                + _clinvarStarsHtml(reviewStatus, ctx);
    return _classificationChip([CLINVAR_SOMATIC_CHIP_SLOT, chip.css], inner,
                               "ClinVar somatic: " + titleParts.join(", "));
}


// Renderer-only column (no row key of its own) - reads its member columns, which ride along hidden
// @see CompositeColumnMember
VariantGridFormat.classifications = (_value, type, rowData, ctx) => {
    // Always four chips in the same order, empty ones included - germline pair then somatic pair.
    // The .cs-chips grid gives each a fixed share of the cell, so they line up down the column
    // whatever the column is set to
    const chips = _clinvarChip(rowData, ctx)
        + _internalClassificationChip("Internally Classified (Germline)", "allele-origin-G", GERMLINE_CHIP_SLOT,
            rowData["max_internal_classification"], rowData["internally_classified"],
            GERMLINE_CLASSIFICATION_CHIPS, GERMLINE_CLASSIFICATION_BOXES)
        + _clinvarSomaticChip(rowData, ctx)
        + _internalClassificationChip("Internally Classified (Somatic)", "allele-origin-S", SOMATIC_CHIP_SLOT,
            rowData["max_internal_somatic_classification"], rowData["internally_classified_somatic"],
            SOMATIC_CLASSIFICATION_CHIPS, SOMATIC_CLASSIFICATION_BOXES);
    return `<span class='cs-chips'>${chips}</span>`;
};


// The generic composite cell: the first member is what the cell reads as, every other non-blank
// member goes on hover as "label: value". A group whose members are all blank draws nothing, which
// is what lets the eye skip down a sparse column.
// Members (path, label, the client renderer the member carries standalone, and its separator where
// it holds multiple values) come from the members render kwarg
// @see _composite_column_kwargs in snpdb/grid_columns/custom_columns.py
function _compositeMemberHtml(member, value, type, rowData, ctx) {
    if (member.renderer) {
        // The member's own separator, not the composite's
        const memberCtx = {...ctx, kwargs: {...(ctx && ctx.kwargs), separator: member.separator}};
        return eval(member.renderer)(value, type, rowData, memberCtx);
    }
    return escapeHtml(_memberText(member, value));
}

function _memberText(member, value) {
    return member.separator ? splitValues(value, member.separator).join(VALUE_JOIN) : String(value);
}

// The separator declared by the composite member reading `path`
function _memberSeparator(ctx, path) {
    const members = (ctx && ctx.kwargs && ctx.kwargs.members) || [];
    const member = members.find(m => m.path === path);
    return member && member.separator;
}

function _isBlank(value) {
    return value == null || value === '';
}

VariantGridFormat.composite = (_value, type, rowData, ctx) => {
    const members = (ctx && ctx.kwargs && ctx.kwargs.members) || [];
    if (!members.length) {
        return '';
    }
    const detail = [];
    for (let i = 1; i < members.length; ++i) {
        const value = rowData[members[i].path];
        if (!_isBlank(value)) {
            detail.push(`${members[i].label}: ${_memberText(members[i], value)}`);
        }
    }
    const headline = rowData[members[0].path];
    if (_isBlank(headline) && !detail.length) {
        return '';
    }
    // A group with detail but no headline still has to be hoverable, so it keeps a muted placeholder
    const cell = _isBlank(headline)
        ? `<span class='composite-empty-headline'>&middot;</span>`
        : _compositeMemberHtml(members[0], headline, type, rowData, ctx);
    const headlineText = _isBlank(headline) ? '' : _memberText(members[0], headline);
    const title = [`${members[0].label}: ${headlineText}`].concat(detail).join(' \u00b7 ');
    return `<span class='composite-cell' title='${escapeHtml(title)}'>${cell}</span>`;
};


// Impact drives the dot colour; the consequence text is what the cell reads as
const IMPACT_DOT_CSS = {
    'HIGH': 'impact-high',
    'MODERATE': 'impact-moderate',
    'LOW': 'impact-low',
    'MODIFIER': 'impact-modifier',
};

// Impact + Consequence in one cell. Sorts by consequence, its headline member.
VariantGridFormat.impactConsequence = (_value, type, rowData, ctx) => {
    const rawConsequence = rowData["variantannotation__consequence"];
    const impact = rowData["variantannotation__impact"];
    if (rawConsequence == null && impact == null) {
        return '';
    }
    // VEP gives every consequence the variant hits, joined - "missense_variant, splice_region_variant"
    const consequence = rawConsequence == null
        ? '' : splitValues(rawConsequence, _memberSeparator(ctx, "variantannotation__consequence")).join(VALUE_JOIN);
    const dotCss = IMPACT_DOT_CSS[impact] || 'impact-unknown';
    const title = [impact, consequence].filter(v => v !== '' && v != null).join(' · ');
    const impactWord = impact == null ? '' : `<span class='ic-line2'>${escapeHtml(impact)}</span>`;
    return `<span class='impact-consequence' title='${escapeHtml(title)}'>`
         + `<i class='impact-dot ${dotCss}'></i>${escapeHtml(consequence)}`
         + `${impactWord}</span>`;
};

// The sample's whole call in one cell. The zygosity text comes through the server side formatter
// (Zygosity.CHOICES); everything else - AF, AD, DP, GQ, PL and the sample filters - rides along on
// this sample's own hidden columns, named by the samplePrefix render kwarg.
// @see get_grid_genotype_columns_and_overrides in analysis/grids.py
const ZYGOSITY_GLYPHS = {
    'HET': '<circle cx="8" cy="8" r="6.5" fill="none" stroke="currentColor" stroke-width="1.6"/>'
         + '<path d="M8 1.5 A6.5 6.5 0 0 1 8 14.5 Z" fill="currentColor"/>',
    'HOM_ALT': '<circle cx="8" cy="8" r="6.5" fill="currentColor"/>',
    'REF': '<circle cx="8" cy="8" r="6.5" fill="none" stroke="currentColor" stroke-width="1.6"/>',
    '?': '<circle cx="8" cy="8" r="6.5" fill="none" stroke="currentColor" stroke-width="1.6" stroke-dasharray="2 2"/>',
};
const ZYGOSITY_GLYPH_CSS = {
    'HET': 'zyg-het',
    'HOM_ALT': 'zyg-hom-alt',
    'REF': 'zyg-ref',
    '?': 'zyg-unknown',
};

function _hasSampleValue(value) {
    return value != null && value !== '' && value !== '.';
}

// GQ and PL are only worth a reader's attention when they're bad, so the cell marks them rather
// than spelling them out - the numbers are in the cell's title. GQ is a minimum (higher is better),
// PL the phred likelihood of the called genotype, which reads the other way.
// @see VARIANT_GRID_GENOTYPE_QUALITY_THRESHOLDS
const GENOTYPE_QUALITIES = [
    {key: 'gq', field: 'samples_genotype_quality', label: 'GQ', higherIsBetter: true},
    {key: 'pl', field: 'samples_phred_likelihood', label: 'PL', higherIsBetter: false},
];

function _genotypeQualityCss(value, thresholds, higherIsBetter) {
    if (!thresholds) {
        return 'zyg-qual-good';
    }
    if (higherIsBetter ? value <= thresholds.bad : value >= thresholds.bad) {
        return 'zyg-qual-bad';
    }
    if (higherIsBetter ? value < thresholds.marginal : value > thresholds.marginal) {
        return 'zyg-qual-marginal';
    }
    return 'zyg-qual-good';
}

VariantGridFormat.sampleZygosity = (zygosity, type, rowData, ctx) => {
    const glyph = ZYGOSITY_GLYPHS[zygosity];
    if (!glyph) {
        // '.' - the sample has no call for this variant at all
        return '';
    }
    const prefix = (ctx && ctx.kwargs && ctx.kwargs.samplePrefix) || '';
    const sampleValue = (field) => rowData[`${prefix}${field}`];
    const detail = [zygosity];

    // Allele frequency, then the alt reads over the total behind it
    const af = sampleValue('samples_allele_frequency');
    const ad = sampleValue('samples_allele_depth');
    const dp = sampleValue('samples_read_depth');
    let values = '';
    if (_hasSampleValue(af)) {
        values += `<span class='zyg-af'>${escapeHtml(String(af))}</span>`;
        detail.push(`AF ${af}`);
    }
    const depths = [];
    if (_hasSampleValue(ad)) {
        depths.push(ad);
        detail.push(`AD ${ad}`);
    }
    if (_hasSampleValue(dp)) {
        depths.push(dp);
        detail.push(`DP ${dp}`);
    }
    if (depths.length) {
        values += `<span class='zyg-depth'>${escapeHtml(depths.join('/'))}</span>`;
    }

    const thresholds = (ctx && ctx.extra && ctx.extra.genotypeQuality) || {};
    let quality = '';
    for (const q of GENOTYPE_QUALITIES) {
        const value = sampleValue(q.field);
        if (!_hasSampleValue(value)) {
            continue;
        }
        const css = _genotypeQualityCss(parseFloat(value), thresholds[q.key], q.higherIsBetter);
        quality += `<i class='zyg-qual ${css}' title='${q.label} ${escapeHtml(String(value))}'></i>`;
        detail.push(`${q.label} ${value}`);
    }

    // Passing every filter is the common case and says nothing - only a call that failed one shows
    const filters = sampleValue('samples_filters');
    let filtersHtml = '';
    if (_hasSampleValue(filters) && String(filters).toUpperCase() !== 'PASS') {
        filtersHtml = `<i class='zyg-filtered fa-solid fa-triangle-exclamation' `
                    + `title='${escapeHtml(`Filtered: ${filters}`)}'></i>`;
        detail.push(`FT ${filters}`);
    }

    return `<span class='sample-zygosity' title='${escapeHtml(detail.join(' · '))}'>`
         + `<svg class='zyg-glyph ${ZYGOSITY_GLYPH_CSS[zygosity]}' viewBox='0 0 16 16'>${glyph}</svg>`
         + (values ? `<span class='zyg-values'>${values}</span>` : '')
         + quality + filtersHtml + '</span>';
};


// SpliceAI: the max delta score with a dot at the standard 0.2 / 0.5 / 0.8 cutoffs, and the position
// offset of whichever of the four predictions the max came from. Reads the eight spliceai_pred_*
// members riding along; sorts on spliceai_max_ds, its headline member.
const SPLICEAI_PREDICTIONS = [
    ["ag", "acceptor gain"],
    ["al", "acceptor loss"],
    ["dg", "donor gain"],
    ["dl", "donor loss"],
];

function _spliceaiThresholdCss(maxDs) {
    if (maxDs >= 0.8) {
        return 'spliceai-high';
    }
    if (maxDs >= 0.5) {
        return 'spliceai-moderate';
    }
    if (maxDs >= 0.2) {
        return 'spliceai-low';
    }
    return 'spliceai-none';
}

VariantGridFormat.spliceai = (_value, type, rowData) => {
    const maxDs = rowData["variantannotation__spliceai_max_ds"];
    if (maxDs == null || maxDs === '') {
        return '';
    }
    const detail = [];
    for (const [suffix, label] of SPLICEAI_PREDICTIONS) {
        const ds = rowData[`variantannotation__spliceai_pred_ds_${suffix}`];
        if (ds == null) {
            continue;
        }
        const dp = rowData[`variantannotation__spliceai_pred_dp_${suffix}`];
        detail.push(`${label} ${ds}${dp == null ? '' : ` @ ${dp}`}`);
    }
    const title = detail.length ? detail.join(', ') : `SpliceAI max Δ ${maxDs}`;
    return `<span class='spliceai' title='${escapeHtml(title)}'>`
         + `<i class='spliceai-dot ${_spliceaiThresholdCss(maxDs)}'></i>${escapeHtml(maxDs)}</span>`;
};


// MaxEntScan: the drop from the reference splice site score, with the raw ref -> alt scores in the
// title. Negative means the variant weakens the site, which is the direction that matters.
VariantGridFormat.maxentscan = (_value, type, rowData) => {
    const percentDiffRef = rowData["variantannotation__maxentscan_percent_diff_ref"];
    if (percentDiffRef == null || percentDiffRef === '') {
        return '';
    }
    const ref = rowData["variantannotation__maxentscan_ref"];
    const alt = rowData["variantannotation__maxentscan_alt"];
    const diff = rowData["variantannotation__maxentscan_diff"];
    const parts = [];
    if (ref != null && alt != null) {
        parts.push(`ref ${ref} → alt ${alt}`);
    }
    if (diff != null) {
        parts.push(`diff ${diff}`);
    }
    const title = parts.length ? parts.join(', ') : `${percentDiffRef}% of the reference score`;
    const weakened = Number(percentDiffRef) < 0 ? ' maxentscan-weakened' : '';
    return `<span class='maxentscan${weakened}' title='${escapeHtml(title)}'>${escapeHtml(percentDiffRef)}%</span>`;
};


// Mastermind's three granularities in one cell - cDNA · cDNA+protein · AA change - linked to the
// article list through the MMID3 ride-along.
VariantGridFormat.mastermind = (_value, type, rowData) => {
    const counts = [rowData["variantannotation__mastermind_count_1_cdna"],
                    rowData["variantannotation__mastermind_count_2_cdna_prot"],
                    rowData["variantannotation__mastermind_count_3_aa_change"]];
    if (counts.every(c => c == null || c === '')) {
        return '';
    }
    const labels = ["cDNA", "cDNA + protein", "AA change"];
    const title = counts.map((c, i) => `${labels[i]}: ${c == null ? 0 : c}`).join(', ');
    const text = counts.map(c => c == null ? 0 : c).join(' · ');
    const body = `<span class='mastermind-counts' title='${escapeHtml(title)}'>${escapeHtml(text)}</span>`;
    const mmid3 = rowData["variantannotation__mastermind_mmid3"];
    if (!mmid3) {
        return body;
    }
    const url = `https://mastermind.genomenon.com/detail?mutation=${encodeURIComponent(String(mmid3).split("&")[0])}`;
    return `<a class='mastermind-link' href='${escapeHtml(url)}' target='_blank'>${body}</a>`;
};


// ALoFT's whole call in one cell - the prediction, a dot for how damaging it is and the probability
// behind it. A low confidence call (ALoFT's p >= 0.05) recedes rather than reading as a firm answer.
// Reads the probabilities, the confidence flag and the chosen transcript from its members;
// sorts on the prediction, its headline member.
// The choices are expanded server side, so the value here is 'Tolerant' / 'Recessive' / 'Dominant'
const ALOFT_PREDICTIONS = {
    'dominant': {css: 'aloft-dominant', field: 'variantannotation__aloft_prob_dominant'},
    'recessive': {css: 'aloft-recessive', field: 'variantannotation__aloft_prob_recessive'},
    'tolerant': {css: 'aloft-tolerant', field: 'variantannotation__aloft_prob_tolerant'},
};
const ALOFT_PROBABILITIES = [
    ["Tol", "variantannotation__aloft_prob_tolerant"],
    ["Rec", "variantannotation__aloft_prob_recessive"],
    ["Dom", "variantannotation__aloft_prob_dominant"],
];

VariantGridFormat.aloft = (_value, type, rowData) => {
    const prediction = rowData["variantannotation__aloft_pred"];
    if (prediction == null || prediction === '') {
        return '';
    }
    const called = ALOFT_PREDICTIONS[String(prediction).toLowerCase()];
    // null is no call on confidence either way, which is not the same as a call ALoFT doubts
    const highConfidence = rowData["variantannotation__aloft_high_confidence"];
    const detail = [prediction];
    if (highConfidence != null) {
        detail.push(highConfidence ? 'high confidence' : 'low confidence');
    }

    const probabilities = ALOFT_PROBABILITIES
        .filter(([, field]) => rowData[field] != null && rowData[field] !== '')
        .map(([label, field]) => `${label} ${rowData[field]}`);
    if (probabilities.length) {
        detail.push(probabilities.join(' / '));
    }
    const transcript = rowData["variantannotation__aloft_ensembl_transcript"];
    if (transcript) {
        detail.push(transcript);
    }

    // The probability of the call itself - the other two, and the full precision, are on hover.
    // ALoFT calls the highest of the three, so this one never rounds away to nothing
    const calledProbability = called && rowData[called.field] != null && rowData[called.field] !== ''
        ? Number(rowData[called.field]) : NaN;
    const probabilityHtml = isFinite(calledProbability)
        ? `<span class='aloft-prob'>${calledProbability.toFixed(2)}</span>` : '';
    const confidenceCss = highConfidence === false ? ' aloft-low-confidence' : '';
    return `<span class='aloft${confidenceCss}' title='${escapeHtml(detail.join(' · '))}'>`
         + `<i class='aloft-dot ${called ? called.css : 'aloft-unknown'}'></i>`
         + `${escapeHtml(String(prediction))}${probabilityHtml}</span>`;
};


// Pathogenicity predictions as a segmented meter - one segment per tool that made a call, damaging
// first, with the call each tool made on hover. Sorts on the damaging count, its headline member.
const PREDICTION_TOOL_CALLS = [
    ["SIFT", "variantannotation__sift"],
    ["Polyphen2 HVAR", "variantannotation__polyphen2_hvar_pred_most_damaging"],
    ["Mutation Taster", "variantannotation__mutation_taster_pred_most_damaging"],
    ["Mutation Assessor", "variantannotation__mutation_assessor_pred_most_damaging"],
    ["FATHMM", "variantannotation__fathmm_pred_most_damaging"],
    ["MetaLR rankscore", "variantannotation__metalr_rankscore"],
];

VariantGridFormat.predictions = (_value, type, rowData) => {
    const damaging = Number(rowData["variantannotation__predictions_num_pathogenic"]) || 0;
    const benign = Number(rowData["variantannotation__predictions_num_benign"]) || 0;
    const total = damaging + benign;
    if (!total) {
        return '';
    }
    const segments = '<i class="pred-damaging"></i>'.repeat(damaging)
                   + '<i class="pred-benign"></i>'.repeat(benign);
    const calls = PREDICTION_TOOL_CALLS
        .filter(([, field]) => rowData[field] != null && rowData[field] !== '')
        .map(([label, field]) => `${label} ${rowData[field]}`);
    const title = [`${damaging} damaging \u00b7 ${benign} benign of ${total} prediction tools`]
        .concat(calls).join(', ');
    return `<span class='predictions' title='${escapeHtml(title)}'>${segments}`
         + `<span class='predictions-count'>${damaging}/${total}</span></span>`;
};


// The population frequencies other than gnomAD - 1000 Genomes, UK10K, TOPMed - as the highest of them
// with which one it came from, each on hover. The members this version annotates (and their labels)
// come from the members render kwarg; their values arrive formatted server side, in the deployment's
// unit or percent, so they compare as numbers and print as they are. Sorts on the headline member.
VariantGridFormat.popFreqOther = (_value, type, rowData, ctx) => {
    const members = (ctx && ctx.kwargs && ctx.kwargs.members) || [];
    const detail = [];
    let highest = null;
    for (const member of members) {
        const value = rowData[member.path];
        if (_isBlank(value)) {
            continue;
        }
        detail.push(`${member.label} ${value}`);
        if (highest === null || parseFloat(value) > parseFloat(highest.value)) {
            highest = {source: member.label.replace(/ AF$/, ''), value};
        }
    }
    if (highest === null) {
        return '';
    }
    return `<span class='pop-freq-other' title='${escapeHtml(detail.join(' · '))}'>`
         + `${escapeHtml(highest.value)}<span class='pop-freq-source'>${escapeHtml(highest.source)}</span></span>`;
};


// UniProt accession, linked to the entry. The uniprot cell reuses it for its headline
const UNIPROT_ACCESSION = "variantannotation__transcript_version__gene_version__hgnc__uniprot__accession";
const UNIPROT_FUNCTION = "variantannotation__transcript_version__gene_version__hgnc__uniprot__function";

VariantGridFormat.uniprotLink = (accession) => {
    if (_isBlank(accession)) {
        return '';
    }
    const url = `https://www.uniprot.org/uniprotkb/${encodeURIComponent(String(accession))}/entry`;
    return `<a class='uniprot-link' href='${escapeHtml(url)}' target='_blank'>${escapeHtml(accession)}</a>`;
};

// The gene's UniProt entry in one cell - the linked accession with the function summary running on
// after it, and the function, pathway, tissue specificity and Reactome pathways on hover
VariantGridFormat.uniprot = (_value, type, rowData, ctx) => {
    const members = (ctx && ctx.kwargs && ctx.kwargs.members) || [];
    const accession = rowData[UNIPROT_ACCESSION];
    const detail = [];
    for (const member of members) {
        const value = rowData[member.path];
        if (member.path !== UNIPROT_ACCESSION && !_isBlank(value)) {
            detail.push(`${member.label}: ${value}`);
        }
    }
    if (_isBlank(accession) && !detail.length) {
        return '';
    }
    const func = rowData[UNIPROT_FUNCTION];
    const funcHtml = _isBlank(func) ? '' : `<span class='uniprot-function'>${escapeHtml(func)}</span>`;
    return `<span class='uniprot' title='${escapeHtml(detail.join(' · '))}'>`
         + `${VariantGridFormat.uniprotLink(accession)}${funcHtml}</span>`;
};


// Conservation as a row of dots, one per score this version annotates, in member order - filled where
// the score reaches its conserved threshold, muted below it, hollow where the position has no score.
// The cell carries a detail table that the hover handler below floats under it: each score between
// its (faded) minimum and maximum. Ranges and thresholds come from ctx.extra.conservation
// (@see VariantAnnotation.CONSERVATION_SCORES); sorts on PhyloP 100 way, its headline member.
function _conservationRange(value) {
    return value == null ? '' : String(Math.round(value * 100) / 100);
}

// Scores come out of a float4 column, so 0.314 arrives as 0.314000010490417
function _conservationScore(value) {
    const score = Number(value);
    return isFinite(score) ? String(parseFloat(score.toPrecision(3))) : String(value);
}

VariantGridFormat.conservation = (_value, type, rowData, ctx) => {
    const members = (ctx && ctx.kwargs && ctx.kwargs.members) || [];
    const scales = (ctx && ctx.extra && ctx.extra.conservation) || {};
    let dots = '';
    let rows = '';
    let anyScore = false;
    for (const member of members) {
        const value = rowData[member.path];
        const scale = scales[member.path] || {};
        let css = 'cons-none';
        let valueHtml = `<span class='cons-blank'>–</span>`;
        if (!_isBlank(value)) {
            anyScore = true;
            const conserved = scale.conserved != null && parseFloat(value) >= scale.conserved;
            css = conserved ? 'cons-conserved' : 'cons-not-conserved';
            valueHtml = escapeHtml(_conservationScore(value));
        }
        dots += `<i class='cons-dot ${css}'></i>`;
        rows += `<tr><th>${escapeHtml(member.label)}</th>`
              + `<td class='cons-range'>${escapeHtml(_conservationRange(scale.min))}</td>`
              + `<td class='cons-value'><i class='cons-dot ${css}'></i>${valueHtml}</td>`
              + `<td class='cons-range'>${escapeHtml(_conservationRange(scale.max))}</td></tr>`;
    }
    if (!anyScore) {
        return '';
    }
    return `<span class='conservation'>${dots}<table class='conservation-detail'>${rows}</table></span>`;
};

// The conservation cell's hover - its detail table, floated below the cell. A panel a reader opened
// on purpose (a header sort menu) stays put
$(document).on('mouseenter', '.conservation', function() {
    if (FloatingPanel.panel && !FloatingPanel.panel.hasClass('conservation-detail')) {
        return;
    }
    FloatingPanel.show($(this).children('.conservation-detail').clone(), this);
}).on('mouseleave', '.conservation', function() {
    if (FloatingPanel.panel && FloatingPanel.panel.hasClass('conservation-detail')) {
        FloatingPanel.hide();
    }
});


// Zygosity counts in one cell - hom · het, with the rest in the title. Sorts on the het count, the
// column it lives on. The row key prefix comes from the countPrefix render kwarg: this database's
// global counts by default, a cohort node's own counts where it names its own.
// @see CohortNode._get_node_extra_columns
const DB_ZYGOSITY_COUNT_PREFIX = "global_variant_zygosity__";

VariantGridFormat.dbZygosityCounts = (_value, type, rowData, ctx) => {
    const prefix = (ctx && ctx.kwargs && ctx.kwargs.countPrefix) || DB_ZYGOSITY_COUNT_PREFIX;
    const hetCount = rowData[`${prefix}het_count`];
    const hom = rowData[`${prefix}hom_count`];
    if ((hetCount == null || hetCount === '') && (hom == null || hom === '')) {
        return '';
    }
    // A cohort node counts calls it made, so it has no unknown count - it drops out of the title
    const counts = {hom: hom, het: hetCount,
                    ref: rowData[`${prefix}ref_count`],
                    unknown: rowData[`${prefix}unk_count`]};
    const title = Object.entries(counts)
        .filter(([, v]) => v != null && v !== '')
        .map(([k, v]) => `${v} ${k}`).join(', ');
    const cell = ["hom", "het"].map(k => {
        const v = counts[k] == null || counts[k] === '' ? 0 : counts[k];
        const zero = Number(v) ? '' : ' db-count-zero';
        return `<span class='db-count${zero}'>${escapeHtml(v)}</span>`;
    }).join(' · ');
    return `<span class='db-zygosity-counts' title='${escapeHtml(title)}'>${cell}</span>`;
};


// The record's VCF FILTER. Nearly every row passed, and a column of 'PASS' is a column of nothing -
// so a pass fades almost out and the codes a record failed read as ordinary text.
// @see CohortMixin._get_node_extra_columns
VariantGridFormat.vcfFilters = (filters) => {
    if (filters == null || filters === '' || filters === '.') {
        return '';
    }
    const text = String(filters);
    if (text.toUpperCase() === 'PASS') {
        return `<span class='vcf-filter-pass' title='Passed every VCF filter'>PASS</span>`;
    }
    return `<span class='vcf-filter-failed' title='${escapeHtml(`Filtered: ${text}`)}'>${escapeHtml(text)}</span>`;
};


// hgvs_c / hgvs_p / hgvs_g are "ACCESSION:change" or "ACCESSION(SYMBOL):change". The change has to
// start with an HGVS kind so gene-level fusion nomenclature ("BCR::ABL1") isn't split as one
const HGVS_REGEX = /^([^:(]+)(?:\(([^)]+)\))?:([cgmnopr]\..+)$/;
const HGVS_NOT_CALCULATED = "HGVS not calculated due to length";  // VariantAnnotation.SV_HGVS_TOO_LONG_MESSAGE
const REPRESENTATIVE_MAX_ALLELE_BASES = 10;   // longer ref/alt collapse to "[52bp]"
const REPRESENTATIVE_MAX_HGVS_CHARS = 40;     // longer g.HGVS (inserted sequence spelled out) fall through to the coordinate

function _splitHgvs(hgvs) {
    const m = hgvs ? String(hgvs).match(HGVS_REGEX) : null;
    return m ? {accession: m[1], symbol: m[2] || null, change: m[3]} : null;
}

function _formatAllele(seq) {
    seq = seq == null ? '' : String(seq);
    return seq.length > REPRESENTATIVE_MAX_ALLELE_BASES ? `[${seq.length}bp]` : seq;
}

// The contig the way a reader says it - 'chr8' rather than the bare '8' the contig is named
function _contigLabel(chrom) {
    if (chrom == null || chrom === '') {
        return '';
    }
    const name = String(chrom);
    return name.toLowerCase().startsWith('chr') ? name : `chr${name}`;
}

function _formatBases(bases) {
    if (bases >= 1e6) {
        return (bases / 1e6).toFixed(2) + " Mb";
    }
    if (bases >= 1e3) {
        return (bases / 1e3).toFixed(1) + " kb";
    }
    return bases + " bp";
}

// The cascade from the mockup's "Representative variant" card. Returns {html, title}; title is the
// plain string for the link tooltip. Every key may be undefined on a row cached before these fields
// existed - each step checks and falls through, ending at the VariantGrid id.
function _representativeVariantLabel(variantId, rowData) {
    const chrom = rowData["locus__contig__name"];
    const position = rowData["locus__position"];
    const ref = rowData["locus__ref__seq"];
    const alt = rowData["alt__seq"];
    const svlen = rowData["svlen"];
    const hgvsC = rowData["variantannotation__hgvs_c"];
    const hgvsP = rowData["variantannotation__hgvs_p"];
    const hgvsG = rowData["variantannotation__hgvs_g"];
    const symbol = rowData["variantannotation__symbol"];

    // 1. c.HGVS on the representative transcript - gene symbol leads, transcript is in the expanded row
    if (hgvsC && hgvsC !== HGVS_NOT_CALCULATED) {
        const c = _splitHgvs(hgvsC);
        if (!c) {
            // Gene-level (fusion) nomenclature has no accession prefix - show it whole
            return {html: `<span class='rv-hgvs'>${escapeHtml(hgvsC)}</span>`, title: hgvsC};
        }
        const gene = c.symbol || symbol;
        let html = gene ? `<span class='rv-gene'>${escapeHtml(gene)}</span> ` : '';
        html += `<span class='rv-hgvs'>${escapeHtml(c.change)}</span>`;
        const p = _splitHgvs(hgvsP);
        if (p) {
            html += ` <span class='rv-hgvs-p'>${escapeHtml(p.change)}</span>`;
        }
        // Two line rows move the accession and the protein change below the name (global.scss picks
        // which of the two shows) - the identifiers a reader wants next, without a second column
        html += `<span class='rv-line2'>${escapeHtml(c.accession)}`
              + (p ? ` <span class='rv-hgvs-p'>${escapeHtml(p.change)}</span>` : '') + '</span>';
        return {html: html, title: hgvsC + (hgvsP ? " " + hgvsP : "")};
    }
    // 2. g.HGVS (no transcript - intergenic, or annotation not run yet for this variant). The contig
    // leads the way the gene symbol does above, and the accession it stands in for moves to line 2
    const g = _splitHgvs(hgvsG);
    if (g && g.change.length <= REPRESENTATIVE_MAX_HGVS_CHARS) {
        const contig = _contigLabel(chrom);
        const html = (contig ? `<span class='rv-gene'>${escapeHtml(contig)}</span> ` : '')
                   + `<span class='rv-hgvs'>${escapeHtml(g.change)}</span>`
                   + `<span class='rv-line2'>${escapeHtml(g.accession)}</span>`;
        return {html: html, title: hgvsG};
    }
    // 3/4. Coordinate - symbolic as a span with the type, otherwise ref>alt with long alleles collapsed
    if (chrom != null && position != null && alt != null) {
        if (String(alt).startsWith("<")) {
            const size = Math.abs(svlen || 0);
            const end = position + size;
            const svType = String(alt).slice(1, -1);
            const html = `<span class='rv-hgvs'>${escapeHtml(chrom)}:${position}-${end}</span> <b>${escapeHtml(svType)}</b>`
                       + (size ? ` <span class='rv-sub'>${_formatBases(size)}</span>` : '');
            return {html: html, title: `${chrom}:${position}-${end} ${alt}`};
        }
        const title = `${chrom}:${position} ${ref}>${alt}`;
        return {html: `<span class='rv-hgvs'>${escapeHtml(chrom)}:${position} ${escapeHtml(_formatAllele(ref))}&gt;${escapeHtml(_formatAllele(alt))}</span>`,
                title: title};
    }
    // 5. Last resort while annotation is still running - the id the old details box linked to
    return {html: `<span class='rv-sub'>v ${variantId}</span>`, title: `VariantGrid variant ${variantId}`};
}

// Mandatory Variant column. Reads the members riding along hidden - @see CompositeColumnMember.
// Markup contract: .variant_id-container[variant_id] > input.variant-select (analysis only)
//                  + a.variant-link (details; grid.js swaps its href to the full URL on right click).
// Clicking the row itself expands it - @see variantGridRowDetail in grid.js
VariantGridFormat.representativeVariant = (variantId, type, rowData, ctx) => {
    const parts = [];
    if (_isNodeVisible(ctx)) {
        parts.push(`<input type='checkbox' class='variant-select' variant_id='${variantId}'>`);
    }
    const label = _representativeVariantLabel(variantId, rowData);
    const detailsUrl = `javascript:load_variant_details(${variantId});`;
    parts.push(`<a class='variant-link rv-label' title='${escapeHtml(label.title)}' href='${detailsUrl}' orig_href='${detailsUrl}'>${label.html}</a>`);
    return `<span class='variant_id-container' variant_id='${variantId}'>${parts.join('')}</span>`;
};


// This is driven entirely off variantTags (not passed through SQL->the grid)
// This is so we can add/remove tags without wrecking cache
VariantGridFormat.tags = (tagsCellValue, type, rowData) => {
    const variantId = rowData['id'];
    let tagHtml = "";
    const aWin = getAnalysisWindow();
    const readOnly = aWin.variantTagsReadOnly || !inAnalysis();

    if (!readOnly) {
        tagHtml += "<a class='show-tag-autocomplete' variant_id='" + variantId + "' href='javascript:showTagAutocomplete(" + variantId + ")'><span class='add-variant-tag' title='Tag variant..'></span></a>";
    }

    const tagList = (aWin.variantTags || {})[variantId];
    if (tagList) {
        const sortedTags = sortVariantTags(aWin, tagList);
        for (let i=0 ; i<sortedTags.length ; ++i) {
            const tag = sortedTags[i];
            tagHtml += getVariantTagHtml(variantId, tag, readOnly);
        }
    }
    // Wrapped so the set can lift out of the clipped cell as one thing on hover - @see .grid-tags
    return tagHtml ? `<span class='grid-tags'>${tagHtml}</span>` : "";
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
    return tagGlobalHtml ? `<span class='grid-tags'>${tagGlobalHtml}</span>` : "";
};


VariantGridFormat.clinvarLink = (clinvar_variation_id) => {
    let clinvar_string = '';
    if (clinvar_variation_id) {
        clinvar_string = "<a title='View ClinVar entry in new window' target='_blank' href='http://www.ncbi.nlm.nih.gov/clinvar/variation/" + clinvar_variation_id + "'>" + clinvar_variation_id + "</a>";
    }
    return clinvar_string;
};


VariantGridFormat.cosmicLink = (cosmic_ids, type, rowData, ctx) => {
    const COSMIC_PREFIX = "COSV";
    const COSMIC_LEGACY_PREFIX = "COSM";

    let cosmic_string = '';
    if (cosmic_ids) {
        const cosmic_ids_list = splitValues(cosmic_ids, _separator(ctx));
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

        cosmic_string = cosmic_links.join(VALUE_JOIN);
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
                // The gene list node's own badge glyph - the child this makes is a GeneListNode
                geneLinkString = "<a class='grid-link' title='Gene list node for " + geneSymbol + "' href='" + filterGeneLink + "'>"
                               + "<svg class='gene-list-node-icon' width='14' height='14' viewBox='0 0 18 20'><use href='#node-icon-gene-list'></use></svg></a>";
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
    if (gnomadFilteredCellValue != null) {
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


// The whole gnomAD read on the variant in one cell: the overall AF, the popmax AF with the
// population it came from, and the Pass/Fail link on the cell's right edge. The per-population
// frequencies, allele counts, homozygotes and FAFs are on hover.
// The AFs are formatted server side (unit -> percent) so the CSV matches the grid. Sorts on
// gnomad_af, its headline member.
const GNOMAD_DETAIL = [
    ["AC", "variantannotation__gnomad_ac"],
    ["AN", "variantannotation__gnomad_an"],
    ["Hom alt", "variantannotation__gnomad_hom_alt"],
    ["Hemi", "variantannotation__gnomad_hemi_count"],
    ["PopMax AC", "variantannotation__gnomad_popmax_ac"],
    ["PopMax AN", "variantannotation__gnomad_popmax_an"],
    ["PopMax hom alt", "variantannotation__gnomad_popmax_hom_alt"],
    ["AFR", "variantannotation__gnomad_afr_af"],
    ["AMR", "variantannotation__gnomad_amr_af"],
    ["ASJ", "variantannotation__gnomad_asj_af"],
    ["EAS", "variantannotation__gnomad_eas_af"],
    ["FIN", "variantannotation__gnomad_fin_af"],
    ["MID", "variantannotation__gnomad_mid_af"],
    ["NFE", "variantannotation__gnomad_nfe_af"],
    ["OTH", "variantannotation__gnomad_oth_af"],
    ["SAS", "variantannotation__gnomad_sas_af"],
    ["XY AF", "variantannotation__gnomad_xy_af"],
    ["XY AC", "variantannotation__gnomad_xy_ac"],
    ["XY AN", "variantannotation__gnomad_xy_an"],
    ["non-PAR", "variantannotation__gnomad_non_par"],
    ["FAF95", "variantannotation__gnomad_faf95"],
    ["FAF99", "variantannotation__gnomad_faf99"],
    ["FAF95 max", "variantannotation__gnomad_fafmax_faf95_max"],
    ["FAF99 max", "variantannotation__gnomad_fafmax_faf99_max"],
    ["gnomAD2 AF", "variantannotation__gnomad2_liftover_af"],
];

VariantGridFormat.gnomad = (_value, type, rowData, ctx) => {
    const af = rowData["variantannotation__gnomad_af"];
    const popmaxAf = rowData["variantannotation__gnomad_popmax_af"];
    const population = rowData["variantannotation__gnomad_popmax"];
    const filteredLink = VariantGridFormat.gnomadFiltered(rowData["variantannotation__gnomad_filtered"],
                                                          type, rowData, ctx);
    const detail = [];
    if (af != null && af !== '') {
        detail.push(`AF ${af}`);
    }
    if (popmaxAf != null && popmaxAf !== '') {
        detail.push(`PopMax ${popmaxAf}${population ? ` (${population})` : ''}`);
    }
    for (const [label, field] of GNOMAD_DETAIL) {
        const value = rowData[field];
        if (value != null && value !== '') {
            detail.push(`${label} ${value}`);
        }
    }
    if (!detail.length) {
        return filteredLink ? `<span class='gnomad-af-cell'>${filteredLink}</span>` : '';
    }

    // At or above the import 'common' filter's AF the variant is one everybody carries - mute it so
    // the rare frequencies down the column are the ones that catch the eye
    const commonAf = ctx && ctx.extra && ctx.extra.commonGnomadAf;
    const common = commonAf != null && parseFloat(af) >= commonAf ? ' gnomad-af-common' : '';
    let html = af == null || af === '' ? '' : `<span class='gnomad-af${common}'>${escapeHtml(af)}</span>`;
    if (popmaxAf != null && popmaxAf !== '') {
        html += `<span class='gnomad-popmax-af'>${escapeHtml(popmaxAf)}</span>`;
        if (population) {
            const codes = (ctx && ctx.extra && ctx.extra.gnomadPopulationCodes) || {};
            html += `<span class='gnomad-popmax-pop'>${escapeHtml(codes[population] || population)}</span>`;
        }
    }
    // The Pass/Fail goes to the cell's right edge, so they line up down the column whatever width
    // the frequencies beside them take
    return `<span class='gnomad-af-cell' title='${escapeHtml(detail.join(' \u00b7 '))}'>${html}${filteredLink}</span>`;
};


VariantGridFormat.clinGenAlleleId = (cellValue) => {
    // warning: doesn't use settings.CLINGEN_ALLELE_REGISTRY_DOMAIN as static JS
    if (!cellValue) {
        return "";
    }
    const url = `http://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_caid?caid=${cellValue}`;
    return `<a href="${url}" target="_blank">${cellValue}</a>`;
};


VariantGridFormat.dbsnp = (dbsnp_rs_ids, type, rowData, ctx) => {
    function buildDBSNPLink(dbsnp_id) {
        return "<a title='View dbSNP in new window' target='_blank' href='https://www.ncbi.nlm.nih.gov/snp/" + dbsnp_id + "'>" + dbsnp_id + "</a>";
    }
    return splitAndLink(dbsnp_rs_ids, _separator(ctx), buildDBSNPLink);
};


VariantGridFormat.pubMed = (pubmed, type, rowData, ctx) => {
    function buildPubMedLink(pubmed_id) {
        return "<a title='View PubMed article in new window' target='_blank' href='https://pubmed.ncbi.nlm.nih.gov/" + pubmed_id + "'>" + pubmed_id + "</a>";
    }
    return splitAndLink(pubmed, _separator(ctx), buildPubMedLink);
};


VariantGridFormat.ontologyTerms = (ontology_terms) => {
    function buildOntologyTermLink(ontology_term) {
        const termSlug = ontology_term.split(" ")[0].replace(":", "_");
        const url = Urls.ontology_term(termSlug);
        return "<a title='View Ontology Term in new window' target='_blank' href='" + url + "'>" + ontology_term + "</a>";
    }
    return splitAndLink(ontology_terms, " | ", buildOntologyTermLink);
};


VariantGridFormat.masterMind = (value, type, rowData, ctx) => {
    function buildMasterMindLink(mmid3) {
        return "<a title='View MasterMind in new window' target='_blank' href='https://mastermind.genomenon.com/detail?mutation=" + mmid3 + "'>" + mmid3 + "</a>";
    }
    return splitAndLink(value, _separator(ctx), buildMasterMindLink);
};


