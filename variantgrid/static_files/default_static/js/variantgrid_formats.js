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
// room for
const CLINVAR_SOMATIC_TIER_CHIPS = {
    'tier_1': {text: 'I', css: 'scs-tier_1'},
    'tier_1_or_2': {text: 'I/II', css: 'scs-tier_1_or_2'},
    'tier_2': {text: 'II', css: 'scs-tier_2'},
    'tier_3': {text: 'III', css: 'scs-tier_3'},
    'tier_4': {text: 'IV', css: 'scs-tier_4'},
};

// Each chip sits in a fixed slot (global.scss) so the four origins line up down the column - a
// reader scans one origin without their eye tracking sideways
const CLINVAR_CHIP_SLOT = 'cs-chip-clinvar';
const GERMLINE_CHIP_SLOT = 'cs-chip-germline';
const CLINVAR_SOMATIC_CHIP_SLOT = 'cs-chip-clinvar-somatic';
const SOMATIC_CHIP_SLOT = 'cs-chip-somatic';

function _classificationChip(cssClasses, innerHtml, title, gridColumn) {
    const href = gridColumn ? `javascript:showGridCell("${gridColumn}")` : 'javascript:void(0)';
    return `<a class='cs-chip ${cssClasses.join(' ')}' title='${escapeHtml(title)}' href='${href}'>${innerHtml}</a>`;
}

function _emptyChip(slotCssClass, title) {
    return _classificationChip([slotCssClass, 'cs-chip-empty'], '&mdash;', title, null);
}

function _internalClassificationChip(originLabel, originCssClass, slotCssClass, maxClassification,
                                     classifiedSummary, gridColumn, chipLookup, boxLookup) {
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
                               `${originLabel}: ${summaryLabels.join(' | ')}`, gridColumn);
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
    return _classificationChip([CLINVAR_CHIP_SLOT, chip.css], inner, title, "clinvar__clinical_significance");
}

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
    const gridColumn = oncChip ? "clinvar__highest_oncogenicity" : "clinvar__somatic_tier";
    const inner = `<span class='cs-chip-src'>CV</span>${escapeHtml(chip.text)}`
                + _clinvarStarsHtml(reviewStatus, ctx);
    return _classificationChip([CLINVAR_SOMATIC_CHIP_SLOT, chip.css], inner,
                               "ClinVar somatic: " + titleParts.join(", "), gridColumn);
}


// Renderer-only column (no row key of its own) - reads the hidden fields listed in
// CLASSIFICATIONS_COLUMN_ROW_FIELDS / _ANNOTATIONS (snpdb/grid_columns/custom_columns.py)
VariantGridFormat.classifications = (_value, type, rowData, ctx) => {
    // Always four chips in the same order, empty ones included - germline pair then somatic pair.
    // The .cs-chips grid gives each a fixed share of the cell, so they line up down the column
    // whatever the column is set to
    const chips = _clinvarChip(rowData, ctx)
        + _internalClassificationChip("Internally Classified (Germline)", "allele-origin-G", GERMLINE_CHIP_SLOT,
            rowData["max_internal_classification"], rowData["internally_classified"],
            "max_internal_classification", GERMLINE_CLASSIFICATION_CHIPS, GERMLINE_CLASSIFICATION_BOXES)
        + _clinvarSomaticChip(rowData, ctx)
        + _internalClassificationChip("Internally Classified (Somatic)", "allele-origin-S", SOMATIC_CHIP_SLOT,
            rowData["max_internal_somatic_classification"], rowData["internally_classified_somatic"],
            "max_internal_somatic_classification", SOMATIC_CLASSIFICATION_CHIPS, SOMATIC_CLASSIFICATION_BOXES);
    return `<span class='cs-chips'>${chips}</span>`;
};


// Impact drives the dot colour; the consequence text is what the cell reads as
const IMPACT_DOT_CSS = {
    'HIGH': 'impact-high',
    'MODERATE': 'impact-moderate',
    'LOW': 'impact-low',
    'MODIFIER': 'impact-modifier',
};

// Impact + Consequence in one cell - reads variantannotation__impact from COMPOSITE_COLUMN_ROW_FIELDS
// (snpdb/grid_columns/custom_columns.py). Sorts by consequence, the column it lives on.
VariantGridFormat.impactConsequence = (consequence, type, rowData) => {
    const impact = rowData["variantannotation__impact"];
    if (consequence == null && impact == null) {
        return '';
    }
    const dotCss = IMPACT_DOT_CSS[impact] || 'impact-unknown';
    const title = [impact, consequence].filter(v => v != null).join(' · ');
    const impactWord = impact == null ? '' : `<span class='ic-line2'>${escapeHtml(impact)}</span>`;
    return `<span class='impact-consequence' title='${escapeHtml(title)}'>`
         + `<i class='impact-dot ${dotCss}'></i>${escapeHtml(consequence == null ? '' : consequence)}`
         + `${impactWord}</span>`;
};

// The sample's call as a glyph, with the allele frequency and read depth beside it. The zygosity
// text comes through the server side formatter (Zygosity.CHOICES); the AF and depth ride along on
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

VariantGridFormat.sampleZygosity = (zygosity, type, rowData, ctx) => {
    const glyph = ZYGOSITY_GLYPHS[zygosity];
    if (!glyph) {
        // '.' - the sample has no call for this variant at all
        return '';
    }
    const prefix = (ctx && ctx.kwargs && ctx.kwargs.samplePrefix) || '';
    const named = [["AF", rowData[`${prefix}samples_allele_frequency`]],
                   ["depth", rowData[`${prefix}samples_read_depth`]]].filter(([, v]) => _hasSampleValue(v));
    const title = [zygosity].concat(named.map(([name, v]) => `${name} ${v}`)).join(' · ');
    let html = `<span class='sample-zygosity' title='${escapeHtml(title)}'>`
             + `<svg class='zyg-glyph ${ZYGOSITY_GLYPH_CSS[zygosity]}' viewBox='0 0 16 16'>${glyph}</svg>`;
    if (named.length) {
        html += `<span class='zyg-values'>${escapeHtml(named.map(([, v]) => v).join(' '))}</span>`;
    }
    return html + '</span>';
};


// SpliceAI: the max delta score with a dot at the standard 0.2 / 0.5 / 0.8 cutoffs, and the position
// offset of whichever of the four predictions the max came from. Reads the eight spliceai_pred_*
// ride-alongs; sorts on spliceai_max_ds, the indexed column it lives on.
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

VariantGridFormat.spliceai = (maxDs, type, rowData) => {
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
VariantGridFormat.maxentscan = (percentDiffRef, type, rowData) => {
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
VariantGridFormat.mastermind = (cdnaCount, type, rowData) => {
    const counts = [cdnaCount,
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


// Pathogenicity predictions as a segmented meter - one segment per tool that made a call, damaging
// first. Reads variantannotation__predictions_num_benign; sorts on the damaging count.
VariantGridFormat.predictions = (numPathogenic, type, rowData) => {
    const damaging = Number(numPathogenic) || 0;
    const benign = Number(rowData["variantannotation__predictions_num_benign"]) || 0;
    const total = damaging + benign;
    if (!total) {
        return '';
    }
    const segments = '<i class="pred-damaging"></i>'.repeat(damaging)
                   + '<i class="pred-benign"></i>'.repeat(benign);
    const title = `${damaging} damaging · ${benign} benign of ${total} prediction tools`;
    return `<span class='predictions' title='${escapeHtml(title)}'>${segments}`
         + `<span class='predictions-count'>${damaging}/${total}</span></span>`;
};


// This database's zygosity counts in one cell - hom · het, with ref and unknown in the title.
// Sorts on the het count, the column it lives on.
VariantGridFormat.dbZygosityCounts = (hetCount, type, rowData) => {
    const hom = rowData["global_variant_zygosity__hom_count"];
    if ((hetCount == null || hetCount === '') && (hom == null || hom === '')) {
        return '';
    }
    const counts = {hom: hom, het: hetCount,
                    ref: rowData["global_variant_zygosity__ref_count"],
                    unknown: rowData["global_variant_zygosity__unk_count"]};
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


// gnomAD popmax AF with the population it came from - reads variantannotation__gnomad_popmax.
// The AF itself is formatted server side (unit -> percent) so the CSV matches the grid.
VariantGridFormat.gnomadPopmax = (popmaxAf, type, rowData) => {
    if (popmaxAf == null || popmaxAf === '') {
        return '';
    }
    const population = rowData["variantannotation__gnomad_popmax"];
    let html = `<span class='gnomad-popmax-af'>${escapeHtml(popmaxAf)}</span>`;
    if (population) {
        html += ` <span class='gnomad-popmax-pop'>${escapeHtml(population)}</span>`;
    }
    return html;
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
    // 2. g.HGVS (no transcript - intergenic, or annotation not run yet for this variant)
    const g = _splitHgvs(hgvsG);
    if (g && g.change.length <= REPRESENTATIVE_MAX_HGVS_CHARS) {
        return {html: `<span class='rv-hgvs'>${escapeHtml(g.change)}</span> <span class='rv-sub'>${escapeHtml(chrom)}</span>`,
                title: hgvsG};
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

// Mandatory Variant column. Reads VARIANT_COLUMN_ROW_FIELDS (snpdb/grid_columns/custom_columns.py).
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


// gnomAD AF with the gnomAD Pass/Fail link beside it - reads variantannotation__gnomad_filtered from
// COMPOSITE_COLUMN_ROW_FIELDS (snpdb/grid_columns/custom_columns.py). The AF is formatted server side
// (unit -> percent) so the CSV matches the grid. Sorts by frequency, the column it lives on.
VariantGridFormat.gnomadAf = (af, type, rowData, ctx) => {
    const filteredLink = VariantGridFormat.gnomadFiltered(rowData["variantannotation__gnomad_filtered"],
                                                          type, rowData, ctx);
    if (af == null || af === '') {
        return filteredLink;
    }
    // At or above the import 'common' filter's AF the variant is one everybody carries - mute it so
    // the rare frequencies down the column are the ones that catch the eye
    const commonAf = ctx && ctx.extra && ctx.extra.commonGnomadAf;
    const common = commonAf != null && parseFloat(af) >= commonAf ? ' gnomad-af-common' : '';
    const afHtml = `<span class='gnomad-af${common}'>${escapeHtml(af)}</span>`;
    return filteredLink ? `${afHtml} ${filteredLink}` : afHtml;
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
