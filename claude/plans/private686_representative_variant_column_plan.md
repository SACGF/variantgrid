# Representative Variant column for the variant grids (variantgrid_private#686)

Written by Claude Fable 5 (claude-fable-5), 2026-08-31.
Implemented and revised by Claude Opus 5 (claude-opus-5), 2026-08-31 - see "As built" at the end.

Mockup: `claude/mockups/private686_variant_grid_scannable.html` (open it in a browser - the "Proposed default
grid" table and the "Anatomy" cards are the visual spec for the cell and the chips). Where the mockup shows a
popover, this plan uses the project's standard expandable row instead (see Step 3c).

## Goal

The mandatory `Variant` column currently renders five 16px boxes (details, ClinVar, germline, somatic, IGV)
and the grid then needs four more columns (chrom / position / ref / alt) before the reader knows *which*
variant a row is. Replace that with one cell that names the variant the way a curator would say it -
`TP53 c.743G>A p.(Arg248Gln)` - picked by a preference cascade and sorted in genomic order. Expanding the row
(click, as on the other DataTables grids) shows the less-scanned identifiers: transcript, protein, genomic
HGVS, coordinate in every build, ClinGen allele, VariantGrid id, the other transcripts, and the Details /
new tab / IGV actions. The classification boxes move into their own `Classifications` column of chips.

Scope of this plan (the mockup's "first slice" plus the Classifications column that gives the displaced
boxes a home):

1. Hidden row fields for the two renderers
2. `Classifications` column (new optional `VariantGridColumn` + chip renderer)
3. `Variant` column renderer (cascade) replacing `detailsLink`, plus the expandable row
4. Genomic sort for the `Variant` column
5. Data migration: new column, system collections, cache bump
6. Tests, cleanup, verification

The other glyph columns in the mockup (Sample zygosity glyphs, Impact·Consequence, gnomAD muting,
Predictions meter, SpliceAI dot, two-line rows, header sort menus) are separate follow-up issues - each is
an independent renderer on columns that already exist.

## Where it applies

Every grid built on `AbstractVariantGrid` (`snpdb/grids.py:492`) gets all of this, because the `id` override,
the hidden row fields, the sort branch and the expand renderer live on the base class and the collection:

| Grid | Page | Mounted at |
|---|---|---|
| `analysis.grids.VariantGrid` | analysis node grid (and `ExportVariantGrid` for CSV/VCF) | `analysis/views/views_grid.py` |
| `variantopedia.grids.AllVariantsGrid` | Variantopedia "All variants" | `variantopedia/urls.py:91` |
| `variantopedia.grids.NearbyVariantsGrid` | variant details "nearby" tabs | `variantopedia/urls.py:85-88` |
| `variantopedia.grids.TaggedVariantGrid` | "Variants with tags" (variant-centric, `tags_global` shows every tag) | `variantopedia/urls.py:95` |
| `genes.grids.GeneSymbolVariantsGrid` | gene symbol page variants tab | `genes/urls.py:77` |

`GeneSymbolVariantsGrid._get_non_gene_fields` drops fields containing `__transcript_version__` / `__gene__`;
none of the new hidden fields match, so they survive. The tag-centric `VariantTagsColumns` table
(`variantopedia/grids.py`, a `DatatableConfig`) is a different table with its own renderer and is untouched.

## How the grid works today (read these first)

The variant grids are the `library.jqgrid` engine served to a DataTables client through a thin adapter.
Nothing in this plan changes the engine - the feature lives in colmodel overrides, one sort branch, one
small view and JS.

| What | Where |
|---|---|
| Column catalogue (`VariantGridColumn`), collections (`CustomColumnsCollection`, `CustomColumn`), mandatory `variant` column | `snpdb/models/models_columns.py` |
| Collection → grid fields/overrides, plus the hidden fields `detailsLink` needs (`ID_FORMATTER_REQUIRED_FIELDS`) | `snpdb/grid_columns/custom_columns.py:15-78` |
| Per-allele classification subqueries (`internally_classified`, `max_internal_classification`, somatic twins, `tags_global`) | `snpdb/grid_columns/custom_columns.py:81` `get_variantgrid_extra_annotate` |
| `AbstractVariantGrid` - the `id` override (`formatter: 'detailsLink'`, width 110), standard formatters, `get_queryset`, `get_datatable_extra` | `snpdb/grids.py:492-618` |
| Analysis node grid `VariantGrid` (sort limits, packed genotype sort, `ExportVariantGrid`) | `analysis/grids.py:76-440` |
| Engine sort (`_sort_items`: `F(sidx)` + secondary `sortname` + `-pk`), `get_colmodels`, `field_to_colmodel` | `library/jqgrid/jqgrid.py:275-303, 434-457, 547` |
| Adapter: formatter name → JS renderer map, colmodel → DataTables column (`hidden` → `visible:false`), definition JSON, `order[0][column]` → `sidx` | `library/django_utils/jqgrid_datatable_adapter.py:33-51, 54-72, 132-176, 208-232` |
| Client renderers (`VariantGridFormat.detailsLink` etc.) | `variantgrid/static_files/default_static/js/variantgrid_formats.js` |
| Page helpers: `load_variant_details`, `createGridLink`, `createIgvUrl`/`create_igv_link`, `showGridCell`, `gridCompleteExtra` (right-click → full URL on `a.variant-link`, reads `.variant_id-container[variant_id]`) | `variantgrid/static_files/default_static/js/grid.js:88-110, 227-240, 243-297, 385-395` |
| Analysis page: `input.variant-select[variant_id]` checkbox wiring | `analysis/templates/analysis/node_data/node_data_grid.html:12-54` |
| Renderer wrapper: `(data, type, row, ctx)` with `ctx.extra` = `grid.get_datatable_extra()`, `ctx.kwargs` = colmodel `formatter_kwargs` | `variantgrid/static_files/default_static/js/datatable_definition.js:309-330` |
| **Row expansion**: `defn.expandClientRenderer` → first column gets `toggle-link`; `setupClientExpend` toggles `row.child(expandFn(row.data()))` on row click, prefetches after a 500ms hover | `datatable_definition.js:333-336, 498-570`; arrow/cursor CSS `global.scss:456, 575-600, 1128-1140, 1197` |
| Existing ajax row-detail pattern: `DatatableConfig._row_expand_ajax(view_name, id_field, expected_height)` → `TableFormat.expandAjax` → `loadAjaxBlock(dom, Urls[view](id))` | `snpdb/views/datatable_view.py:319`, `datatable_definition.js:850-890`, e.g. `snpdb/views/views.py:907` `manual_variant_entry_collection_detail` |
| Node grid rows and definition are `cache_page`d for a week; `ccc_version_id` is in the URL | `analysis/views/views_grid.py:71-130` |
| Grid CSS (`table.variantgrid-datatable`, `.variant_id-container`, `.grid-link*`, `.c-pill`/`.cs-*` colours) | `variantgrid/static_files/default_static/css/global.scss:883-1020, 3204-3290` |
| Data-migration pattern for new columns + collection edits | `snpdb/migrations/0123_new_vg_columns_and_custom_columns.py` |
| Deployment check that every `VariantGridColumn.variant_column` is a real queryset path | `variantgrid/deployment_validation/column_check.py` (`SPECIAL_COLUMNS`) |

Facts established while researching, which the design relies on:

- `VariantGridColumn` rows of interest: `variant` → `id`; `chrom` → `locus__contig__name`; `position` →
  `locus__position`; `ref` → `locus__ref__seq`; `alt` → `alt__seq`; `svlen` → `svlen`; `hgvs_c` →
  `variantannotation__hgvs_c`; `hgvs_p` → `variantannotation__hgvs_p`; `hgvs_g` → `variantannotation__hgvs_g`;
  `symbol` → `variantannotation__symbol`; `clinvar_review_status` → `clinvar__review_status`.
- System collections (`user IS NULL`): `Default columns`, `Minimalism`, `All columns`. All start
  `variant, chrom, position, ref, alt, ...`; `Default columns` also has `svlen, hgvs_g` right after.
- `hgvs_c` in the DB comes in two shapes: `NM_199135.4:c.466A>C` and `NM_016329.4(SFMBT1):c.-130-31769_-130-30020del`.
  `hgvs_p` is `NP_954586.4:p.Ile156Leu`; `hgvs_g` is `NC_000009.11:g.70918333A>C`. Structural variants can
  have `hgvs_c = NULL` or the sentinel `VariantAnnotation.SV_HGVS_TOO_LONG_MESSAGE` (`"HGVS not calculated due to length"`).
  Gene-level (fusion) variants carry VICC nomenclature in `hgvs_c` with no `accession:` prefix.
- Symbolic variants: `alt__seq` is `<DEL>`/`<DUP>`/`<INS>`/..., `svlen` is the signed length (`-1750` for a
  1750bp deletion), position + |svlen| is the end.
- Contig ids are **not** in genome-build order for GRCh38 (`MT` is contig id 25, shared with GRCh37), so
  "sort by contig id" is wrong for the last chromosome. `GenomeBuild.standard_contigs` is ordered by
  `genomebuildcontig__order`. A `CASE` over those ids annotated onto a node `.values()` queryset and used in
  `order_by` was run against a real GRCh38 node and works.
- DataTables here is 1.13.1. Its `fnGetData` is `s&&e ? s(a,e,t,n) : a` - the renderer **is** called when the
  row has no key for the column (that is how the `tags` column works: `queryset_field=False`, renderer only).
  DataTables 2.x changed this, so a future upgrade needs `defaultContent` on renderer-only columns.
- Hidden colmodels are exported in CSV (`grid_export_csv` walks every colmodel), so the labels of hidden
  fields matter: `locus__ref__seq` and `alt__seq` would both come out as `seq`.
- `setupClientExpend` binds the toggle to the whole `tr` and only skips rows without `odd`/`even` classes.
  Variant grid rows are full of links and a checkbox, so it needs a guard for clicks that start on a control.
- Every `AbstractVariantGrid` subclass sets `self.genome_build`; all but the analysis `VariantGrid` also set
  `self.annotation_version` (the analysis one has `node.analysis.annotation_version`).
- `Variant.all_build_variants` (`models_variant.py:769`) gives the same mutation in every build via the
  allele; `Variant.allele.clingen_allele` gives the ClinGen id; `VariantTranscriptAnnotation` holds every
  transcript's annotation per `VariantAnnotationVersion`.

## Design decisions

- **One renderer, client side.** Everything the cascade needs is already in the row for other columns
  (or is one more hidden field on a join that is already there). No per-row fetch to draw the grid.
- **Cascade** (from the mockup's "Representative variant" card): c.HGVS (gene symbol leads, `p.` change muted
  after it) → g.HGVS (change only, chrom muted) → coordinate with long alleles collapsed to `[52bp]`,
  symbolic as `1:155000000-155120000 DEL · 120 kb` → `v 48812331`.
- **Click the name = details; click the row = expand.** The name is the details link exactly as the old
  details box was (`load_variant_details`; right-click gives the full URL via the existing `gridCompleteExtra`
  wiring). Clicking anywhere else on the row expands it, using `DataTableDefinition`'s existing
  `expandClientRenderer` mechanism, so the variant grids behave like every other expandable grid on the
  site (arrow in the first column, one row open at a time). The expanded row is a small server-rendered
  fragment for the annotation version the grid is showing, plus a client-built action bar
  (Details / Open in new tab / IGV) so IGV keeps following `ANALYSIS_SETTINGS.show_igv_links`.
- **Sort is genomic** whatever the cell shows: `(genome build contig order, contig id, position, ref, alt, pk)`.
  Server side only - `sidx == "id"` is the trigger, so the analysis default sort column can be set to
  `variant` and it just works.
- **Classifications is an ordinary optional column**, renderer-only like `tags`/`tags_global`. It reads the
  same hidden fields `detailsLink` reads today plus `clinvar__review_status` for the stars. It goes into every
  system collection right after `variant`.
- **Collections:** `chrom`, `position`, `ref`, `alt`, `hgvs_g` leave `Default columns` and `Minimalism`
  (they stay in `All columns` and remain in the chooser for everyone). `hgvs_c` / `hgvs_p` / `gene_symbol`
  stay - the full transcript-qualified HGVS is what people copy out, and the gene column carries the
  filter-child link. User-cloned collections are left alone.
- **Cache:** the renderers tolerate rows cached before this change (missing keys → cascade falls through to
  coordinate / id; chips treat `undefined` like "not classified"), and `CACHE_VERSION` is bumped so the
  week-long node grid definition cache stops referencing `VariantGridFormat.detailsLink`.

## Step 1 - Hidden row fields (server)

`snpdb/grid_columns/custom_columns.py`. Replace `ID_FORMATTER_REQUIRED_FIELDS` / `ID_FORMATTER_REQUIRED_ANNOTATIONS`
with module-level constants (they are referenced from tests and the renderers' comments):

```python
from snpdb.models import CustomColumn, CustomColumnsCollection, VariantGridColumn

# Row fields the client renderers read whichever columns the collection shows - they ride along hidden.
# @see VariantGridFormat.representativeVariant / VariantGridFormat.classifications in variantgrid_formats.js
VARIANT_COLUMN_ROW_FIELDS = [
    "locus__contig__name",
    "locus__position",
    "locus__ref__seq",
    "alt__seq",
    "svlen",
    "variantannotation__symbol",
    "variantannotation__hgvs_c",
    "variantannotation__hgvs_p",
    "variantannotation__hgvs_g",
]
CLASSIFICATIONS_COLUMN_ROW_FIELDS = [
    "clinvar__highest_pathogenicity",
    "clinvar__clinical_significance",
    "clinvar__review_status",
]
# These come from get_variantgrid_extra_annotate rather than the model, so they are colmodel-only
# (not in the values() queryset)
CLASSIFICATIONS_COLUMN_ROW_ANNOTATIONS = [
    "internally_classified",
    "max_internal_classification",
    "internally_classified_somatic",
    "max_internal_somatic_classification",
]
```

and at the end of `get_custom_column_fields_override_and_sample_position`:

```python
    hidden_fields = VARIANT_COLUMN_ROW_FIELDS + CLASSIFICATIONS_COLUMN_ROW_FIELDS
    # A hidden model field still gets a CSV header - use the catalogue label rather than the model
    # verbose_name (locus__ref__seq and alt__seq are both 'seq')
    hidden_labels = dict(VariantGridColumn.objects.filter(variant_column__in=hidden_fields)
                         .values_list("variant_column", "label"))
    for field in hidden_fields + CLASSIFICATIONS_COLUMN_ROW_ANNOTATIONS:
        if field not in fields:
            fields.append(field)
            if field in CLASSIFICATIONS_COLUMN_ROW_ANNOTATIONS:
                ov = get_overrides([field], [{}], model_field=False, queryset_field=False)[field]
            else:
                ov = {}
                if label := hidden_labels.get(field):
                    ov["label"] = label
            ov["hidden"] = True
            override[field] = ov
```

Notes:
- `variantannotation__*` fields add no join - every variant grid already joins `variantannotation` (the
  base querysets filter on annotation version). `locus__ref__seq` / `alt__seq` add two PK joins to
  `snpdb_sequence`; `svlen` is on `Variant`.
- These bypass the VEP columns-version gate on purpose (same as today's required fields); `hgvs_c`/`hgvs_p`
  have no `min_columns_version`, `hgvs_g` is computed at insert time.
- Fields like `variantannotation__hgvs_c` that *are* in the collection are simply not appended (the
  `if field not in fields` guard) and keep their visible colmodel - the renderer reads the same row key.

## Step 2 - Classifications column

### 2a. Catalogue row + standard override

The `VariantGridColumn` (inserted by the migration in Step 5):

```python
{'grid_column_name': 'classifications',
 'variant_column': 'classifications',
 'annotation_level': 'D',   # DATABASE_LEVEL, like tags
 'width': 170,
 'label': 'Classifications',
 'description': 'ClinVar, internal germline and internal somatic classifications. Hover a chip for the '
                'per-record summary, click it to scroll to the underlying column',
 'model_field': False,
 'queryset_field': False}
```

`AbstractVariantGrid._get_standard_overrides` in `snpdb/grids.py` - alongside `tags_global`:

```python
            'classifications': {
                'model_field': False, 'queryset_field': False,
                'name': 'classifications', 'index': 'classifications',
                'classes': 'no-word-wrap', 'formatter': 'classificationsFormatter', 'sortable': False,
            },
```

`sortable: False` matters: the base `_sort_items` would otherwise try `F("classifications")` and log a
`FieldError` on every header click.

Add `"classifications"` to `SPECIAL_COLUMNS` in `variantgrid/deployment_validation/column_check.py`,
next to `"tags", "tags_global"` - otherwise the deployment check reports the new column as an invalid
queryset path.

### 2b. Stars mapping via the definition, once per grid

`ClinVarReviewStatus.STARS` is a Python dict; ship it in `get_datatable_extra` rather than copying it into
JS. `snpdb/grids.py` `AbstractVariantGrid`:

```python
from annotation.models.models_enums import ClinVarReviewStatus

    def get_datatable_extra(self) -> dict:
        # gnomAD links are per genome build, and the client renderers have no other way to know it
        return {"genomeBuild": self.genome_build.name,
                "clinvarStars": dict(ClinVarReviewStatus.STARS)}
```

`VariantGrid.get_datatable_extra` (analysis) calls `super()` so it inherits this.

### 2c. Adapter mapping

`library/django_utils/jqgrid_datatable_adapter.py` `JQGRID_FORMATTER_TO_CLIENT_RENDERER` - keep it sorted:

```python
    "classificationsFormatter": "VariantGridFormat.classifications",
```

### 2d. Renderer

`variantgrid_formats.js`. `GERMLINE_CLASSIFICATION_BOXES` / `SOMATIC_CLASSIFICATION_BOXES` already map
codes → `{display, label}`; reuse them for titles and add chip text + colour class:

```js
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

function _classificationChip(cssClasses, innerHtml, title, gridColumn) {
    const href = gridColumn ? `javascript:showGridCell("${gridColumn}")` : 'javascript:void(0)';
    return `<a class='cs-chip ${cssClasses.join(' ')}' title='${escapeHtml(title)}' href='${href}'>${innerHtml}</a>`;
}

function _emptyChip(title) {
    return _classificationChip(['cs-chip-empty'], '&mdash;', title, null);
}

function _internalClassificationChip(originLabel, originCssClass, maxClassification, classifiedSummary,
                                     gridColumn, chipLookup, boxLookup) {
    // null = not classified; undefined = a row cached before these fields existed - treat the same
    if (maxClassification == null) {
        return _emptyChip(originLabel + ": not classified");
    }
    const chip = chipLookup[maxClassification] || {text: maxClassification, css: 'cs-none'};
    const records = (classifiedSummary || '').split('|');
    const summaryLabels = records.map(cs => (boxLookup[cs] || {label: cs}).label);
    let inner = escapeHtml(chip.text);
    if (records.length > 1) {
        inner += ` <span class='cs-chip-count'>&times;${records.length}</span>`;
    }
    return _classificationChip([chip.css, originCssClass], inner,
                               `${originLabel}: ${summaryLabels.join(' | ')}`, gridColumn);
}

function _clinvarChip(rowData, ctx) {
    const highestPath = rowData["clinvar__highest_pathogenicity"];
    if (highestPath == null) {
        return _emptyChip("ClinVar: not classified");
    }
    const chip = CLINVAR_PATHOGENICITY_CHIPS[highestPath] || {text: String(highestPath), css: 'cs-none'};
    const starsLookup = (ctx && ctx.extra && ctx.extra.clinvarStars) || {};
    const stars = starsLookup[rowData["clinvar__review_status"]];
    let inner = `<span class='cs-chip-src'>CV</span>${escapeHtml(chip.text)}`;
    if (stars !== undefined) {
        let starHtml = '';
        for (let i = 0; i < 4; ++i) {
            starHtml += `<span class='${i < stars ? '' : 'off'}'>&#9733;</span>`;
        }
        inner += ` <span class='cs-chip-stars'>${starHtml}</span>`;
    }
    const title = "ClinVar: " + (rowData["clinvar__clinical_significance"] || chip.text);
    return _classificationChip([chip.css], inner, title, "clinvar__clinical_significance");
}

// Renderer-only column (no row key of its own) - reads the hidden fields listed in
// CLASSIFICATIONS_COLUMN_ROW_FIELDS / _ANNOTATIONS (snpdb/grid_columns/custom_columns.py)
VariantGridFormat.classifications = (_value, type, rowData, ctx) => {
    return _clinvarChip(rowData, ctx)
        + _internalClassificationChip("Internally Classified (Germline)", "allele-origin-G",
            rowData["max_internal_classification"], rowData["internally_classified"],
            "max_internal_classification", GERMLINE_CLASSIFICATION_CHIPS, GERMLINE_CLASSIFICATION_BOXES)
        + _internalClassificationChip("Internally Classified (Somatic)", "allele-origin-S",
            rowData["max_internal_somatic_classification"], rowData["internally_classified_somatic"],
            "max_internal_somatic_classification", SOMATIC_CLASSIFICATION_CHIPS, SOMATIC_CLASSIFICATION_BOXES);
};
```

For the CV chip at `highest_pathogenicity == 0` (drug response / risk factor / oncogenic-only records) the chip
reads `CV` on grey with the raw `clinical_significance` in the tooltip - the same information the `0` box
gives today. ClinVar somatic / oncogenic chips are a follow-up.

### 2e. Styles

`global.scss` (then hand-apply to `global.css` per CLAUDE.md). Hook the chips into the shared colour block:

```scss
// Clinical significance colours - shared by classification pills, analysis node chips and grid chips
.c-pill, .node-chip, .cs-chip {
	&.scs-none, &.cs-none { ...unchanged... }
```

and next to `.variant_id-container` inside `table.variantgrid-datatable`:

```scss
	.cs-chip {
		display: inline-block;
		height: 16px;
		line-height: 16px;
		padding: 0 5px;
		border-radius: 4px;
		font-size: 10px;
		font-weight: 700;
		color: #333;
		margin-right: 3px;
		vertical-align: middle;
		text-decoration: none;

		&.cs-chip-empty {
			background: transparent;
			border: 1px dashed #ccc;
			color: #bbb;
			font-weight: 500;
		}

		// Same germline / somatic colours as .allele-origin-box, as a left rail
		&.allele-origin-G { box-shadow: inset 3px 0 0 #88cc88; padding-left: 7px; }
		&.allele-origin-S { box-shadow: inset 3px 0 0 #cc99cc; padding-left: 7px; }

		.cs-chip-src { font-size: 9px; opacity: .6; margin-right: 2px; }
		.cs-chip-stars { color: #d9a400; font-size: 9px; letter-spacing: -1px; .off { color: #ccc; } }
		.cs-chip-count { font-size: 9px; background: rgba(0, 0, 0, .08); border-radius: 99px; padding: 0 4px; }
	}
```

The existing `.grid-link-icon.allele-origin-G/S` rules (`global.scss:906-912`) set `background-color` on
`.grid-link-icon` only, so they do not collide with `.cs-chip.allele-origin-G`.

## Step 3 - Variant column renderer and expandable row

### 3a. Colmodel

`snpdb/grids.py` `_get_standard_overrides`, replacing the `id` entry and its comment:

```python
            # The representative variant cell: expand arrow, select checkbox, cascade label (details link)
            'id': {'editable': False, 'width': 280, 'fixed': True, 'formatter': 'representativeVariant',
                   'sorttype': 'int'},
```

280 = the old 110 plus the chrom/position/ref/alt widths (20 + 70 + 40 + 40) that leave the default
collections, so the leading cluster stays the same width (the 16px expand arrow comes out of that).
`sorttype: 'int'` stays - it drives the filter builder's input type for the (pre-existing) "filter by
Variant id" rule.

Adapter map: replace the `detailsLink` line with

```python
    "representativeVariant": "VariantGridFormat.representativeVariant",
```

### 3b. Renderer

`variantgrid_formats.js` - `VariantGridFormat.detailsLink`, `internalClassificationBox` and their comment go;
this takes their place. Keep the `.variant_id-container[variant_id]` wrapper, the `variant-select` checkbox
and the `a.variant-link` class - `grid.js:227-240, 385-395` and `node_data_grid.html:12-54` depend on them.

```js
// hgvs_c / hgvs_p / hgvs_g are "ACCESSION:change" or "ACCESSION(SYMBOL):change"
const HGVS_REGEX = /^([^:(]+)(?:\(([^)]+)\))?:(.+)$/;
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
    if (bases >= 1e6) return (bases / 1e6).toFixed(2) + " Mb";
    if (bases >= 1e3) return (bases / 1e3).toFixed(1) + " kb";
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
    let nodeVisible = _isNodeVisible(ctx);
    const kwargs = ctx && ctx.kwargs;
    if (kwargs) {
        nodeVisible = kwargs.node_visible;
    }
    const parts = [];
    if (nodeVisible) {
        parts.push(`<input type='checkbox' class='variant-select' variant_id='${variantId}'>`);
    }
    const label = _representativeVariantLabel(variantId, rowData);
    const detailsUrl = `javascript:load_variant_details(${variantId});`;
    parts.push(`<a class='variant-link rv-label' title='${escapeHtml(label.title)}' href='${detailsUrl}' orig_href='${detailsUrl}'>${label.html}</a>`);
    return `<span class='variant_id-container' variant_id='${variantId}'>${parts.join('')}</span>`;
};
```

`escapeHtml` lives in `global.js:1158` and is loaded on every page.

The `kwargs.node_visible` branch is carried over from `detailsLink` unchanged - `grep -rn node_visible` to see
whether any node still sets it; if nothing does, drop the three lines.

### 3c. Expandable row

Three parts: the grid declares an expand renderer, the adapter passes it through, and a small view renders
the fragment. This reuses `setupClientExpend` (`datatable_definition.js:498`) exactly as the other
expandable grids do, with two small generic additions to it.

**Grid** (`snpdb/grids.py` `AbstractVariantGrid`):

```python
    def _get_annotation_version(self) -> AnnotationVersion:
        return self.annotation_version

    def get_expand_client_renderer(self) -> str:
        """ Row expansion (@see DataTableDefinition.setupClientExpend) - variantGridRowDetail in grid.js
            fetches variant_grid_row_detail for the annotation version this grid is showing """
        return f"variantGridRowDetail.bind(null, {self._get_annotation_version().pk})"
```

`analysis.grids.VariantGrid` overrides `_get_annotation_version` to return `self.node.analysis.annotation_version`
(every other subclass already sets `self.annotation_version` in `__init__`).

**Adapter** (`jqgrid_datatable_adapter.py` `datatable_definition`), after `"extra"`:

```python
    if expand_client_renderer := getattr(grid, "get_expand_client_renderer", lambda: None)():
        data["expandClientRenderer"] = expand_client_renderer
        # Rows are full of links; the arrow makes the affordance clear enough without a hover prefetch
        # request per row the mouse rests on
        data["expandPrefetch"] = False
```

Add both keys to the module docstring / `datatable_definition` docstring where the other keys are described.

**`DataTableDefinition.setupClientExpend`** (`datatable_definition.js:498-570`) - two generic changes:

```js
            dom.on('click', 'tr', function(event) {
                // A click on a link or control belongs to it (the variant grids have a details link, a
                // checkbox and tag links in every row) - only bare row space toggles
                if ($(event.target).closest('a, input, button, select, label').length) {
                    return;
                }
                const tr = $(this);
                ...unchanged...
            });
            if (this.serverParams.expandPrefetch !== false) {
                dom.on('mouseenter', 'tr', function() { ...unchanged prefetch... });
                dom.on('mouseleave', 'tr', function() { ...unchanged... });
            }
```

The click guard changes behaviour on every expandable grid (a click on a link no longer also toggles the
row). That is the intended behaviour everywhere; call it out in the commit message.

**JS** (`grid.js`, next to `load_variant_details`):

```js
// Expand renderer for the variant grids (@see AbstractVariantGrid.get_expand_client_renderer). The
// identifiers come from the server for the grid's annotation version; the actions are built here so
// IGV follows ANALYSIS_SETTINGS.show_igv_links like the rest of the page
function variantGridRowDetail(annotationVersionId, rowData) {
    const variantId = rowData["id"];
    const detail = $('<div>', {class: 'variant-row-detail'});
    const fragment = $('<div>', {class: 'variant-row-detail-fragment', text: 'Loading...'});
    loadAjaxBlock(fragment, Urls.variant_grid_row_detail(variantId, annotationVersionId));

    const actions = $('<div>', {class: 'variant-row-detail-actions'});
    actions.append($('<a>', {href: `javascript:load_variant_details(${variantId});`, text: 'Details'}));
    actions.append($('<a>', {href: Urls.view_variant(variantId), target: '_blank', text: 'Open in new tab'}));
    const chrom = rowData["locus__contig__name"];
    const position = rowData["locus__position"];
    if (chrom != null && position != null) {
        const igvUrl = createIgvUrl(`${chrom}:${position}`, 'getBams');  // null unless the analysis shows IGV links
        if (igvUrl) {
            actions.append($('<a>', {href: igvUrl, text: 'IGV'}));
        }
    }
    detail.append(fragment, actions);
    return detail;
}
```

`getBams` is the analysis page's function (`createIgvUrl` takes its *name* and appends `()`), the same
string `detailsLink` passed to `create_igv_link`.

**URL** (`variantopedia/urls.py`, next to `view_variant`):

```python
    path('variant_grid_row_detail/<int:variant_id>/<int:annotation_version_id>',
         views.variant_grid_row_detail, name='variant_grid_row_detail'),
```

`Urls.variant_grid_row_detail` is generated by django-js-reverse - run `python3 manage.py collectstatic_js_reverse`
after adding the URL (the `TableFormat.expandAjax` error text says the same).

**View** (`variantopedia/views.py`):

```python
def variant_grid_row_detail(request, variant_id: int, annotation_version_id: int):
    """ The expanded row under a variant grid row - identifiers for the grid's annotation version.
        @see variantGridRowDetail in grid.js """
    variant = get_object_or_404(Variant, pk=variant_id)
    annotation_version = get_object_or_404(AnnotationVersion, pk=annotation_version_id)
    vav = annotation_version.variant_annotation_version
    variant_annotation = VariantAnnotation.objects.filter(variant=variant, version=vav).first()
    transcript_annotations = VariantTranscriptAnnotation.objects.filter(variant=variant, version=vav) \
        .exclude(hgvs_c__isnull=True).order_by("transcript_version__transcript_id").select_related("transcript_version")
    context = {
        "variant": variant,
        "variant_annotation": variant_annotation,
        "transcript_annotations": transcript_annotations[:VARIANT_GRID_ROW_DETAIL_MAX_TRANSCRIPTS],
        "build_variants": variant.all_build_variants,
        "clingen_allele": variant.allele.clingen_allele if variant.allele else None,
    }
    return render(request, "variantopedia/variant_grid_row_detail.html", context)
```

Check `VariantTranscriptAnnotation` ordering fields against the model (`annotation/models/models.py:2396`) -
the intent is a stable order with the representative transcript first; if that needs a `Case`, sort in
Python over the (small) list instead. `VARIANT_GRID_ROW_DETAIL_MAX_TRANSCRIPTS = 20` as a module constant.

**Template** (`variantopedia/templates/variantopedia/variant_grid_row_detail.html`) - a two-column layout,
no page chrome:

```django
<div class="variant-row-detail-columns">
  <table class="variant-row-detail-table">
    {% if variant_annotation.hgvs_c %}<tr><th>Transcript</th><td class="mono">{{ variant_annotation.hgvs_c }}</td></tr>{% endif %}
    {% if variant_annotation.hgvs_p %}<tr><th>Protein</th><td class="mono">{{ variant_annotation.hgvs_p }}</td></tr>{% endif %}
    {% if variant_annotation.hgvs_g %}<tr><th>Genomic</th><td class="mono">{{ variant_annotation.hgvs_g }}</td></tr>{% endif %}
    {% for build_variant in build_variants %}
      <tr><th>{{ build_variant.any_genome_build }}</th><td class="mono">{{ build_variant.full_string }}</td></tr>
    {% endfor %}
    {% if clingen_allele %}<tr><th>ClinGen</th><td class="mono">{{ clingen_allele }}</td></tr>{% endif %}
    <tr><th>VariantGrid</th><td class="mono">v {{ variant.pk }}</td></tr>
  </table>
  {% if transcript_annotations %}
  <table class="variant-row-detail-table">
    <tr><th colspan="2">Transcripts</th></tr>
    {% for ta in transcript_annotations %}
      <tr><td class="mono">{{ ta.hgvs_c }}</td><td class="mono text-muted">{{ ta.hgvs_p|default:"" }}</td></tr>
    {% endfor %}
  </table>
  {% endif %}
</div>
```

Django auto-escapes `>` in HGVS. `variant_annotation` may be `None` (reference / unannotated variants) -
the `{% if %}` guards handle it; `build_variants` always has the variant itself.

### 3d. Styles

`global.scss` inside `table.variantgrid-datatable`, replacing the `.variant_id-container` block and its comment:

```scss
	// Representative variant cell: expand arrow (td.toggle-link:before) · checkbox · label (details link)
	.variant_id-container {
		display: inline-flex;
		align-items: center;
		gap: 4px;
		min-width: 0;
		max-width: calc(100% - 16px);  // the arrow's width + margin

		a.rv-label {
			overflow: hidden;
			text-overflow: ellipsis;
			white-space: nowrap;
			text-decoration: none;
		}
	}

	.rv-gene {
		font-weight: 600;
	}

	.rv-hgvs {
		font-family: ui-monospace, Menlo, Consolas, monospace;
		font-size: 11px;
	}

	.rv-hgvs-p,
	.rv-sub {
		color: #777;
		font-size: 11px;
	}

	// Only the arrow cell invites a click - the rest of the row is links
	&.expandable td.toggle-link {
		cursor: pointer;
	}

	// The child row's td spans every column of a table that can be thousands of px wide - keep the
	// detail pinned to the visible left edge of the scroll body
	.variant-row-detail {
		position: sticky;
		left: 0;
		width: max-content;
		max-width: 90vw;
		padding: 6px 10px;
		font-size: 12px;
		white-space: normal;

		.variant-row-detail-columns {
			display: flex;
			gap: 32px;
			align-items: flex-start;
		}

		.variant-row-detail-table th {
			color: #777;
			font-weight: 500;
			padding-right: 10px;
			white-space: nowrap;
			text-align: left;
		}

		.variant-row-detail-table .mono {
			font-family: ui-monospace, Menlo, Consolas, monospace;
			font-size: 11px;
		}

		.variant-row-detail-actions {
			margin-top: 6px;
			padding-top: 6px;
			border-top: 1px solid #e1e5ea;
			display: flex;
			gap: 12px;
		}
	}
```

`table.datatable.expandable tr.odd/even { cursor: pointer }` (`global.scss:1197`) targets `table.datatable`,
which the variant grids are not (`variantgrid-datatable`), so the row-wide pointer does not apply here.
`#node-data-container a.grid-link { text-align: center; font-weight: bold }` (`global.scss:892`) still
serves tag links; leave it.

## Step 4 - Genomic sort for the Variant column

`snpdb/grids.py` `AbstractVariantGrid`:

```python
from django.db.models import Case, IntegerField, Q, QuerySet, Value, When

    VARIANT_COLUMN_SIDX = "id"  # the mandatory Variant column's colmodel name/index

    def _genomic_order_by(self, items: QuerySet, descending: bool) -> QuerySet:
        """ The Variant column sorts in genome build order whatever it displays. Contig order is a CASE over
            the build's standard contigs rather than a join through GenomeBuildContig - MT is shared between
            builds and the join would return the row once per build. Non-standard contigs sort after, by id """
        whens = [When(locus__contig_id=contig_id, then=Value(i))
                 for i, contig_id in enumerate(self.genome_build.standard_contigs.values_list("pk", flat=True))]
        items = items.annotate(_contig_order=Case(*whens, default=Value(len(whens)), output_field=IntegerField()))
        fields = ["_contig_order", "locus__contig_id", "locus__position", "locus__ref__seq", "alt__seq", "pk"]
        prefix = "-" if descending else ""
        return items.order_by(*[f"{prefix}{f}" for f in fields])

    def _sort_items(self, items, sidx, sord):
        if sidx == self.VARIANT_COLUMN_SIDX:
            return self._genomic_order_by(items, descending=sord == "desc")
        return super()._sort_items(items, sidx, sord)
```

Check on each subclass:
- `analysis.grids.VariantGrid._sort_items` returns `super()._sort_items(items, None, sord)` when sorting is
  disabled and otherwise only intercepts `":"` packed genotype indexes, then calls `super()` → this branch. ✔
- `ExportVariantGrid.sort_items` returns items untouched (export order is set by `paginate_items`). ✔
- `AllVariantsGrid._sort_items` replaces the whole thing with `DEFAULT_ORDER_BY`. ✔ (nothing to do)
- `NearbyVariantsGrid`, `TaggedVariantGrid`, `GeneSymbolVariantsGrid` inherit this. All set `self.genome_build`
  before `super().__init__`. ✔
- `standard_contigs` is a `cached_property` queryset ordered by `genomebuildcontig__order` (`models_genome.py:145-153`);
  the `values_list` is one small query per sort request.
- `.annotate()` after `.values()` adds `_contig_order` to every row dict. The DataTables client ignores keys
  with no column; CSV export writes colmodel keys only. Harmless.

The analysis default sort (`Analysis.default_sort_by_column`, a `VariantGridColumn` FK) set to `variant`
gives `sortname = "id"`; `datatable_order_from_config` finds colmodel index 0 and the first page request
carries `order[0][column]=0` → `sidx=id` → this branch. The `sorting_disabled` guard still applies.

## Step 5 - Migration, cache, deployment check

`snpdb/migrations/<next>_classifications_column_and_representative_variant.py` (check `ls snpdb/migrations | tail -1`
for the number and the dependency):

```python
from django.db import migrations
from django.db.models import F

from library.django_utils import bulk_insert_class_data

# Raw coordinate columns the Variant column now carries - they stay in the catalogue and in 'All columns'
_REPLACED_BY_VARIANT_COLUMN = ["chrom", "position", "ref", "alt", "hgvs_g"]
_KEEP_EVERYTHING = "All columns"

_CLASSIFICATIONS_COLUMN = {
    'grid_column_name': 'classifications',
    'variant_column': 'classifications',
    'annotation_level': 'D',
    'width': 170,
    'label': 'Classifications',
    'description': 'ClinVar, internal germline and internal somatic classifications. Hover a chip for the '
                   'per-record summary, click it to scroll to the underlying column',
    'model_field': False,
    'queryset_field': False,
}


def _add_columns(apps, _schema_editor):
    bulk_insert_class_data(apps, "snpdb", [("VariantGridColumn", [_CLASSIFICATIONS_COLUMN])])
    VariantGridColumn = apps.get_model("snpdb", "VariantGridColumn")
    VariantGridColumn.objects.filter(pk="variant").update(
        width=280,
        description="The variant named by c.HGVS, g.HGVS or coordinate - whichever is available. "
                    "Click the name for details, the row for transcript / HGVS / coordinates / IGV. "
                    "Sorts in genomic order")


def _update_system_collections(apps, _schema_editor):
    CustomColumnsCollection = apps.get_model("snpdb", "CustomColumnsCollection")
    CustomColumn = apps.get_model("snpdb", "CustomColumn")

    for ccc in CustomColumnsCollection.objects.filter(user__isnull=True):
        cc_qs = CustomColumn.objects.filter(custom_columns_collection=ccc)
        variant_cc = cc_qs.get(column_id="variant")
        cc_qs.filter(sort_order__gt=variant_cc.sort_order).update(sort_order=F("sort_order") + 1)
        CustomColumn.objects.create(custom_columns_collection=ccc, column_id="classifications",
                                    sort_order=variant_cc.sort_order + 1)
        if ccc.name != _KEEP_EVERYTHING:
            cc_qs.filter(column_id__in=_REPLACED_BY_VARIANT_COLUMN).delete()
        # Historical models skip CustomColumn.save()/delete(), which is what normally bumps the version
        # the node grid definition cache is keyed on
        ccc.version_id += 1
        ccc.save()


def _remove_columns(apps, _schema_editor):
    VariantGridColumn = apps.get_model("snpdb", "VariantGridColumn")
    VariantGridColumn.objects.filter(pk="classifications").delete()  # cascades to CustomColumn


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '<previous>'),
    ]

    operations = [
        migrations.RunPython(_add_columns, reverse_code=_remove_columns),
        migrations.RunPython(_update_system_collections, reverse_code=migrations.RunPython.noop),
    ]
```

Notes:
- `CustomColumn` has `unique_together (collection, column)`; `bulk_insert_class_data` + `create` are fine
  because the column is new.
- A fresh test database runs `0002_initial_data` (which still seeds `chrom, position, ...`) and then this
  migration, so tests see the final shape. Leave `0002` alone.
- `settings.DEFAULT_COLUMNS_NAME` is `'Default columns'`; the filter is on `user__isnull=True` rather than
  the name so `Minimalism` is covered too.

`variantgrid/settings/components/default_settings.py`: `CACHE_VERSION = 47` → `48`. Without this, node grid
definitions cached up to a week ago for **user-cloned** collections (whose `ccc_version_id` did not change)
still name `VariantGridFormat.detailsLink`, and `eval` of that in `datatable_definition.js:315` throws.

`variantgrid/deployment_validation/column_check.py`: add `"classifications"` to `SPECIAL_COLUMNS` (Step 2a).

## Step 6 - Tests

Python (`analysis/tests/test_representative_variant_column.py`, extend `GridExportTestCase` from
`analysis/tests/test_grid_export.py` - it already creates variants on contigs `3, 1, 10, 1, 2, 10` plus one on
`2`, deliberately out of both pk and string order):

```python
from analysis.grids import VariantGrid
from analysis.tests.test_grid_export import GridExportTestCase
from library.django_utils import FakeRequest
from library.django_utils.jqgrid_datatable_adapter import datatable_definition
from snpdb.grid_columns.custom_columns import (
    CLASSIFICATIONS_COLUMN_ROW_ANNOTATIONS, CLASSIFICATIONS_COLUMN_ROW_FIELDS, VARIANT_COLUMN_ROW_FIELDS,
)
from snpdb.models import CustomColumnsCollection


class RepresentativeVariantColumnTest(GridExportTestCase):
    def setUp(self):
        super().setUp()
        self.node = self._sample_node()
        self.grid = VariantGrid(self.user, self.node)

    def _colmodels_by_name(self) -> dict:
        return {cm["name"]: cm for cm in self.grid.get_colmodels()}

    def test_variant_column_sorts_in_genome_build_order(self):
        qs = self.grid._sort_items(self.node.get_queryset(), "id", "asc")
        coords = [(v.locus.contig.name, v.locus.position) for v in qs]
        self.assertEqual(coords, [("1", 1000), ("1", 2000), ("2", 5000), ("2", 9000),
                                  ("3", 3000), ("10", 500), ("10", 1000)])

    def test_variant_column_sort_descending(self):
        qs = self.grid._sort_items(self.node.get_queryset(), "id", "desc")
        self.assertEqual([(v.locus.contig.name, v.locus.position) for v in qs][:2], [("10", 1000), ("10", 500)])

    def test_renderer_fields_ride_along_hidden(self):
        colmodels = self._colmodels_by_name()
        self.assertEqual(colmodels["id"]["formatter"], "representativeVariant")
        for field in VARIANT_COLUMN_ROW_FIELDS + CLASSIFICATIONS_COLUMN_ROW_FIELDS + CLASSIFICATIONS_COLUMN_ROW_ANNOTATIONS:
            self.assertIn(field, colmodels)
        # Not in the default collection, so hidden, and labelled from the catalogue for the CSV header
        self.assertTrue(colmodels["locus__ref__seq"]["hidden"])
        self.assertEqual(colmodels["locus__ref__seq"]["label"], "Reference")
        self.assertEqual(colmodels["alt__seq"]["label"], "Alt")

    def test_rows_carry_renderer_fields(self):
        request = FakeRequest(user=self.user)
        request.GET = {"rows": "1", "page": "1"}
        row = self.grid.get_data(request)["rows"][0]
        for field in VARIANT_COLUMN_ROW_FIELDS + CLASSIFICATIONS_COLUMN_ROW_FIELDS:
            self.assertIn(field, row)
        self.assertNotIn("classifications", row)  # renderer-only column, like tags

    def test_classifications_column_is_not_sortable(self):
        cm = self._colmodels_by_name()["classifications"]
        self.assertEqual(cm["formatter"], "classificationsFormatter")
        self.assertFalse(cm["sortable"])

    def test_definition_declares_row_expansion(self):
        definition = datatable_definition(self.grid)
        annotation_version_id = self.node.analysis.annotation_version_id
        self.assertEqual(definition["expandClientRenderer"], f"variantGridRowDetail.bind(null, {annotation_version_id})")
        self.assertFalse(definition["expandPrefetch"])

    def test_system_default_collection(self):
        columns = list(CustomColumnsCollection.get_system_default().customcolumn_set
                       .order_by("sort_order").values_list("column_id", flat=True))
        self.assertEqual(columns[:2], ["variant", "classifications"])
        for removed in ["chrom", "position", "ref", "alt", "hgvs_g"]:
            self.assertNotIn(removed, columns)


class VariantGridRowDetailViewTest(GridExportTestCase):
    def test_row_detail_renders(self):
        client = Client()
        client.force_login(self.user)
        variant = self.variants[0]
        url = reverse("variant_grid_row_detail", kwargs={"variant_id": variant.pk,
                                                         "annotation_version_id": self.analysis.annotation_version_id})
        response = client.get(url)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, f"v {variant.pk}")
        self.assertContains(response, variant.full_string)
```

`test_rows_carry_renderer_fields` may need the same setup the other `GridExportTestCase` tests use - copy from
`AdapterParamTranslationTest` in `analysis/tests/test_jqgrid_datatable_adapter.py` if `get_data` complains.
The view test's variants have no `VariantAnnotation` (fake annotation version only), which is exactly the
"unannotated" branch of the template - good, that is the path most likely to be wrong.

Existing tests to update/run:
- `analysis/tests/test_jqgrid_datatable_adapter.py::AdapterColumnsTest.test_column_shape` uses the
  `detailsLink` name literally - switch it to `representativeVariant` / `VariantGridFormat.representativeVariant`.
- `analysis/tests/test_grid_export.py` - the CSV header gains `Reference`, `Alt`, `SVLEN`, `Symbol`, `c.HGVS`,
  `p.HGVS`, `HGVS genomic`, `Classifications` columns where they weren't in the collection. Its tests look
  columns up by name so they should pass; run the module.
- `analysis/tests/test_node_grid_sort_limit.py`, `test_node_grid_datatable_endpoints.py`,
  `snpdb/tests/test_vcf_export_columns.py`, and the deployment column check
  (`python3 manage.py deployment_check` or the function directly) - all must stay green.
- A URL test for `variant_grid_row_detail` in the variantopedia `URLTestCase` module if one lists the variant URLs.

```bash
python3 manage.py test --keepdb analysis.tests.test_representative_variant_column analysis.tests.test_jqgrid_datatable_adapter analysis.tests.test_grid_export analysis.tests.test_grid_export_vcf analysis.tests.test_node_grid_sort_limit analysis.tests.test_node_grid_datatable_endpoints snpdb.tests.test_vcf_export_columns
```

There is no JS test runner (`package.json` only has eslint); run
`npx eslint variantgrid/static_files/default_static/js/variantgrid_formats.js variantgrid/static_files/default_static/js/grid.js variantgrid/static_files/default_static/js/datatable_definition.js`
and do the manual checks below.

## Step 7 - Cleanup

- Delete `VariantGridFormat.detailsLink`, `internalClassificationBox` and the `detailsLink` adapter entry.
  `GERMLINE_CLASSIFICATION_BOXES` / `SOMATIC_CLASSIFICATION_BOXES` stay (chip tooltips use `.label`).
- `createGridLink` / `create_igv_link` in `grid.js` are still used by other renderers (tags, gene filter link);
  leave them.
- Update the `.view-details-link` / `.igv-link` CSS only if nothing else references them (`grep -rn` first -
  `igv-link` is used on the variant details page).
- Update `claude/research/analysis.md` if it describes the `detailsLink` boxes.
- `pylint` via `./scripts/linting/run_pylint.sh` on touched Python files.

## Things that could go wrong

1. **Stale node grid caches.** Rows and definitions are `cache_page`d for a week. Rows lacking the new keys
   must render (the cascade and chips handle `undefined`); definitions naming `detailsLink` must not be
   served → bump `CACHE_VERSION`. Test by loading an analysis whose grid was viewed before the change.
2. **Migration versions.** Queryset `.delete()`/`create()` on historical `CustomColumn` bypass the model's
   `increment_version()`; without the explicit `version_id += 1` the node grid config URL keeps the old
   `ccc_version_id` and hits the stale cache even after the cache bump has expired. Keep the explicit bump.
3. **`hgvs_c` shapes.** `NM_x(SYMBOL):c.` and `NM_x:c.` both occur; fusions have no `accession:` at all;
   the SV sentinel string is not an HGVS. The regex + sentinel check cover all four - keep a row of each in
   the manual checks.
4. **Long inserted sequences in `hgvs_g`** (`g.123_124insACGT…` with the full sequence) blow the cell width -
   hence `REPRESENTATIVE_MAX_HGVS_CHARS` falls through to the coordinate step, which collapses alleles to `[Nbp]`.
5. **Row click vs. the row's own controls.** Without the `closest('a, input, ...')` guard in
   `setupClientExpend`, ticking the select checkbox or clicking a tag also expands the row. The guard is
   generic and changes every expandable grid (for the better); watch the eventlog / classification grouping
   grids still expand on a bare-row click.
6. **`odd`/`even` row classes.** `setupClientExpend` only toggles rows carrying them. DataTables adds them by
   default (`stripeClasses`), but the variant grids use their own `tableClass`; if rows come out without
   `odd`/`even`, pass `stripeClasses: ['odd', 'even']` in `dtParams` (`datatable_definition.js:214`).
7. **Child row width.** The child `td` spans every column of a `table-layout: fixed` table that can be
   thousands of px wide inside `.dataTables_scrollBody`. `.variant-row-detail { position: sticky; left: 0 }`
   keeps the content at the visible left edge - verify by scrolling right with a row open. If sticky misbehaves
   inside the DataTables scroll body, fall back to `width: <visible width>` set from `$('.dataTables_scrollBody').width()`
   in `variantGridRowDetail`.
8. **Hover prefetch.** Left on, every row the mouse rests on for 500ms would fetch the detail view. The
   adapter sets `expandPrefetch: false` for these grids; keep it that way.
9. **Annotation version in the URL.** `variantGridRowDetail` is bound to the grid's annotation version at
   definition time, so an analysis on an older annotation version shows the HGVS it was analysed with (the
   same version the grid columns show). `get_object_or_404(AnnotationVersion)` guards a bad id.
10. **`sortable: False` on `classifications`.** Without it `F("classifications")` raises `FieldError` (caught
    and logged as a warning, and the sort silently becomes `-pk`).
11. **Deployment check** flags `classifications` unless it is in `SPECIAL_COLUMNS`.
12. **CSV header labels.** Hidden model fields are exported; without the catalogue label lookup the header has
    two `seq` columns. `variant_id` in `GridExportTestCase` tests is the `id` column's *label override* from
    the export path - leave that alone.
13. **Contig order.** Sorting by `locus__contig_id` alone puts GRCh38 `MT` first. The `CASE` over
    `standard_contigs` is the fix; the `locus__contig_id` tiebreaker only orders the non-standard contigs.
14. **`sorting_disabled` nodes.** `VariantGrid._sort_items` returns before reaching the genomic branch on
    big nodes (`ANALYSIS_GRID_SORT_MAX_ROWS`) - expected, and the "sorting disabled" banner explains it.
15. **Right-click "open in new tab"** on the label relies on `a.variant-link` + `.variant_id-container[variant_id]`
    (`grid.js:227-240`). Keep both class names exactly.
16. **`showGridCell` targets a hidden column** for the chips when the underlying column isn't in the
    collection - same as today's boxes (it scrolls nowhere). Acceptable; the tooltip carries the summary.
17. **Width.** `table-layout: fixed` - the `id` colmodel width is the width the cell gets. 280 fits
    `GBA1 c.1226A>G p.(Asn409Ser)` plus the arrow at 11-12px; longer p.HGVS (frameshifts) ellipsize, which is
    fine because the tooltip and expanded row hold the full strings.
18. **`Urls.variant_grid_row_detail` missing** until `collectstatic_js_reverse` runs - the page then shows
    "Method or URL not configured" style errors from `loadAjaxBlock`. Part of the deploy, not the code.
19. **`escapeHtml` on the `title` attribute** - HGVS contains `>`; single quotes in gene-level nomenclature
    are escaped as `&#39;`. Always route through `escapeHtml`, never string-concatenate raw values.
20. **Dual-screen analysis mode** - `getAnalysisWindow()` is the *other* window for some helpers;
    `variantGridRowDetail` only uses `Urls`, `loadAjaxBlock` and `createIgvUrl` (which already goes through
    `getAnalysisWindow`), so nothing new crosses windows.

## Manual verification

Analysis page, a node with a few hundred rows, `Default columns`:
- Leading columns are `Variant` (with the expand arrow), `Classifications`, `Tags`; no chrom/position/ref/alt/hgvs_g.
- SNV row reads `GENE c.… p.(…)`; hover shows full HGVS; click loads details in the editor pane; right-click →
  open in new tab lands on the variant page.
- Clicking bare row space (or the arrow) expands the row: transcript/protein/genomic, coordinate per build,
  ClinGen, `v id`, the transcript list, and Details / Open in new tab / IGV (IGV only when `show_igv_links`).
  Clicking another row closes the first; clicking the same row closes it; the checkbox and tag links do
  **not** expand the row; scrolling right keeps the detail visible; paging closes it.
- Rows for: an intergenic variant (g.HGVS step), a `<DEL>`/`<DUP>` (span + type + size), a long insertion
  (`[Nbp]`), a variant with the SV sentinel in `hgvs_c`, a reference (`=`) variant, an unannotated variant
  (`v id`; its expanded row shows just coordinates and id).
- Header click on `Variant` sorts `1, 2, …, 22, X, Y, MT` then position; second click reverses.
- Classifications: ClinVar chip with stars, `×n` germline count, somatic tier chip, dashed placeholders;
  tooltip text matches what the old boxes said.
- Checkbox selection still marks the node dirty and survives paging (`selectedVariants`).
- CSV export of the node: header sane, rows unchanged apart from the new columns.
- Variantopedia "All variants", a "Variants with tags" page, a gene symbol page grid and a nearby-variants
  tab all render the new cell, chips and expansion (no checkbox, no IGV action outside an analysis).
- `Minimalism` and a user-cloned collection (still containing chrom/position…) both render.
- One other expandable grid (e.g. the eventlog) still expands on row click and its links still work.

## Follow-ups (separate issues)

- Sample zygosity glyph column with header sort menu (zygosity / AF / depth)
- Merge gnomAD AF and gnomAD filtered into one column, the way popmax AF and popmax are merged -
  `variantannotation__gnomad_af` carries the value, `variantannotation__gnomad_filtered` rides along
  hidden for the filtered marker (it has its own `gnomadFilteredFormatter` today)
- gnomAD "common" AF muting (the value + popmax pairing is built - see "As built")
- Predictions meter, SpliceAI threshold dot
- Two-line rows as a per-user toggle
- ClinVar somatic / oncogenic chips
- Expanded row: overlapping genes for SVs, links to the classification records behind the chips

## As built

Steps 1-7 are implemented and tested. Where the code differs from the steps above, the code is right and
this section says why.

### Extra scope, added during implementation

- **Leading cluster is `variant, classifications, tags, tags_global`.** The migration reorders every
  system collection so the two tag columns sit with the classifications, matching the mockup's grid.
- **`svlen` leaves the default collections** with the coordinate fields. It rides along hidden, so the
  Variant cell still sizes symbolic variants and the CSV/VCF exports keep their SVLEN column.
- **Impact · Consequence is one cell** - `impactConsequenceFormatter` on `variantannotation__consequence`
  draws a coloured impact dot before the consequence text, reading `variantannotation__impact` from the
  hidden ride-along. The column sorts by consequence; `impact` stays in the catalogue and `All columns`.
- **gnomAD PopMax is one cell** - `gnomadPopmaxFormatter` on `variantannotation__gnomad_popmax_af` puts
  the population beside the frequency, reading `variantannotation__gnomad_popmax` from the hidden
  ride-along. The frequency keeps its server side unit->percent formatting so the CSV matches the grid,
  so the renderer only appends the population. The column sorts by frequency. Pairing the two popmax
  fields is the confirmed intent, in place of the mockup card's `gnomad_af` + `gnomad_popmax`.
- **`_get_standard_overrides` merges `af_override` per key** rather than replacing, so a column carrying
  both a client formatter and a server side one (`gnomad_popmax_af`) keeps both.

### Details the steps left open

- **`HGVS_REGEX` requires an HGVS kind after the colon** (`/^([^:(]+)(?:\(([^)]+)\))?:([cgmnopr]\..+)$/`).
  Gene-level fusion nomenclature is `BCR::ABL1` (@see `fusion_canonical_str` in `genes/models.py`), which a
  looser pattern splits into accession `BCR` + change `:ABL1`. With the kind required, fusions fall to the
  "show it whole" branch.
- **`representativeVariant` reads node visibility from `ctx.extra`** - `formatter_kwargs` has no server side
  writer, so the cascade through `_isNodeVisible(ctx)` is the whole story.
- **The row detail orders transcripts by `hgvs_c` in SQL**, slices, then hoists the representative
  transcript to the front in Python (a stable order before the slice, representative first after it).
- **Classification chips sit in a `.cs-chips` grid** with proportional slots
  (`grid-template-columns: 74fr 52fr 60fr`). The three origins line up down the column at whatever width
  the column is set to. Shares rather than fixed pixels: the cell is `overflow: hidden`, so fixed widths
  wider than the column silently cut the last chip off the right edge.
- **`_update_system_collections` is idempotent** - it recomputes the order from the collection's current
  state and creates the column only when absent, so a reverse-then-forward converges on a database that
  already ran an earlier revision of this migration.
- **A `variant_grid_row_detail` case was added to `variantopedia/tests/test_urls.py`.** That fixture has a
  fully annotated variant, so it covers the template's annotated branch while
  `analysis/tests/test_representative_variant_column.py` covers the unannotated one.

### State to know about

- Everything is edits only; nothing is committed. The working tree also carries a separate in-progress
  ClinVar-node change (`clinvar_node.py`, `forms_nodes.py`, `node_display.py`, `clinvarnode_editor.html`,
  `test_clinvar_node.py`, `analysis_nodes.css/js`, `analysis/migrations/0127_*`) that belongs to the user -
  a commit here selects only the files this plan names.
- The dev database has run migration `snpdb/0225` at its current content (reversed and re-applied).
- `analysis.tests.test_serializers.test_node_counts_round_trip` fails on clean `master` as well -
  `analysis/analysis_import_export.py` reads `NodeCountSettings.built_in_filter`, which snpdb migration
  0224 removed. Unrelated to this work.
- `create_igv_link` in `grid.js` has no callers now that `detailsLink` is gone.

### Verified

Renderers were exercised under node against real GRCh38 rows for every `hgvs_c` shape, the symbolic and
long-insertion branches, rows cached before the new fields existed, and every chip combination. The genomic
sort was run against a real 1086-row GRCh38 node (ascending starts at chr1, descending at MT). The row
detail view was rendered for a real annotated variant. `python3 manage.py test --keepdb` over the modules in
Step 6 plus `variantopedia.tests.test_urls`: 77 tests, OK. eslint: 0 errors. pylint: 9.96/10.

Checked in a browser against the "Manual verification" list above and working.
