# Convert the variant grids (AbstractVariantGrid family) from jqGrid to DataTables

## Scope

The last jqGrid grids that render *variants* with user-configurable columns:

| Grid class | File | Page(s) |
|---|---|---|
| `VariantGrid` | `analysis/grids.py:75` | Analysis node grid (also reused by the Sample page `sample_variants_tab.html`) |
| `AllVariantsGrid` | `variantopedia/grids.py:75` | Variants page (`variants.html`) |
| `NearbyVariantsGrid` | `variantopedia/grids.py:187` | Nearby variants page (5 grid instances) |
| `TaggedVariantGrid` | `variantopedia/grids.py:301` | Variant tags page (variant-centric) |
| ~~`VariantTagsGrid`~~ | ~~`variantopedia/grids.py:219`~~ | **DONE (Stage 1)** — now `VariantTagsColumns` in `variantopedia/grids.py`, serving the variant tags page + `base_related_analyses.html` |
| `GeneSymbolVariantsGrid` | `genes/grids.py:125` | Gene symbol page |

Staying as they are (converted later, separately): `ExportVariantGrid` (server-only, powers Celery CSV/VCF export — no browser involvement), `NodeColumnSummaryGrid` / `NodeOntologyGenesGrid` / `NodeGeneDiseaseClassificationGenesGrid` (DataFrameJqGrid popups), `AbstractSkippedAnnotationGrid`, `CohortSampleListGrid`, `GenesGrid`, `QCGeneCoverageGrid`, and the seqauto/pathtests/upload jqGrids. `library/jqgrid/` server code also stays: 16+ other grids use it, and `FakeFilterGrid` (`analysis/models/nodes/filters/filter_node.py:26`) uses its filter-rule→Q parser as FilterNode's persistence format.

## Recommendation: staged, using an adapter — the jqGrid *server classes stay as the engine*, the *client* converts to DataTables

### Why staged

- The client surface is huge and disjoint: 6 pages, ~16 custom formatters, the whole analysis editor wiring (selection, tagging, IGV, filter-child links, deferred load, resize, component registry), plus three different export paths. One PR touching all of it cannot be validated incrementally and would sit unmergeable for weeks.
- The pages have independent risk profiles. Nearby-variants is a read-only page with five plain grids; the analysis node grid carries week-long response caching, per-user locks, deferred loading and the FilterNode editor. Landing the simple pages first proves the shared layer (adapter + ported formatters) with real traffic before the analysis work starts.
- Each stage leaves the app fully working — the jqGrid and DataTables systems already coexist everywhere else in the codebase.

### Why an adapter rather than a native `DatatableConfig` rewrite

The server classes encode years of hard-won performance behaviour, each with regression tests tied to their current API:

- `get_known_count()` — stored `NodeCount` skips `COUNT(*)` (`analysis/grids.py:180`), `EXPLAIN`-based approximate counts on All Variants (`variantopedia/grids.py:153`) — issues #1279, #1722
- `sorting_disabled()` + packed-genotype sort decoding (`analysis/grids.py:141,373`) — issue #1651
- Forced genomic-order pagination on All Variants (`variantopedia/grids.py:147`) — issue #1663
- `variantallele` join restricted to the grid's build (`snpdb/grids.py:552`) — issue #1626
- Server-side formatters shared verbatim with CSV/VCF export (`grid_export.py` consumes `get_colmodels()` + `get_items()` + `iter_format_items()`)

The current `DatatableConfig`/`DatabaseTableView` has none of the needed hooks (it always runs `qs.count()`, caps pages at 100, has fixed per-construction column sets, no server-side export, no known/approximate counts). Rebuilding all of the above in DataTables idiom duplicates the logic *and* the tests, for zero user-visible benefit — the user-facing problem is the jqGrid **JavaScript**: a dead upstream library, inconsistent UI, client-side CSV hacks.

So: a thin protocol adapter serves the existing grid classes to the standard DataTables client. Column building, querysets, sorting, counting, filtering and formatting all stay where they are, fully covered by existing tests. `grid_export.py`, `NodeColumnSummaryGrid`, `node_debug`, and the node editor forms keep consuming `VariantGrid` unchanged.

A later native rewrite (folding the engine into `DatatableConfig` proper and retiring `library/jqgrid`) remains possible as a pure server-side refactor once every jqGrid page is converted — it becomes invisible to users and is much safer done then.

---

## Stage 1 — VariantTagsGrid: straight native `DatatableConfig` conversion — **DONE**

First cab off the rank — a self-contained conversion following the established pattern (same shape as the already-converted `VariantTagCountsColumns` / `VariantTagDetailColumns` in the same file), which also forces the first ports of the tag/variant client renderers that later stages reuse.

1. **`VariantTagsColumns(DatatableConfig[VariantTag])`** in `variantopedia/grids.py`:
   - `get_initial_queryset`: `VariantTag.get_for_build(genome_build)` (build from query param), restricted through `allele__variantallele__genome_build` and annotated with `variant_string` via `Variant.annotate_variant_string(..., path_to_variant="allele__variantallele__variant__")` — same query the jqGrid class builds today (`variantopedia/grids.py:238-281`).
   - `filter_queryset`: the current `extra_filters` keys become plain query params (`analysis_ids` via `Analysis.filter_for_user`, `gene`, `tag`, `tags`, `user`), plus the `UserGridConfig.get(user, 'Variant Tags').show_group_data` logic with the explicit-user override.
   - Columns: `variant_string` (server renderer packages `variant__id` + genome build; client renderer is the ported `formatVariantTagFirstColumn` incl. the "New Classification" button for `RequiresClassification`), genome build, gene symbol (ported `geneSymbolNewWindowLink`), tag pill (ported `formatVariantTag` — reads the tag-colour CSS already emitted by `{% render_tag_styles_and_formatter %}`), analysis link (ported `formatAnalysis`), username, created timestamp, and a delete action column via `render_delete` → `group_permissions_object_delete` (replacing jqGrid `delete_row`). Default sort `variant_string` asc.
   - VariantTag is ~412k rows at SA Path: keep querysets subquery-shaped as they are today and leave the CSV on the existing server-side `variant_tags_export` streaming view (toolbar button via the `data-toolbar` mechanism) — the DataTables client CSV button fetches rows through the 100-row page cap and is wrong for this table.
2. **Templates**: `variant_tags.html` and `base_related_analyses.html` swap the `{% jqgrid %}` tag for `<table data-datatable-url=...>` with a `data-datatable-data` function carrying the current tag/user/analysis filter state (replacing `jqgrid_config_get_parameters_func` / `init_func`); tag-click and show-all handlers call `table.ajax.reload()`. Those handlers currently drive *both* grids on the page, so during this stage they fan out to both systems: `table.ajax.reload()` for the new DataTable and the existing `trigger("reloadGrid")` for the still-jqGrid `TaggedVariantGrid`, with the shared filter state (tag/user/show-all) held in one place both reads pull from (the `data-datatable-data` function and `getVariantTagsConfigUrl` respectively). Slightly fiddly, low risk; the dual path collapses in Stage 3.
3. **Wiring**: `DatabaseTableView.as_view(column_class=VariantTagsColumns)` URL replaces the `variant_tags_grid` JQGridView entry. `variant_tags_export` (`variantopedia/views.py:918`) currently subclasses the jqGrid class, so it moves onto the new config's queryset: a small streaming CSV view (chunked like `grid_export_csv`) over the same annotated queryset with an explicit `order_by("variant_string")`.
4. **Tests**: port the `VariantTagsGrid` assertions in `variantopedia/tests/test_tagged_variant_grid.py` (`_tags_grid_variant_ids`, user filter, show_group_data override) to the config; switch the `test_urls.py` entry to the datatable harness.
5. **Delete `VariantTagsGrid`** and its URL/imports once both templates are converted.

The variant tags page then runs one DataTable (tag-centric) beside one jqGrid (`TaggedVariantGrid`, variant-centric) until Stage 3 converts the rest of the page — the two systems already coexist on other pages.

### Stage 1 outcome — what actually landed

- `VariantTagsColumns(DatatableConfig[VariantTag])` in `variantopedia/grids.py`, `variant_tags_datatable` URL, `VariantTagsGrid` and its `variant_tags_grid` URL deleted.
- `variant_tags_export` (`variantopedia/views.py`) is now a `StreamingHttpResponse` over `filter_queryset(get_initial_queryset()).order_by("variant_string")`, chunked via `StashFile` + `EXPORT_ROWS_PER_CHUNK`, reached from a `data-toolbar` CSV button on both pages.
- Ported client renderers live in `datatables_client_renderers.js`: `renderVariantTagVariant`, `renderVariantTagAnalysis`, `renderVariantTagPill`, `renderGeneSymbolNewWindow`. URLs are built server-side through `url_if_visible` (renamed from `_url_if_visible` in `snpdb/grids.py`, now shared cross-app), so renderers stay free of `url_name_visible` guards.
- `variant_tags.html` holds the shared `gridParams` filter state that both systems read: `variantTagsDatatableFilter` (DataTables `data` hook, JSON-encoding list values) and `getVariantTagsGridExtraFilters` (jqGrid `extra_filters`); `reloadGrid()` fans out to `trigger("reloadGrid")` + `.DataTable().ajax.reload()`. This collapses in Stage 3.
- Dead jqGrid formatters removed: `formatVariantTag`/`formatAnalysis`/`formatVariantString` (`render_tag_styles_and_formatter.html`, now styles-only — the tag name is stale, rename it in Stage 5), `formatVariantTagFirstColumn` + `classifyAndCloseButton` (`grid.js`), and the matching entry in `variant_details_link_grid.html`.
- Tests: ported assertions plus row-payload / classify-URL / streaming-CSV coverage in `variantopedia/tests/test_tagged_variant_grid.py`; `test_urls.py` switched to `_test_datatables_grid_urls_contains_objs`.

**Carried into Stage 2:** `DatabaseTableView` runs both a `recordsTotal` and a `recordsFiltered` count per draw, where the jqGrid ran one — on VariantTag (412k rows at SA Path) that is the conversion's main cost. The adapter's known/approximate-count work is the place to give `DatatableConfig` a way to skip or supply the unfiltered count.

### Stage 1 reference — where everything lives

- jqGrid formatters to port: `formatVariantTag` + `formatAnalysis` + `formatVariantString` in `analysis/templates/analysis/tags/render_tag_styles_and_formatter.html:15-35` (uses `getVariantTagHtml`); `formatVariantTagFirstColumn` in `variantgrid/static_files/default_static/js/grid.js:649-658` (helper `classifyAndCloseButton` at `:643`); `geneSymbolNewWindowLink` at `grid.js:458-460` (impl `_geneSymbolLink` `:430-451`). Ported versions go in `variantgrid/static_files/default_static/js/datatables_client_renderers.js` with DataTables `(data, type, row)` signature; `TableFormat.*` examples in `datatable_definition.js:387+`.
- `{% jqgrid %}` usages to replace: `variantopedia/templates/variantopedia/variant_tags.html:159` (plus its `init_func`, `getVariantTagsConfigUrl`, CSV navButtons at `:20-31,51-58`, tagClick/showAll handlers) and `analysis/templates/analysis/tags/base_related_analyses.html:180` (`getRelatedVariantTagsConfigUrl`, `downloadNodeCSV`).
- Auto-wiring: a `<table data-datatable-url=... data-datatable-data=...>` is picked up by `global.js:347-358`; toolbar elements move in via `[data-toolbar="#tableId"]` (`datatable_definition.js` `setupDom`).
- Old URL entry: `variantopedia/urls.py:89-90` (`variant_tags_grid`, `JQGridView` with `delete_row=True`); export view `variant_tags_export` at `variantopedia/views.py:908-925` (uses `JQGridView.export_grid_as_csv`; chunked CSV writer `grid_export_csv` in `library/django_utils/jqgrid_view.py:139-158` is reusable for the new streaming export view).
- Delete-column pattern: `DatatableConfig.render_delete` + `TableFormat.deleteRow` client renderer, as used by e.g. `CohortListColumns` (`snpdb/grids.py`).
- Existing tests to port: `variantopedia/tests/test_tagged_variant_grid.py` (`_tags_grid_variant_ids` at L60, user-filter L79, show_group_data L108); URL test entry in `variantopedia/tests/test_urls.py:129-141`. DataTables test patterns: `seqauto/tests/test_sequencing_run_datatable.py` (config unit test via `RequestFactory` + `resolver_match`), `_test_datatables_grid_urls_contains_objs` in `library/django_utils/unittest_utils.py:196`.
- Tag colour CSS comes from `{% render_tag_styles_and_formatter %}` already on both pages — the ported renderers keep emitting `<span class='grid-tag tagged-<tag>'>` markup so it applies unchanged.

## Stage 2 — Adapter infrastructure + ported client formatters — **DONE**

New server pieces (suggested: `library/django_utils/jqgrid_datatable_adapter.py`):

1. **`JqGridDatatableView(View)`** — replaces `JQGridView` for converted grids, keeping its constructor conventions (`create_grid_from_request`: URL kwargs + `extra_filters` JSON → grid ctor).
   - `?dataTableDefinition=1` → definition JSON built from `grid.get_config(as_json=False)` / `get_colmodels()`:
     - per column: `data` = colmodel `name`, `label`, `orderable` = `sortable` (default True), `visible` = `not hidden` (respecting `hide_non_admin`), `className` from `classes`/width hints, `render` = mapped client formatter name (see formatter table below), `headerTitle` passthrough.
     - table level: `pageLength` from `UserGridConfig` rows (jqGrid `rowNum`), `order` from `sortname`/`sortorder`, `downloadUrl` (see exports), and an `extra` dict for grid-specific metadata (`analysisNode.visible`, `create_filter_child_links`, …).
   - Data requests → translate DataTables params to the jqGrid engine: `start`/`length` → `page`/`rows`; `order[0][column]` index → that colmodel's `index` (this is what makes packed-genotype sort work for free — the index string `"cohortgenotype_N:i:column"` round-trips exactly as today) → `sidx`/`sord`; optional `filters` param passed through untouched. Call `grid.get_items(request)` (wrapping params, e.g. via a lightweight request shim or by accepting the translated params directly), respond `{draw, recordsTotal, recordsFiltered, data}` with both records values from the paginator count (known/approximate counts flow through). Pass through `approximate_records` for pager display.
   - `sanitize`/`limit_value_size` equivalents from `DatabaseTableView.prepare_results` where sensible; rows keep their `.values()` field-name keys so column `data` keys match.
2. **Caching compatibility**: the data endpoint must stay cacheable (`NodeGridHandler` uses `@cache_page(WEEK_SECS)`). The client's `ajax.data` hook deletes `draw` (and any other varying params) before send — same idea as today's `deleteNdParam` (`grid.js:789`) — and the response's `draw` is restored client-side from the request. Params must serialize in a stable order.
3. **Rows-per-page persistence — built as a general DataTables feature** (users like it): `DatatableConfig` gains an optional `grid_name`; when set, the definition JSON carries `pageLength` from `UserGridConfig.get(user, grid_name).rows` plus a flag telling `datatable_definition.js` to POST `set_user_row_config` on length change (server value takes precedence over the localStorage default, which stays as the fallback for tables without a `grid_name`). The adapter sets `grid_name` from the jqGrid caption, so converted variant grids keep their existing per-user row counts (and the export `rows=0` contract stays in sync); existing DataTables can opt in one by one later.

New client pieces:

4. **`variant_grid_formats.js`** (or extend `datatables_client_renderers.js`) — port the formatters from `grid.js`/`variant_details_link_grid.html` to DataTables `(data, type, row)` signature. jqGrid formatters that read `options.colModel` metadata instead read the definition's `extra` block (closed over at table-build time):

   | jqGrid formatter | Notes |
   |---|---|
   | `detailsLink` | reads hidden row fields (`locus__contig__name`, `clinvar__*`, `internally_classified*`) — these already arrive in the row payload; `isNodeVisible` reads `extra.analysisNode.visible` |
   | `tagsFormatter`, `tagsGlobalFormatter` | keep reading `window.variantTags` (deliberate — keeps the cached grid response tag-independent) |
   | `clinvarLink`, `cosmicLink`, `formatDBSNP`, `formatPubMed`, `formatOntologyTerms`, `omimLink`, `formatMasterMindMMID3`, `formatMavedbUrnLinks`, `formatClinGenAlleleId`, `gnomadFilteredFormatter`, `unitAsPercentFormatter` | mechanical ports |
   | `geneSymbolLink`, `geneSymbolNewWindowLink` | filter-child link gated on `extra.analysisNode.visible` |
   | `formatVariantTagFirstColumn`, `formatVariantTag`, `formatAnalysis` | variant-tags page |
   | `linkFormatter` `formatter_kwargs` | becomes per-column `renderKwargs` in the definition |

5. **DataTableDefinition support** for: definition-supplied `pageLength`, `extra` metadata handed to renderers, `deferLoading` flag (needed for the analysis deferred-load behaviour in Stage 4), the `draw`-stripping ajax hook, and "~N" approximate-count pager text.

Tests: adapter unit tests (definition JSON shape, param translation incl. packed-genotype `index` round-trip, known-count passthrough with zero `COUNT(*)` queries, draw handling). The existing grid-class tests keep passing untouched.

### Stage 2 outcome — what actually landed

**Adapter** — `library/django_utils/jqgrid_datatable_adapter.py`:

- `JqGridDatatableView(View)` with class attrs `grid`, `csv_download`, `scroll_x` (default True),
  `defer_loading` (default False), `cache_stable_params` (default True). It keeps `JQGridView`'s URL shape:
  register it on the existing `<slug:op>` path and the client points at `.../handler/` while
  `op == "download"` still streams `grid_export_request` (raw values, `rows=0`).
- `?dataTableDefinition=1` returns the definition; anything else returns data. Both are encoded with
  `library.jqgrid.jqgrid.json_encode` (strict, `allow_nan=False`) rather than `JsonResponse`.
- Definition keys: `columns`, `order`, `extra`, `gridName`, `pageLength`, `lengthMenu`, `ajaxType`
  (`"GET"`), `cacheStableParams`, `deferLoading`, `approximateCount`, `downloadUrl`, `csvName`,
  `scrollX`, `searchBoxEnabled`, `downloadCsvButtonEnabled`. Per column: `data` (colmodel `name`),
  `label`, `orderable` (from `sortable`), `visible` (from `hidden`), `className`, `render`, `width`,
  `headerTitle`, `renderKwargs` (from `formatter_kwargs`).
- Every colmodel becomes a column, hidden ones included, so a client column index maps back to its
  colmodel. `datatable_columns_from_colmodels()` and `datatable_order_from_config()` are module level
  and unit tested on their own.
- Data envelope: `{recordsTotal, recordsFiltered, data, approximateRecords?, draw?}`. Both counts come
  from the one paginator count, so `get_known_count()` flows through with zero `COUNT(*)`. Rows keep
  their `.values()` field-name keys and pass through `DatabaseTableView.sanitize_value` /
  `limit_value_size`. `draw` is echoed only when the request carried one.
- `translate_params()`: `length`→`rows` (0 = no paging), `start//length + 1`→`page`,
  `order[0][column]`→`colmodels[i]["index"]`→`sidx` + `order[0][dir]`→`sord`. `filters` / `_search`
  and any page params pass through untouched, so the packed-genotype index string round-trips exactly.
- `JqGrid.get_datatable_extra()` is the grid-wide metadata hook (jqGrid formatters read the equivalent
  off `options.colModel`). `AbstractVariantGrid` returns `{"genomeBuild": ...}`; `VariantGrid` adds
  `analysisNode: {visible}` and `sortingDisabled`.
- `downloadUrl` is filled in automatically when `csv_download=True` and the URL has an `op` kwarg, and
  the client renders a toolbar CSV button from it carrying the table's live ajax params — Stage 3 pages
  get the server streaming download without page-level wiring.

**Rows per page** — `DatatableConfig.grid_name` (default None) is the general opt-in; when set,
`DatabaseTableView.json_definition()` emits `gridName`/`pageLength`/`lengthMenu` from
`UserGridConfig.get_rows_and_selections(user, grid_name)` (moved off `JqGridUserRowConfig`, which now
calls it), and the client POSTs `set_user_row_config` on length change alongside its localStorage write.
The server value wins over the localStorage default. The adapter opts in automatically for any
`JqGridUserRowConfig` grid, keyed on the caption.

**Client** — `variantgrid/static_files/default_static/js/variantgrid_formats.js`, namespace
`VariantGridFormat`, loaded from `uicore/page/base.html`:

- Renderer signature is `(data, type, row, ctx)` where `ctx = {extra, kwargs}` — `extra` is the
  definition's grid-wide block, `kwargs` the column's `renderKwargs`. `DataTableDefinition` closes
  over both at table build time.
- The formatter bodies **moved** here from `grid.js` (they are not duplicated). `grid.js` keeps a
  marked `$.fn.fmatter` shim that delegates to `VariantGridFormat.*`, mapping jqGrid's
  `(cellvalue, options, rowObject)` onto the ctx object — so the pages still on `{% jqgrid %}` run the
  same implementations. Shared page helpers (`createGridLink`, IGV, tag helpers,
  `load_variant_details`, `inAnalysis`) stayed in `grid.js` and are called from both.
- `gnomadFiltered` reads the genome build from `ctx.extra.genomeBuild`, falling back to
  `ANALYSIS_SETTINGS` — so it works on a converted page with no `ANALYSIS_SETTINGS` shim.
- The jqGrid → renderer mapping is `JQGRID_FORMATTER_TO_CLIENT_RENDERER` in the adapter module; a new
  formatter needs an entry there and a `VariantGridFormat` member.
- `snpdb/templates/jqgrid/variant_details_link_grid.html` and
  `variantopedia/templates/variantopedia/grids/gene_variants_grid.html` no longer register formatters
  (`grid.js` already does); both are now bare `{% extends "jqgrid/jqgrid.html" %}` and get deleted when
  their last page converts.

**`datatable_definition.js`** gained: definition-supplied `pageLength`/`lengthMenu`, `ajaxType`,
`deferLoading` (built as `deferLoading: 0`), `cacheStableParams`, `approximateCount`, the renderer
`ctx` argument, `headerTitle` on the `<th>`, the rows-per-page POST, and the `downloadUrl` toolbar
button. Everything is gated on a definition flag, so existing DataTables are unchanged.

**Cache stability**: with `cacheStableParams` the client's ajax `data` hook returns a minimal,
key-sorted object — `start`, `length`, `order[0][*]`, `search[value]` where searching, plus whatever
the page's own data function added — and drops `draw` and the whole `columns[...]` array (which alone
would be ~600 params on a wide cohort grid). `draw` is stashed on the definition instance and restored
onto the response by `ajax.dataSrc`. Verified: identical grid state produces a byte-identical
querystring.

**Tests**: `analysis/tests/test_jqgrid_datatable_adapter.py` (17) — column/order mapping, param
translation including the packed-genotype round-trip, envelope shape, known-count with zero `COUNT(*)`,
draw-only-when-sent, `pageLength` from `UserGridConfig`. All existing grid tests pass untouched.

**Naming**: new names spell it `variantgrid`, one word (hence `variantgrid_formats.js`).

## Stage 3 — Variantopedia + gene pages — **DONE**

Convert in this order within the one Stage 3 PR (each lands as its own commit so it can be reviewed and reverted independently):

1. **Nearby variants** (`nearby_variants.html`, 5 grids) — simplest: straight conversion of the `{% jqgrid %}` tags to `<table data-datatable-url=...>` pointing at the `handler` op of `JqGridDatatableView.as_view(grid=NearbyVariantsGrid, csv_download=True)`. The download button comes from the definition's `downloadUrl`.
2. **Gene symbol page** (`view_gene_symbol.html`) — `extra_filters` (hotspot protein position, tag clicks) move to a `data-datatable-data` function + `table.ajax.reload()`; CSV switches from the client-side `JSON2CSV` hack to the server streaming download (`csv_download=True` on the URL — raw values, unlimited rows, matches export-consistency rule in `snpdb/grids.py:486`).
3. **All Variants** (`variants.html`) — filter panel's `getExtraFilters` becomes the `data-datatable-data` function (mutating `data` in place, as `variantTagsDatatableFilter` does); the "~N" pager text and the CSV button both come through the definition already.
4. **`TaggedVariantGrid`** (`variant_tags.html`) — completes the variant tags page (`VariantTagsGrid` already converted in Stage 1); row delete via an action column posting to the existing delete endpoint; CSV keeps the dedicated `tagged_variant_export` view.

Per-column ad-hoc filtering (jqGrid `searchGrid` dialog: `search=True` on All Variants, `showFilter()` on the gene page): replaced by the filter-builder component built in Stage 4 for FilterNode. These pages ship with their page-level filters only; the column filter dialog arrives with Stage 4 (resolved decision 1).

Cleanup as pages convert: their `{% jqgrid %}` usage, `ANALYSIS_SETTINGS` shims, and grid-specific `init_func`s go, along with the `JSON2CSV` client CSV helpers. `variant_details_link_grid.html` (now a bare `{% extends %}`) is deleted once `variant_tags.html` is the last page off it in Stage 3, and `grids/gene_variants_grid.html` is already unreferenced.

Tests: switch `variantopedia/tests/test_urls.py` and `genes/tests/test_urls.py` entries from `_test_jqgrid_urls_contains_objs` to the datatable harness; keep the grid-logic tests (they exercise the classes directly).

### Stage 3 outcome — what actually landed

Four `JqGridDatatableView` URLs (`nearby_variants_grid` / `nearby_gene_variants_grid`,
`gene_symbol_variants_grid`, `all_variants_grid`, `tagged_variant_grid`) and four converted pages.
Every page now renders `<table data-datatable-url="…/handler/">`; no `{% jqgrid %}` variant grid
is left.

- **Nearby variants** — 5 plain tables, no page filter state. CSV comes from the definition's
  `downloadUrl` (`csv_download=True`).
- **Gene symbol page** — `geneVariantsDatatableFilter` sends `geneVariantsGridParams` as
  `extra_filters`; hotspot clicks, `tagClick` and "Show all" all call `reloadGeneVariantsTable()`
  (`ajax.reload()`, which resets to page 1). CSV moved off the client `JSON2CSV` hack onto the
  server streaming download. The jqGrid `searchGrid` "Filter grid..." link is gone - the column
  filter dialog arrives with Stage 4 (resolved decision 1).
- **All Variants** — `allVariantsDatatableFilter` sends the filter panel's `gridExtraFilters` as
  `extra_filters`; `filterGrid()` reloads the table. Per-column search and the client CSV are
  replaced by the page filters and the server download respectively.
- **Variant tags** — `taggedVariantDatatableFilter` sends the shared `gridParams` as a single
  `extra_filters` blob (the grid class's own contract), beside `variantTagsDatatableFilter`'s plain
  params for the tag-centric table. `reloadGrid()` now reloads two DataTables - the Stage 1 fan-out
  to `trigger("reloadGrid")` is gone. CSV stays on `tagged_variant_export`, reached from a
  `data-toolbar` button, because that view forces the genomic-order sort.

**Layout** — the converted pages are 87 columns wide with cells holding hundreds of links, which the
default DataTables sizing cannot carry:

- DataTables only reads `columns.width` when `autoWidth` is on, and `autoWidth` sizes to content -
  so the widths the adapter sent were being ignored (the 110px variant column rendered at 64px, so
  `detailsLink`'s somatic box wrapped to a second line) and a PubMed cell with 40 ids rendered
  2267px wide, pushing every later column off screen.
- The adapter now sends a width for *every* column (colmodel width, else jqGrid's own 150 default)
  plus `tableClass: "variantgrid-datatable"`. `datatable_definition.js` puts the width on the `<th>`
  itself and sets the table's width to their sum, so `table-layout: fixed` engages - it is inert
  while the table width is `auto`.
- `.variantgrid-datatable` in `global.scss` clips each cell to one line (`overflow: hidden`,
  `nowrap`, ellipsis) with jqGrid's 2px/4px cell padding, flexes `.variant_id-container` so
  `detailsLink`'s floated boxes stay on one row (ClinVar, germline, somatic, IGV), drops the
  `external-link` trailing icon inside cells, and trims DataTables' 26px sort-arrow header gutter to
  16px. Rows are 20px, the way jqGrid rendered them.

Also:

- `translate_params` only sets `rows`/`page` when the request actually carried a `length`, so a
  bookmarked grid URL pages off the grid's own `rowNum` instead of asking for everything (which
  `JqGrid.get_data` cannot serve - it dereferences a null paginator).
- `grid.js` registers its `$.fn.fmatter` shim from `$(document).ready` - `$.fn.fmatter` does not
  exist until `include_jqgrid.js` has run, which the jqGrid pages do from the body, after grid.js.
  Before this the Stage 2 shim silently no-opped (pre-Stage-2 the formatters were global functions,
  which jqGrid falls back to). The analysis node grid is the page that needed it.
- `showGridCell` (used by `detailsLink`'s ClinVar / internal-classification boxes) matches the
  adapter's `dt-<column>` class as well as jqGrid's `aria-describedby`.
- Deleted: `snpdb/templates/jqgrid/variant_details_link_grid.html` and
  `variantopedia/templates/variantopedia/grids/gene_variants_grid.html` (both already bare
  `{% extends %}`), the gene page's dead `$(window).resize` jqGrid `setGridWidth` handler and its
  `<table id="grid">`/`<div id="pager">` leftovers.
- `TaggedVariantGrid` gets no delete action column: its URL carried `delete_row=True` but
  `{% jqgrid %}` was called with the default `delete=False`, so the page never showed one - and the
  grid's model is `Variant`, so `JqGridUserRowConfig.delete_row` would have targeted the variant
  rather than the tag. Tag deletion stays on the tag-centric table.
- `JSON2CSV` / `download_grid_json_as_csv` stay: seqauto, pathtests, patients and the node-editor
  popup grids still use them.
- Tests: `variantopedia`/`genes` URL tests moved to the datatable harness (plus a nearby-variants
  entry); two new adapter tests cover the no-`length` default and the `downloadUrl` reversal.

**Browser harness** (`scratchpad/harness/`, throwaway) - the layout work above was measured rather
than guessed, and Stage 4's wide-cohort prototype wants the same rig:

- A static page that loads the real `variantgrid_formats.js` / `grid.js` / `datatable_definition.js`
  / `global.css` from `static_files/default_static`, is handed a definition JSON dumped from the
  grid class and a canned rows array, and builds a real `DataTableDefinition`. `python3 -m
  http.server` + `google-chrome --headless=new --dump-dom --screenshot --virtual-time-budget`. No
  server, no login. Add a `<pre>` and write `offsetWidth`/`getComputedStyle` readings into it - the
  numbers are what found both width bugs.
- For the real page, a session works around login: Django's session backend here is
  `cache`, so mint one with `import_module(settings.SESSION_ENGINE).SessionStore`, then serve a
  one-line page from the harness port that sets `document.cookie = "sessionid=..."` and redirects to
  `localhost:8000`. Cookies ignore the port, so the two servers share it.

### Stage 3 open item — the gene page hotspot graph

The user reports the plotly lollipop graph on the gene page has disappeared. Investigated and
**not reproduced**: on the dev box the container renders, `load_hotspot_graph()` fires, the ajax
completes, and `HotspotGraphView` answers with its own "No variants" branch -
`_pick_transcripts` returns `None` for every gene tried across both builds (a local annotation gap:
transcript version never pairs with a Pfam accession), so this DB cannot draw the graph at all.
`git diff | grep -i hotspot` touches nothing but the kept click handlers. Two open readings:

1. The graph itself is genuinely gone on their deployment - needs a gene symbol + genome build that
   used to draw one, then re-run the browser check above against it.
2. What they miss is the jqGrid `searchGrid` column-filter dialog (the "Filter grid..." link), which
   is how you would filter that grid on gnomAD AF. **Reading 2 is now addressed**: the Stage 4
   filter-builder was built first and the gene page has its "Filter grid..." button back, offering 77
   fields including gnomAD AF. If the graph itself is still missing on their deployment, reading 1
   stands and wants a gene symbol + genome build that used to draw one.

## Stage 4 — Analysis node grid — **DONE**

The big one. All of `node_data_grid.html` + `grid.js` `setupGrid` machinery moves to DataTables:

1. **Endpoints**: `NodeGridConfig` gains the definition JSON (same URL shape — it already varies by `node_id`/`node_version`/`extra_filters`, so client definition caching stays correct; `ccc_version_id` is in the postData/hash as today). `NodeGridHandler` speaks the DataTables envelope via the adapter, keeping `@cache_page`, the per-user cache lock, the shareable-node redirect, and `_NODE_GRID_ALLOWED_PARAMS` (extended with the translated param names).
2. **Grid setup JS**: replace `setupGrid`/`loadNodeGridData` with a `DataTableDefinition` build inside `node_data_grid.html`; per-node `postData` (`node_id`, `version_id`, `ccc_id`, `ccc_version_id`, `extra_filters`, `zygosity_samples_hash`, FilterNode `filters`) becomes the ajax `data` function. Deferred load: build with `deferLoading` and trigger the first fetch from `gridPaneShown` — preserving `node_grid_auto_load_max_variants` behaviour. Keep the loading overlay, `resizeGrid` (now DataTables column adjust), the component registry rendezvous, error envelope handling (`errors` / `non_fatal` → `deleteNodesFromDOM`), and the single-active-request XHR tracking.
3. **Row interactions**: selection checkboxes (`selectVariant`, re-check from `selectedVariants` on draw callback), tag add/remove (`showTagAutocomplete`, `tagClickHandler`), variant-details open, IGV links, `createFilterChild` from gene symbol cells — all move onto DataTables `drawCallback`/delegated events; the formatters were already ported in Stage 2.
4. **Sorting/pager**: `sorting_disabled` flows through the definition (`orderable: false` on every column, no initial order) — banner already exists; known counts flow through the adapter.
5. **Exports**: buttons keep hitting `node_grid_export`; `nodeGridExportInfo` builds params from the table's ajax `data` function instead of jqGrid `postData` (server contract unchanged — export drops paging/sort anyway, and `views_grid.py` already ignores the client `version_id`).
6. **FilterNode editor**: the jqGrid `searchGrid` dialog is replaced by a small standalone filter-builder component (fields + types + choices come from the definition JSON; operator list matches `JqGrid.get_q`'s ops) that emits the same `filters` JSON the server already persists and parses via `FakeFilterGrid` — existing saved FilterNodes keep working, and the same component serves the ad-hoc column filtering on the Stage 3 pages.
7. **Sample page reuse** (`sample_variants_tab.html`) and `base_related_analyses.html` follow automatically since they call the same config/handler URLs.
8. **Column summary / node editor forms / node_debug**: untouched — they consume the `VariantGrid` class server-side, which is unchanged.

Tests: keep `test_grid_export*`, `test_grid_known_count`, `test_node_grid_sort_limit`, `test_node_grid_auto_load` (class-level, unchanged); add endpoint tests for the definition/data envelope and the cache-key stability (no `draw` in cached URL).

### Stage 4 outcome — what actually landed

The node grid is a DataTable. `setupGrid`/`loadNodeGridData`/`setRowChangeCallbacks`/`deleteNdParam`
are gone from `grid.js`; nothing on the analysis page calls jqGrid any more.

**Grid CSS, measured first** (browser harness, 166 sample cohort = 1103 columns):

- `box-sizing: border-box` on `.variantgrid-datatable` cells. Without it every column came out
  **exactly 20px wider** than its colmodel (4px left + the 16px sort gutter), so the 45450px table
  rendered at 67470px and the 110px variant column at 130px. Now all 1101 visible widths land exact.
- The sort gutter is dropped on `th.sorting_disabled` - on the 25px per sample genotype columns it
  was most of the column, and there is no arrow to make room for.
- DataTables draws the row separator on every cell but the first row's, so row 1 measured 22px and
  the rest 23px. A transparent 1px border on the first row plus `line-height: 16px` makes **every
  row exactly 21px** - the predictable height the auto-fit follow-up needs.

**Endpoints** (item 1) — `NodeGridConfig` now answers `datatable_definition(grid, defer_loading=True)`
plus a `postData` block; `NodeGridHandler` answers `datatable_data(...)`. Both keep `@cache_page` +
`vary_on_cookie`, the per-user lock, the shareable-node redirect and the error/non-fatal branches.
`_NODE_GRID_ALLOWED_PARAMS` gained `start`/`length`/`order[0][column]`/`order[0][dir]` and lost
`_filters` (only `deleteNdParam` ever set it). The adapter's definition/data/translate/row helpers
moved to module level so these two views use them without a `JqGridDatatableView`.

**Two URLs, one table** — the definition is per node version (it caches with the node), the data
endpoint is one URL for the whole analysis, so `DataTableDefinition` gained `definitionUrl`
(defaulting to `url`, and keyed on for the definition cache). Every other table is unchanged.

**Grid setup** (item 2) — `setupNodeGrid(config_url, handler_url, ...)` in `grid.js` returns a promise
resolving once the table is built. `DataTableDefinition` gained the hooks it needs: `onDefinition`
(the config endpoint answers node errors instead of a table), `onData` (`non_fatal`/`deleted_nodes`
&rarr; `deleteNodesFromDOM`), `onBeforeSend` (the single-active-request XHR handle - it chains
`$.ajaxSettings.beforeSend` so the global CSRF header survives) and `onLoadError`.

- Deferred load: `deferLoading` builds the table and holds the query; `loadNodeGridData()` is
  `ajax.reload()`. `setup()` skips its `columns.adjust().draw()` on a deferred table - that draw *is*
  the fetch being held back.
- A deferred grid fires `gridComplete` once with no rows, the way the old `datatype: 'local'` init did.
  Without it `registerComponent(unique_code, GRID)` never fires and the editor's `everythingLoaded`
  sits behind its overlay until the user clicks "Show grid" - found by driving both paths in the harness.
- `nodeGridHasData()` (has the table ever had an ajax response) replaces the `datatype !== 'local'` test.

**Row interactions** (item 3) — `draw.dt` runs `gridComplete` + `gridCompleteExtra`, which is where
selection checkboxes, tag click handlers and the variant-link right-click handling get rebound. The
renderers were already ported in Stage 2, so the cells themselves needed nothing.

**Sorting/pager** (item 4) — `sorting_disabled` flows through the definition as `orderable: false` on
every column with no initial order, and known counts through the envelope (the pager reads
"Showing 1 to 10 of 4,920 entries" off a stored `NodeCount`, no `COUNT(*)`). `resizeGrid` now sets
`max-height` on `.dataTables_scrollBody` and calls `columns.adjust()`.

**Exports** (item 5) — three toolbar buttons built by `node_data_grid.html` once the table is
(CSV, Canonical transcript CSV where the analysis has a collection, VCF), plus the two placeholder
links, all still through `registerNodeGridDownloadButton`. `nodeGridExportInfo` builds its params
from `dataTable.ajax.params()`, falling back to the definition's `postData` when the grid was never
loaded - the placeholder offers CSV/VCF on a node whose rows the user never fetched. It drops
`start`/`order[0][*]` and sets `rows=0`/`length=0`; the server contract is unchanged.

**Filter builder** (item 6, built first) — `variantgrid_filter_builder.js`, `VariantGridFilterBuilder`.

- Fields, types and choices come from the definition's `filterBuilder` block
  (`datatable_filter_fields_from_colmodels`: colmodel `search`, `stype`/`searchoptions.value`,
  `sorttype`). Columns the engine can't filter are left out - `search=False`, and packed genotype
  columns, whose `index` is a sample position rather than a field path. A 1103 column cohort grid
  offers 96 fields; the gene page offers 77, gnomAD AF among them.
- Operations come from `FILTER_OPERATIONS`, lifted out of `format_operation` in `library/jqgrid/jqgrid.py`
  so `get_q`, `FilterNodeItem` and the builder all read one list.
- It emits `{"groupOp": ..., "rules": [{"field", "op", "data"}]}` - byte for byte what `FakeFilterGrid`
  parses and what `FilterNodeView.post` writes `FilterNodeItem` rows from. A saved rule on a column
  the grid no longer offers keeps its own option, so loading and re-saving a FilterNode doesn't drop it.
- `filternode_editor.html` mounts it against the node grid's definition and POSTs to the same URL;
  the jqGrid `searchGrid` dialog and its `jqGridFilterSearch`/`jqGridFilterReset` bindings are gone.
- Every adapter served grid gets a "Filter grid..." toolbar button and panel for free
  (`JqGridDatatableView.filter_builder`, default on), which is the Stage 3 pages' column filtering back.
- **The analysis node grid is the exception**: `filterBuilderToolbar: false` sends the fields and
  operations down for the FilterNode editor but shows no button on the grid itself - filtering an
  analysis is what FilterNode is for.

**Toolbar placement** - `DataTableDefinition`'s dom string gained a `.bottom-toolbar` slot inside the
existing `<"bottom">` flex row, and `setupDom` moves `[data-toolbar-bottom="#tableId"]` into it. The
node grid's CSV / Canonical transcript CSV / VCF buttons sit under the grid on the pager's line,
between "Showing 1 to 10 of 4,920 entries" and the page controls. `data-toolbar` (above the grid) is
unchanged for everything else.

**Compact controls** - default Bootstrap sizing put **60px** under the table: `.dataTables_info` and
`.dataTables_length` stacked as blocks (33px + 27px) above 32px round pager buttons - nearly three
21px rows the user couldn't see. The definition's `compactControls` puts `dt-compact` on the wrapper,
which lays `.showing` out as a flex row and pins every control to 22px. Measured **60.3px -> 22.8px**,
buying back about 1.7 rows. The adapter sets it for all its grids; native `DatatableConfig` tables
don't send it and are unchanged.

The node editor's tab strip was trimmed to match (`analysis_nodes.css`): `.bottom-pane-tab` (the
horizontal mode Editor/Grid/Summary/Doc bar) **31px -> 24px**, and the vertical mode jQuery UI strip
`#node-editor-tabs > ul.ui-tabs-nav` **42.5px -> 26px** - jQuery UI sizes its tab links in em
(`.5em/1em`), which cost 42px of pane before any content. Together with the grid controls that is
**~45px back in horizontal mode and ~54px in vertical**, so two to three more rows on screen.

**Sample page** (item 7) — `sample_variants_tab.html` calls `setupNodeGrid` with both URLs; its
`<div id='pager-...'>` is gone. `base_related_analyses.html` builds its export params itself and was
already on the datatable for its tags grid, so it needed nothing.

**Untouched** (item 8) — column summary, the node editor forms and `node_debug` consume `VariantGrid`
server-side, which is unchanged. `library/jqgrid` stays as the engine.

**Evidence**

- Live endpoints on the dev box: the config answers 111 columns / `deferLoading: true` / postData
  carrying node_id, version_id, ccc_id, ccc_version_id, extra_filters, zygosity_samples_hash; the
  handler answers `recordsTotal: 4920` with no `draw`. Three identical requests returned the same
  sha256 in ~8ms - the week-long `@cache_page` is intact.
- The client's request URL is exactly the key-sorted param set, matching the hand-built one byte for byte:
  `?ccc_id=1&ccc_version_id=0&extra_filters=default&length=10&node_id=238&start=0&version_id=3&zygosity_samples_hash=...`
- A saved FilterNode (`locus__contig__name eq MT`) round trips: its rules come back on the definition,
  the builder offers both the field and the operation, and the handler returns MT rows only.
- Tests: `analysis/tests/test_node_grid_datatable_endpoints.py` (9 - definition shape, deferred rows,
  postData, envelope, paging, column filters, and the cache key: every client param whitelisted,
  identical state gives an identical URL, `draw` never reaches it) and 5 filter builder tests in
  `test_jqgrid_datatable_adapter.py` (25 total). 1066 tests pass across analysis/snpdb/variantopedia/
  genes/library. pylint unchanged from master.

**Unrelated crash found while verifying** — `major_operation` (`library/django_utils/major_operation.py`)
released its per-user slot with a bare `cache.decr`, which Redis raises on for a missing key. The
counter's TTL (`MAJOR_OPERATION_SLOT_EXPIRE_SECONDS`, 10 min) is set when the key is created and never
extended, so an operation still running as it lapses crashed on the way out — and because the release
is in a `finally`, the 500 replaced whatever the operation had actually produced. Since a node grid
gets a 5 minute `statement_timeout`, the slow operations this limit exists to protect are exactly the
ones that tripped it. Now released through a `_release_slot()` helper that tolerates the key having
gone, with two tests that fail against the old code. Predates this project — `views_grid.py`'s
`with major_operation(...)` is untouched by the conversion.

**Browser harness** (`scratchpad/harness/`, throwaway) — rebuilt for this stage. `dump_definition.py`
writes a real definition from a `CohortNode` over the 166 sample cohort, `make_rows.py` cans matching
rows, `index.html` builds a real `DataTableDefinition` from them and writes `offsetWidth`/computed
style readings into a `<pre>`; `node_grid.html` does the same for `setupNodeGrid` against a live
`node_grid_config`/`node_grid_handler` capture, with the analysis page's globals stubbed - that is
what found the deferred-load `registerComponent` bug. `python3 -m http.server` + `google-chrome
--headless=new --dump-dom --screenshot`. The real pages are reachable too: mint a session with
`SessionStore` (the backend is `cache`), then `login.html?sid=...&next=...` sets the cookie and
redirects to `localhost:8000` - cookies ignore the port.

**Follow-up goal (post-conversion): auto-fit rows to the grid pane height.** Deferred until after the main conversion, but Stage 4 keeps the hooks open — nothing extra to build now:
- Runtime page-length changes are native DataTables (`table.page.len(n).draw()`), and the adapter already passes `length` through as jqGrid `rows`, so any value works server-side.
- The computation slots into the existing `resizeGrid` hook: `floor(paneHeight / rowHeight)` needs a predictable row height, which the Stage 4 CSS now gives it — **every analysis row is exactly 21px**, measured in the harness across a 1103 column grid.
- **Round the fitted value to the standard buckets (10/15/20/25/50/100)** — a raw fit like 23 rows makes `length` vary per user pane size and fragments the week-long `node_grid_handler` cache.
- Auto-fit becomes a mode alongside the `UserGridConfig` fixed value (e.g. `rows=0` or a flag meaning "auto"), decided when the feature is built.

## Stage 5 — Cleanup

Stages 1-4 already removed most of what this stage was written to remove — `setupGrid` and its
machinery, `variant_details_link_grid.html`, `grids/gene_variants_grid.html`, and the FilterNode
`searchGrid` plumbing are gone, and all six grids' URLs already point at the adapter. What is
actually left:

- **`$.fn.fmatter` shim in `grid.js`** (marked `LEGACY-JQGRID`) — the only reason the converted
  formatters are still reachable under their jqGrid names. Delete it once no `{% jqgrid %}` page
  renders a variant column: check `snpdb/templates/jqgrid/jqgrid.html` consumers first.
- **`FilterNode.get_extra_grid_config`'s `search = True`** (marked `LEGACY-JQGRID`) — opened jqGrid's
  dialog; the DataTables builder reads the rules off `postData` instead.
- **jQuery UI / jqGrid assets on the converted pages** — `include_jqgrid.js` is now only pulled in by
  `snpdb/templates/jqgrid/jqgrid.html` and `patients/tags/phenotype_entry_tag.html`, so the analysis
  page no longer needs it. The node editor's own tabs are still jQuery UI (`#node-editor-tabs`), so
  jQuery UI itself stays.
- **`JSON2CSV` / `download_grid_json_as_csv`** — still used by seqauto, pathtests, patients and the
  node-editor popup grids; goes when they convert.
- `library/jqgrid`, `jqgrid.html` and the jqGrid JS library stay for the ~16 unconverted grids
  (tracked separately — the remaining list is in issue #1462's family). `library/jqgrid` is also the
  adapter's engine, so it outlives the client library either way.

Start with `grep -rn LEGACY-JQGRID` — that is the whole worklist for the client side.

### Stage 5 outcome — what actually landed

Both `LEGACY-JQGRID` items are gone; `LEGACY-JQGRID-ENGINE` on `FILTER_OPERATIONS` stays.

- **`$.fn.fmatter` shim** (`grid.js`) — deleted with `jqGridFormatterContext`. Confirmed first that no
  `{% jqgrid %}` page renders a variant column: the 16 shimmed names are only set in
  `AbstractVariantGrid._get_standard_overrides` and `VariantGrid` (both converted), plus
  `GenesGrid`'s `geneSymbolLink` - and the two pages serving `genes_grid` (`genes.html`,
  `view_contig.html`) register their own `geneSymbolLink` in a `{% block formatter %}` that already
  won, being a later `$(document).ready`. `showGridCell` drops its `aria-describedby` branch; only
  the adapter's `dt-<column>` cells are left to scroll to.
- **`FilterNode.get_extra_grid_config`'s `search = True`** — removed. It only ever reached jqGrid
  through `JqGrid.get_config`, which no browser-facing variant grid goes through any more; the node
  grid's definition comes from `datatable_definition()`. The rules still travel on `postData`.
- **`analysis_downloads.js`'s `.ui-pg-div`** — the download buttons are Bootstrap `<a>` elements now,
  so the nav-button lookup and the `target`/`container` split it existed for are gone.
- **`view_pedigree.html`** — dropped its jqGrid css/js: that page has no grid at all (its related
  data/analyses includes are DataTables). Inspection only - this DB has no `Pedigree` rows.

**Pages still loading jqGrid assets, and why** - all now carry a `{# LEGACY-JQGRID #}` marker:

| Page / template | Why |
|---|---|
| `analysis_includes.html` | node editor tabs (gene coverage, phenotype, MOI) are `{% jqgrid %}` and load by ajax, after ready - `include_jqgrid.js`'s async append would race |
| `view_sample.html`, `view_vcf.html` | the skipped annotation tabs, same ajax timing |
| `view_cohort.html` | the two sample pickers build their own jqGrids |
| `jqgrid.html`, `phenotype_entry_tag.html` | `include_jqgrid.js` + `json2csv.js` for every remaining `{% jqgrid %}` page |

**Browser harness** (`scratchpad/harness/`, throwaway) - rebuilt as a CDP driver rather than
`--dump-dom`: `drive.py` starts headless Chrome with `--remote-debugging-port`, sets the minted
`sessionid` cookie through `Network.setCookie`, navigates, then evaluates a probe file and returns
its value alongside every `Log.entryAdded`/`Runtime.exceptionThrown`. `--virtual-time-budget` was no
use here - it expires before the grid ajax lands, so every grid read as empty. The session is minted
with `SessionStore` as in Stage 3, but `_auth_user_backend` must be `ModelBackend`: picking
`AUTHENTICATION_BACKENDS[0]` gets `AxesStandaloneBackend`, which has no `get_user` and 500s every page.

**Evidence** (live pages, this dev box)

- Converted: All Variants 105 columns / 81,595 records / 10 rows; gene symbol page 87 columns / 1 row;
  variant tags 106 columns / 8 rows beside the 9 row tag table; variant details genotype grid 168 rows.
- Analysis 22, FilterNode 223: grid shows "1 to 10 of 3,148", the builder loads the saved
  `locus__contig__name eq MT` rule, and there is no `searchmodfbox_` dialog anywhere on the page.
- Still jqGrid: `sequencing_stats/data` builds 4 grids and 4 rows; the genes page builds 10 rows with
  `<a href="/genes/view_gene_symbol/A1BG">A1BG</a>` from its own formatter; `$.fn.fmatter` carries 16-18
  names depending on page (the jqGrid built-ins plus each page's own), no variant renderers among them.
  Analysis 42's phenotype node editor still constructs its jqGrid, with `JSON2CSV` present and the
  jQuery UI `#node-editor-tabs` strip intact at 9 tabs.
- Console is clean on every page checked. The gene symbol page's 500 on
  `classification/groups/datatables` reproduces on stashed master - unrelated.
- Tests: 1328 across analysis/snpdb/variantopedia/genes/library/pedigree/patients/seqauto/pathtests.
  The 3 `seqauto.tests.test_extraction_link` errors fail identically on stashed master. pylint 10.00
  on `filter_node.py`.

### jqGrid deprecation markers (every stage)

As each stage touches code, tag anything found to exist *only* for jqGrid with a greppable marker comment — `# LEGACY-JQGRID` in Python, `// LEGACY-JQGRID` in JS, `{# LEGACY-JQGRID #}` in templates — with a short note of what still uses it. The eventual full jqGrid removal then becomes a `grep -r LEGACY-JQGRID` sweep instead of re-discovery. Code the adapter keeps as its engine (the `JqGrid` server classes, `get_q` filter parsing for `FakeFilterGrid`) gets tagged `# LEGACY-JQGRID-ENGINE` instead, meaning: survives until the eventual native rewrite.

### Sequencing against the other jqGrid conversions

The remaining ~16 simple jqGrids (seqauto, pathtests, upload, `GenesGrid`, `QCGeneCoverageGrid`, skipped-annotation grids, `CohortSampleListGrid`, and the DataFrameJqGrid popups) convert natively to `DatatableConfig` and are independent of this project — they share no code with it beyond the files deleted at the very end. They can proceed in parallel or afterwards in any order; the Stage 2 infra here (UserGridConfig rows-per-page, toolbar/export patterns) makes them easier, so doing this project first helps them, and full jqGrid client-JS deletion happens whenever the last page converts, regardless of order.

---

## Resolved decisions

1. **Column filter dialog** ships with Stage 4 (built for FilterNode, then reused for ad-hoc column filtering on the Stage 3 pages); the Stage 3 pages ship with page-level filters only in the interim.
2. **Rows-per-page** persists to `UserGridConfig` and is implemented as a general opt-in DataTables feature (Stage 2 item 3), available to all existing `DatatableConfig` tables.
3. **One PR per stage** — smaller review surface; within Stage 3 the four pages land as separate commits inside the stage PR.

## Risks

- **Analysis response caching**: any param that varies per request (`draw`, unordered dict serialization) silently kills the week-long cache. The adapter tests must assert byte-identical querystrings for identical grid state.
- **CSV consistency**: server download/export must keep emitting raw values (server-side formatters), never rendered HTML — the Stage 3 gene/All-Variants CSV switch is a deliberate behaviour *improvement* (previously rendered-HTML client CSV) worth calling out in the changelog.
- **DataTables with ~50–100 columns** (custom columns + per-sample genotype columns on big cohorts): needs `scrollX`, fixed widths from colmodel `width`, `autoWidth:false`. Prototype early in Stage 4 with a wide cohort analysis.
- **Tag cache coupling**: `tags`/`tags_global` rendering must keep reading `window.variantTags` so tag edits don't require busting the node grid cache.
- **`hide_non_admin` columns**: definition JSON is cached per URL client-side and the node config URL is per-user-agnostic — definition must be computed per-request (it is: `vary_on_cookie` covers the handler; confirm the config endpoint has the same cache treatment).
