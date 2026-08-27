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

## Stage 2 — Adapter infrastructure + ported client formatters

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

## Stage 3 — Variantopedia + gene pages

Convert in this order within the one Stage 3 PR (each lands as its own commit so it can be reviewed and reverted independently):

1. **Nearby variants** (`nearby_variants.html`, 5 grids) — simplest: `search=False`, straight conversion of the `{% jqgrid %}` tags to `<table data-datatable-url=...>` against `JqGridDatatableView.as_view(grid=NearbyVariantsGrid)` URLs. Add a toolbar Download button hitting the existing `op=download` streaming CSV (the `data-toolbar` mechanism already moves page elements into the DataTables toolbar).
2. **Gene symbol page** (`view_gene_symbol.html`) — `extra_filters` (hotspot protein position, tag clicks) move to a `data-datatable-data` function + `table.ajax.reload()`; CSV switches from the client-side `JSON2CSV` hack to the server streaming download (`csv_download=True` on the URL — raw values, unlimited rows, matches export-consistency rule in `snpdb/grids.py:486`).
3. **All Variants** (`variants.html`) — filter panel's `getExtraFilters` becomes the `data-datatable-data` function; approximate-count pager text via the definition passthrough; CSV via existing `op=download`.
4. **`TaggedVariantGrid`** (`variant_tags.html`) — completes the variant tags page (`VariantTagsGrid` already converted in Stage 1); row delete via an action column posting to the existing delete endpoint; CSV keeps the dedicated `tagged_variant_export` view.

Per-column ad-hoc filtering (jqGrid `searchGrid` dialog: `search=True` on All Variants, `showFilter()` on the gene page): replaced by the filter-builder component built in Stage 4 for FilterNode. These pages ship with their page-level filters only; the column filter dialog arrives with Stage 4 (resolved decision 1).

Cleanup as pages convert: their `{% jqgrid %}` usage, `ANALYSIS_SETTINGS` shims, and grid-specific `init_func`s go; `variant_details_link_grid.html` is deleted when the last Stage 3 page lands.

Tests: switch `variantopedia/tests/test_urls.py` and `genes/tests/test_urls.py` entries from `_test_jqgrid_urls_contains_objs` to the datatable harness; keep the grid-logic tests (they exercise the classes directly).

## Stage 4 — Analysis node grid

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

**Follow-up goal (post-conversion): auto-fit rows to the grid pane height.** Deferred until after the main conversion, but Stage 4 keeps the hooks open — nothing extra to build now:
- Runtime page-length changes are native DataTables (`table.page.len(n).draw()`), and the adapter already passes `length` through as jqGrid `rows`, so any value works server-side.
- The computation slots into the existing `resizeGrid` hook: `floor(paneHeight / rowHeight)` needs a predictable row height, so the Stage 4 grid CSS should keep analysis rows fixed-height (single line, ellipsis overflow — matching current behaviour).
- **Round the fitted value to the standard buckets (10/15/20/25/50/100)** — a raw fit like 23 rows makes `length` vary per user pane size and fragments the week-long `node_grid_handler` cache.
- Auto-fit becomes a mode alongside the `UserGridConfig` fixed value (e.g. `rows=0` or a flag meaning "auto"), decided when the feature is built.

## Stage 5 — Cleanup

- Remove the ported formatters and `setupGrid` machinery from `grid.js` (keep `nodeGridExportInfo`/IGV/tag helpers that survive), `variant_details_link_grid.html`, jqGrid CSS/`include_jqgrid.js` from the converted pages, `JSON2CSV` usage, and `{% jqgrid %}` search-dialog plumbing for these pages.
- `JQGridView` URL entries for the six grids point at the adapter view; `op=download` streaming stays.
- `library/jqgrid`, `jqgrid.html`, and the jqGrid JS library remain for the still-unconverted grids (tracked separately — the remaining list is in issue #1462's family).

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
