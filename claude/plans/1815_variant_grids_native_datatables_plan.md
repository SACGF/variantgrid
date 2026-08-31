# Variant grids on the native DataTables engine (#1815)

Written by Claude Fable 5 (claude-fable-5), 2026-08-31

Issue: https://github.com/SACGF/variantgrid/issues/1815

## Goal

One server side grid engine. The variant grids (`AbstractVariantGrid` and its four subclasses) become
`DatatableConfig` subclasses served by `DatabaseTableView`, the way every other table in the project
already is. `library/jqgrid/`, `library/django_utils/jqgrid_datatable_adapter.py` and
`library/django_utils/jqgrid_view.py` are deleted, and the pieces the variant grids genuinely need that
`DatatableConfig` lacks today (column filter rules, known/approximate counts, column widths, grid-wide
renderer context, rendered CSV/VCF export rows) become first-class `DatatableConfig` features that any
grid can use.

## Why

The jqGrid client is gone (#1785) but the variant grids still speak jqGrid on the server: column
definitions are jqGrid `colmodel` dicts (`name`/`index`/`sortable`/`sorttype`/`stype`/`searchoptions`/
`editable`/`editrules`/`searchrules`/`hidden`/`formatter`), the request vocabulary is `rows`/`page`/
`sidx`/`sord`/`_search`/`filters`, the config is `rowNum`/`rowList`/`sortname`/`sortorder`/`postData`,
and the response is `{page, total, rows, records}`. The adapter translates all of that to and from the
DataTables protocol on every request. That is ~1,000 lines whose only job is translation, a second
vocabulary every variant grid feature has to be written in (the in-flight private686 Representative
Variant column is expressed as colmodel overrides), and two engines that have already drifted: null
ordering (`nulls_first` on ascending vs `nulls_last`), count strategy (paginator vs `.count()` twice),
CSV semantics (server-formatted rows vs raw values), page size cap, JSON strictness (NaN handling).

## Today - the pieces and where they live

| Piece | Where |
|---|---|
| jqGrid engine: colmodels, `get_q`, sort, paginate, `iter_format_items`, field introspection | `library/jqgrid/jqgrid.py` (BSD, django-jqgrid derived) |
| `get_overrides` - colmodel dict factory for computed columns | `library/jqgrid/jqgrid_sql.py:2` |
| `JqGridUserRowConfig` - `UserGridConfig` rows per page, `hide_non_admin` | `library/jqgrid/jqgrid_user_row_config.py` |
| Protocol translation both ways, definition JSON, `JqGridDatatableView` | `library/django_utils/jqgrid_datatable_adapter.py` |
| `create_grid_from_request`, `grid_export_request` (rows=0 export), `VARIANT_GRID_LABEL_OVERRIDES` | `library/django_utils/jqgrid_view.py` |
| Native engine: `RichColumn`, `DatatableConfig`, `DatabaseTableView` | `snpdb/views/datatable_view.py:66,224,436` |
| `AbstractVariantGrid` - standard overrides, queryset assembly, `queryset_field` split | `snpdb/grids.py:492` |
| `CustomColumn` → fields + overrides, hidden `detailsLink` fields, per-allele annotations | `snpdb/grid_columns/custom_columns.py:15,81` |
| Packed genotype column annotations | `snpdb/grid_columns/grid_sample_columns.py:38` |
| Node grid: sample/source columns, packed sort (`_sort_items`), known count, sorting cutoff | `analysis/grids.py:76` |
| `ExportVariantGrid` - PK batched, contig ordered, unpaged | `analysis/grids.py:409` |
| Node column contributions (colmodel dicts) | `analysis/models/nodes/analysis_node.py:750,926,952`, `cohort_mixin.py:279-296`, `sources/cohort_node.py:312-331` |
| Node grid endpoints, allowed cache-key params | `analysis/views/views_grid.py:41,69,115` |
| CSV/VCF export iterator, transcript replacement | `analysis/grid_export.py:25,314` |
| FilterNode reusing `get_q`; editor `eq ''` → `nu` rule fix-up | `analysis/models/nodes/filters/filter_node.py:26`, `analysis/views/nodes/node_views.py:188` |
| Column summary form off colmodels | `analysis/views/views.py:655`, `analysis/views/nodes/node_view.py:49`, `analysis/forms/forms.py:299` |
| Other variant grids | `variantopedia/grids.py:79,191,355`, `genes/grids.py:122`, `variantopedia/views.py:963` |
| Client: stable params + `_search`/`filters`, node grid `postData`, FilterNode editor | `datatable_definition.js:137-176`, `grid.js:30-45,413-424`, `filternode_editor.html:24` |
| Client renderers (`renderContext.extra`, `renderKwargs`) | `variantgrid_formats.js` |
| Management command reading `sorttype` | `analysis/management/commands/analysis_grid_find_float_nan_data.py:47` |

The definition JSON the client consumes is already engine-neutral (`columns[].data/label/render/
orderable/visible/className/width/headerTitle`, `order`, `extra`, `filterBuilder`, `gridName`,
`pageLength`, `lengthMenu`, `downloadUrl`, `deferLoading`, `cacheStableParams`, `approximateCount`,
`compactControls`, `tableClass`). `datatable_definition.js` and `variantgrid_formats.js` keep working
as they are; the client change is limited to request param names.

## Target design

### Concept mapping

| jqGrid concept | Native equivalent |
|---|---|
| colmodel `name` / `index` | `RichColumn.name` / `key` (`sort_keys` where they differ) |
| `label`, `hidden`, `sortable`, `width`, `classes`, `headerTitle` | `label`, `visible`, `orderable`, **new** `width`, `css_class`, **new** `header_title` |
| `formatter` (+ `JQGRID_FORMATTER_TO_CLIENT_RENDERER`) | `client_renderer='VariantGridFormat.xxx'` directly |
| `server_side_formatter(row, field)` | `renderer(CellData)` with **new** `csv_rendered=True` |
| `model_field=True` (introspected) | built by a field-path → `RichColumn` factory (label, filter type, choices from the Django field) |
| `model_field=False, queryset_field=True` (annotation) | `RichColumn(key=alias)` with no filter |
| `queryset_field=False` (display only / unpacked) | `RichColumn(key=None, name=..., renderer=..., extra_columns=[packed alias])` |
| `search`/`stype`/`searchoptions`/`sorttype` | **new** `RichColumn.filter: Optional[FilterField]` (type + choices) |
| `hide_non_admin` | unused - dropped |
| `editable`, `editrules`, `searchrules`, `fixed`, `shrinkToFit` | dead - dropped |
| `formatter_kwargs` / `renderKwargs` | unused server side - dropped from the definition and `datatable_definition.js` |
| `extra_config` `rowNum`/`rowList` | `DatatableConfig.grid_name` (already reads `UserGridConfig`) |
| `extra_config` `sortname`/`sortorder` | `RichColumn.default_sort` |
| `extra_config` `postData` | definition `postData` stays, emitted by the node grid config (already a plain key) |
| `get_datatable_extra()` | **new** `DatatableConfig.get_extra()` |
| `get_known_count()` / `KnownCountPaginator` | **new** `DatatableConfig.known_count(qs) -> Optional[int]` |
| `approximate_records` | **new** `DatatableConfig.approximate_count(qs) -> Optional[str]` |
| `get_q(filters)` / `FILTER_OPERATIONS` | **new** `library/django_utils/filter_rules.py` |
| `get_items(rows=0)` for export | **new** `DatatableConfig.iter_export_rows(qs)` |
| `JqGridDatatableView(grid=..., csv_download=True)` | `DatabaseTableView(column_class=...)` + `server_csv_download=True` |
| `<slug:op>` URL segment (`handler`/`download`) | one URL; `?dataTableCsv=1` is the download (`DATATABLE_CSV_PARAM`) |
| `rows`/`page`/`sidx`/`sord` | `length`/`start`/`order[0][column]`/`order[0][dir]` (the client already sends these) |
| `_search=true` + `filters` | `filters` (JSON rules) alone |

### Engine additions (`snpdb/views/datatable_view.py`)

`RichColumn` gains:
- `width: Optional[int]` → definition `width` (px). The variant grids' `table-layout: fixed` needs every
  column to carry one; `DatatableConfig.default_column_width: Optional[int]` (variant grids set 150)
  fills in where a column has none.
- `header_title: Optional[str]` → definition `headerTitle`.
- `filter: Optional[FilterField]` where `FilterField(type: 'text'|'int'|'float'|'date'|'select',
  choices: Optional[dict])`. Columns carrying one appear in the definition's `filterBuilder.fields`.
- `null_order: NullOrder` - `LAST` (today's `DatatableConfig` behaviour, the default) or
  `FIRST_ON_ASC` (today's variant grid behaviour: nulls first ascending, last descending). Used by
  `sort_string`.
- `csv_rendered: bool = False` - when set, the server CSV writes `render_cell()` output for this column
  instead of the raw key value. Variant grids set it wherever they have a renderer (AF as percent,
  choices expanded, packed genotype unpacked), which is what makes the CSV match the grid.

`DatatableConfig` gains, each with a neutral default so existing grids are unaffected:
- `get_extra() -> JsonObjType` → definition `extra` (renderer context: genome build, node visibility,
  sorting-disabled flag).
- `filter_fields() -> list[JsonObjType]` (from columns with `filter`) and
  `apply_filter_rules(qs, rules) -> QuerySet` using `filter_rules.rules_to_q`. `DatabaseTableView`
  reads `filters` off the request (JSON) and applies it after `filter_queryset`. The definition carries
  `filterBuilder: {fields, operations}` when `filter_builder = True`, plus `filter_builder_toolbar`.
- `known_count(qs) -> Optional[int]` - returned as both record counts with no `.count()`.
  `approximate_count(qs) -> Optional[str]` - the `approximateRecords` string, and the definition gets
  `approximateCount: True` when the config sets `approximate_count_enabled`.
- `filter_rules_supplied` property so `known_count` implementations can decline when column filters narrow
  the rows (both `VariantGrid` and `AllVariantsGrid` do this today).
- Definition flags: `table_class`, `compact_controls`, `defer_loading`, `cache_stable_params`,
  `ajax_type` (`GET` keeps `@cache_page` working on the node handler).
- `max_page_length: int = 100` replacing the view's hardcoded `max_display_length`.
- `render_rows(rows: Iterable[dict]) -> Iterator[dict]` - the per-row `{name: render_cell()}` step
  factored out of `prepare_results`, so export and the transcript-replacement iterator use the same code
  as the grid page.
- `iter_export_rows(qs) -> Iterator[dict]` - `render_rows(qs.values(*value_columns).iterator())`.
  `download_csv` uses it; `csv_columns()` returns `{name, label}` for every enabled column with
  `include_in_csv` (rendered or raw per `csv_rendered`). `VARIANT_GRID_LABEL_OVERRIDES` becomes a
  `csv_label_overrides` attribute.

`DatabaseTableView`:
- `config_for_request(request, **kwargs)` receives the URL kwargs so a grid needing them (node, genome
  build) can be built from them; `DatatableConfig.get_query_param` already falls back to
  `resolver_match.kwargs`.
- JSON goes out strict (`allow_nan=False`) and `sanitize_value` maps float NaN → `None`, so a NaN in an
  annotation column renders as blank instead of breaking `JSON.parse` (today only the adapter's
  `json_encode` guards this; `analysis_grid_find_float_nan_data` exists because the data occurs).
- `get_context_data` honours `known_count` / `approximate_count`, and `count_unfiltered=False` grids
  count once.

### `library/django_utils/filter_rules.py` (new)

`FILTER_OPERATIONS`, `format_operation`, `rules_to_q(filters: dict) -> Optional[Q]` (the `get_q` body,
including the `__isnull` exclude-avoidance that keeps FilterNode off a full table scan),
`parse_filters(param: Optional[str]) -> Optional[dict]`. `resolve_field_path(opts, path)` (the
`lookup_foreign_key_field` walk through FK/one-to-one, stopping at a JSONField) moves to
`library/django_utils/__init__.py`. FilterNode and the FilterNode editor view call these directly;
`FakeFilterGrid` goes.

### Variant grid column building (`snpdb/grid_columns/`)

`custom_columns.py` returns `list[RichColumn]` (plus the sample-columns insertion position) instead of
`(fields, override_dicts, position)`:
- `variant_column_rich_column(path, **overrides) -> RichColumn` - the `field_to_colmodel` replacement:
  resolves the field via `resolve_field_path`, sets `label` from `verbose_name`, `filter` from the field
  class (`IntegerField`/`AutoField` → int, `FloatField` → float, `DateTimeField` → date, `BooleanField`
  and `choices` → select), leaves `filter=None` for a `ForeignKey` endpoint, then applies `overrides`.
  JSON key paths: verify no `VariantGridColumn.variant_column` traverses a `JSONField`
  (`add_annotations_for_json_fields`); if one does, keep the `KeyTransform` annotation in
  `AbstractVariantGrid.get_initial_queryset`, otherwise it goes.
- `VARIANT_COLUMN_OVERRIDES: dict[str, dict]` - `_get_standard_overrides` in `RichColumn` kwargs
  vocabulary (`{'clinvar__clinvar_variation_id': {'width': 60, 'client_renderer': 'VariantGridFormat.clinvarLink'}, ...}`),
  the AF percent renderers as `renderer=` + `csv_rendered=True`.
- `ID_FORMATTER_REQUIRED_FIELDS` become `visible=False` columns; `ID_FORMATTER_REQUIRED_ANNOTATIONS`
  and `tags_global` are `key=<alias>` columns whose values come from `get_variantgrid_extra_annotate`.
- `tags` (analysis only): `RichColumn(key=None, name='tags', client_renderer='VariantGridFormat.tags', orderable=False)`.

Node column contributions (`_get_node_extra_columns` / `_get_node_extra_colmodel_overrides`) merge into
one `_get_node_extra_columns() -> list[RichColumn]` on `AnalysisNode`, `CohortMixin` and `CohortNode`
(inheritance through `get_extra_columns` stays, de-duplicated by `RichColumn.__eq__` on name).

Genotype columns (`get_grid_genotype_columns_and_overrides`) become `RichColumn(key=None,
name=f'sample_{pk}_{column}', label=..., width=..., renderer=<unpack>, csv_rendered=True,
extra_columns=[packed alias], sort_keys=[sort alias], orderable=True)`. The `':'`-packed `index` goes:
`VariantGrid.ordering()` annotates the sort alias (array index transform or `Substr`) for whichever
sample column the request orders on, then delegates to `super().ordering()`. `vcf_source` is the same
shape with `orderable=False`.

### The grids

`AbstractVariantGrid(DatatableConfig[Variant])` in `snpdb/grids.py`:
- `__init__(request, ...)` sets `genome_build`, `annotation_version`, `rich_columns`, `grid_name`
  (today's `caption`), `count_unfiltered = False`, `default_column_width = 150`,
  `table_class = 'variantgrid-datatable'`, `compact_controls = True`, `ajax_type = 'GET'`,
  `cache_stable_params = True`, `scroll_x = True`, `server_csv_download = True`,
  `csv_label_overrides = {'id': 'variant_id'}`, `filter_builder = True`.
- `get_initial_queryset()` = `_get_base_queryset()` + the variantallele build filter + global zygosity
  counts + `_get_q()` + `_get_grid_only_annotation_kwargs()` (today's `get_queryset`).
- `get_extra()` → `{"genomeBuild": ...}`.
- `AllVariantsGrid`: `ordering()` returns `DEFAULT_ORDER_BY`; every column `orderable=False`;
  `known_count` + `approximate_count` off the `EXPLAIN` estimate; `extra_filters` via `get_query_json`.
- `NearbyVariantsGrid`, `GeneSymbolVariantsGrid`, `TaggedVariantGrid`: URL kwargs via
  `get_query_param`; `TaggedVariantGrid` keeps `tag_count` as an annotation column
  (`filter=FilterField('int')`); `tagged_variant_export` becomes the grid's own `?dataTableCsv=1`
  download with `locus__position` as `default_sort`, and `variantopedia/views.py:963` goes.
- `VariantGrid(request, node, extra_filters, af_show_in_percent=None)` in `analysis/grids.py`:
  `known_count` as today (node count, `count_is_deterministic`, `ANALYSIS_NODE_STORE_ID_SIZE_MAX`,
  declining when `filter_rules_supplied`); `sorting_disabled()` sets `orderable=False` on every column
  and makes `ordering()` return `-pk` only; `default_sort` from `analysis.default_sort_by_column` when
  sorting is allowed; `get_extra()` adds `analysisNode` and `sortingDisabled`; `post_data` property
  replaces `_set_post_data` / `extra_config['postData']` (FilterNode's `get_extra_grid_config` returns
  `{'filters': ...}` for it); `column(name) -> RichColumn` replaces `get_column_colmodel`.
- `ExportVariantGrid.iter_export_rows()` does the PK-batched, contig-ordered walk and
  `render_rows` per batch. `af_show_in_percent=False` for VCF as today.

### Views and URLs

- `genes/urls.py:77`, `variantopedia/urls.py:85-95`: `DatabaseTableView.as_view(column_class=...)`
  under one path (the `op` segment goes); `major_operation_name` stays a view kwarg. The templates
  (`variants.html:262`, `nearby_variants.html`, `view_gene_symbol.html`, `variant_tags.html`) reference the
  single URL.
- `analysis/views/views_grid.py`: `NodeGridHandler._get_data` builds a `DatabaseTableView`-equivalent
  response through a shared helper (`datatable_response(config, request)` factored out of
  `DatabaseTableView.get_context_data` so the node views, which sit on `NodeJSONViewMixin` for their
  error envelope and locking, reuse it). `NodeGridConfig` returns `config.json_definition()` +
  `postData`. `_NODE_GRID_ALLOWED_PARAMS` becomes `{ccc_id, ccc_version_id, extra_filters, filters,
  length, node_id, order[0][column], order[0][dir], start, version_id, zygosity_samples_hash}`.
- `node_grid_export` / `analysis/grid_export.py`: `grid_params` drop `rows`/`sidx`/`sord`;
  `node_grid_get_export_iterator` uses `grid.iter_export_rows()` and `grid.csv_columns()`;
  `_replace_transcripts_iterator` uses `grid.render_rows()`; `grid_export_csv(columns, rows,
  label_overrides)` keeps its `{name, label}` column shape.
- Column summary (`views.py:655`, `node_view.py:49`, `forms.py:299`) and
  `analysis_grid_find_float_nan_data` read `RichColumn`s (`filter.type == 'float'` for the latter).
- `FakeRequest` in `library/django_utils/__init__.py:367` stays for the Celery export, docstring updated.

### Client

- `datatable_definition.js:158`: send `filters` alone (the `_search` flag goes); `renderKwargs` /
  `col.renderKwargs` removed from `renderContext` (`variantgrid_formats.js` header comment and the two
  `ctx.kwargs` readers updated - `VariantGridFormat.link` reads its url name from the column instead,
  or goes if unused).
- `grid.js:413-424` and `filternode_editor.html:24`: `postData.filters` without `_search`;
  `nodeGridExportInfo` drops `rows`.

## Behaviour to preserve (the regression checklist)

1. Node grid: stored node count → no `COUNT(*)`; extra-filter count; column filters recount;
   non-deterministic small node recounts (`test_grid_known_count`).
2. Sorting cutoff (`ANALYSIS_GRID_SORT_MAX_ROWS`): every column unsortable, requested sort ignored,
   no initial sort, banner (`test_node_grid_sort_limit`).
3. Packed genotype sort: array columns via index transform, string columns via `Substr`; the sort
   survives the request → column index → sort alias round trip.
4. Ascending sort puts nulls first, descending puts them last, secondary sort on the default column,
   `-pk` tiebreak.
5. `AllVariantsGrid`: genomic order regardless of request, planner estimate ≥ 1M shown as "~N",
   non-selective filters match nothing (`test_all_variants_grid`).
6. Filter rules: identical `Q` for a grid column filter and a FilterNode; `isnull` avoids `~Q`;
   filter builder offers the same operations and only filterable fields (model fields, no packed
   columns, no FK endpoints); FilterNode editor loads saved rules from `postData.filters`.
7. CSV export rows equal the grid page cell values (AF in percent, choices expanded, dates local,
   genotype unpacked, `.` for missing); header `id` → `variant_id`; SYLK quoting; VCF header contig
   order and record order; canonical transcript replacement; tags summarised
   (`test_grid_export`, `test_grid_export_vcf`, `test_tagged_variant_grid`).
8. Node export cache: identical grid state → identical `params_hash`; stale `version_id` ignored.
9. `NodeGridHandler` cache key: the client's stable params are exactly the allowed set; `draw` never
   reaches the URL (`test_node_grid_datatable_endpoints`).
10. Definition: every column carries a width; hidden `detailsLink` fields present but invisible;
    `extra.genomeBuild` / `extra.analysisNode`; `pageLength`/`lengthMenu` from `UserGridConfig`;
    `downloadUrl`; `deferLoading` + no toolbar filter button on the node grid.
11. NaN in a float annotation column renders blank.

## Stages

Each stage leaves the test suite green and the grids working.

### Stage 0 - engine additions (no grid moves)
`snpdb/views/datatable_view.py`, new `library/django_utils/filter_rules.py`,
`library/django_utils/__init__.py` (`resolve_field_path`). FilterNode and `node_views.py` switch to
`filter_rules` / `resolve_field_path`. The adapter imports `FILTER_OPERATIONS` from `filter_rules` so
`library/jqgrid/jqgrid.py` shrinks rather than duplicating.

### Stage 1 - column factory and `AbstractVariantGrid`
`snpdb/grid_columns/custom_columns.py`, `snpdb/grids.py`, then `variantopedia/grids.py`,
`genes/grids.py`, their `urls.py` and templates. Port `test_all_variants_grid` and
`test_tagged_variant_grid`. `variantopedia/views.py:963` goes.

### Stage 2 - the node grid
`analysis/grids.py`, the node column contributions, `views_grid.py`, `views.py:602-660`,
`node_view.py`, `forms.py`, `grid_export.py`, the export tasks, the management command.
`grid.js`, `datatable_definition.js`, `filternode_editor.html`. Port `test_grid_known_count`,
`test_node_grid_sort_limit`, `test_node_grid_datatable_endpoints`, `test_grid_export*`,
`test_jqgrid_datatable_adapter` → `analysis/tests/test_variant_grid_datatable.py` (keep the packed
sort round trip, filter builder field/operation tests, known count, download URL, definition shape).

### Stage 3 - delete
`library/jqgrid/`, `jqgrid_datatable_adapter.py`, `jqgrid_view.py`, `library/tests/test_jqgrid_paginator.py`.
Docs: `CLAUDE.md` "Grid/table views", `library/__library_readme.md:19`,
`claude/research/library.md:126`, `claude/plans/circular_references_fixes_plan.md:103`. Bump
`CACHE_VERSION` (node grid definitions and handler responses are cached for a week under the old param
names and envelope).

## Sequencing with private686

`claude/plans/private686_representative_variant_column_plan.md` adds colmodel overrides, a sort branch
and JS on the current engine. Land it first - it is small and self-contained - and carry it across in
Stage 1 (its overrides become `VARIANT_COLUMN_OVERRIDES` entries, its sort branch a `sort_keys`).

## Verification

- `python3 manage.py test --keepdb analysis snpdb variantopedia genes library`
- Run the app: analysis node grid (sort a sample column, column filter, FilterNode editor round trip,
  CSV + canonical CSV + VCF downloads, big node placeholder → deferred load), all variants page
  (estimate shown as ~N, CSV), nearby variants, gene symbol variants with hotspot/tag extra filters,
  tagged variants + CSV, sample variants tab.
- Diff a node CSV and VCF exported before and after the change for the same node (test fixture
  analysis) - byte-identical apart from the date in the filename.
