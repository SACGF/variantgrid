# Retire `library/jqgrid` — the native rewrite (issue #1785 §3)

Follow-on from `variantgrid_to_datatables_plan.md`, whose stages 1-5 all landed (PR #1787). Every
page now renders DataTables and the jqGrid *client* is deleted. What is left is the jqGrid *server*
code, kept alive by three consumers and marked `LEGACY-JQGRID-ENGINE`:

| Consumer | Where | What it uses |
|---|---|---|
| `JqGridDatatableView` + the `AbstractVariantGrid` family | `library/django_utils/jqgrid_datatable_adapter.py`, `snpdb/grids.py`, `analysis/grids.py`, `variantopedia/grids.py`, `genes/grids.py` | the whole `JqGrid` pipeline (colmodels, queryset, filter/sort/paginate/format) |
| `ExportVariantGrid` | `analysis/grids.py` + `analysis/grid_export.py` | `get_items()`, `get_colmodels()`, `iter_format_items()`, `get_field_names()` |
| `FakeFilterGrid` | `analysis/models/nodes/filters/filter_node.py` | `JqGrid.get_q` + `FILTER_OPERATIONS` |

Plus three strays: `JqGrid.lookup_foreign_key_field` (`analysis/views/nodes/node_views.py`),
`KnownCountPaginator` (`library/tests/test_jqgrid_paginator.py`), and `jqgrid_sql.get_overrides`
(`snpdb/grid_columns/custom_columns.py`, `analysis/grids.py`).

## The shape of the rewrite

Today a converted variant page runs **DataTables client → jqGrid protocol → jqGrid engine**. The
adapter translates `start`/`length`/`order[0][column]` into `page`/`rows`/`sidx`/`sord`, calls
`JqGrid.get_data()`, and translates the `{page, total, rows, records}` envelope back. Column metadata
makes the same round trip: Django field → colmodel dict (`formatter`, `sortable`, `hidden`,
`classes`, `sorttype`, `stype`, `searchoptions`, `editrules`, `index`) → DataTables column dict.

The native version deletes the middle layer. **The variant grids become `DatatableConfig` subclasses**
like every other table in the app, and `DatatableConfig` / `RichColumn` / `DatabaseTableView` gain the
hooks the original plan listed as missing. That is the one-line summary of §3: *one table system, not
two.*

The **client contract does not change**. The definition JSON and the data envelope the JS already
consumes (`datatable_definition.js`, `variantgrid_formats.js`, `variantgrid_filter_builder.js`,
`grid.js`) stay key-for-key identical, apart from two deliberate drops noted in Stage E. That is what
makes a server-side refactor of the widest, most performance-critical grid in the app safe: the
browser cannot tell.

### What must survive (from the original plan's "Why an adapter" section)

| Behaviour | Today | Native home |
|---|---|---|
| stored `NodeCount` skips `COUNT(*)` (#1279, #1722) | `VariantGrid.get_known_count` + `KnownCountPaginator` | `DatatableConfig.get_known_count(qs)`, honoured by `DatabaseTableView` |
| `EXPLAIN`-based approximate count on All Variants | `AllVariantsGrid.get_known_count` + `get_data` | same hook + `config.approximate_records` → `approximateRecords` |
| `sorting_disabled()` + packed-genotype sort decoding (#1651) | `VariantGrid._sort_items` decodes `'cohortgenotype_134:1:samples_zygosity'` out of the colmodel `index` | `RichColumn.sort_annotation` (a real Django expression on the column) + `VariantGrid.ordering()` |
| forced genomic-order pagination on All Variants (#1663) | `AllVariantsGrid._sort_items` | `AllVariantsGrid.ordering()` |
| `variantallele` join restricted to the grid's build (#1626) | `AbstractVariantGrid.get_queryset` | `AbstractVariantGrid.get_initial_queryset` |
| server-side formatters shared verbatim with CSV/VCF export | `iter_format_items` / `get_field_formatters` | `RichColumn.server_side_formatter` + `format_rows(columns, items)`, used by both the grid and the export |
| the filter-rule vocabulary shared by `get_q`, `FilterNodeItem` and the client builder | `FILTER_OPERATIONS` in `library/jqgrid/jqgrid.py` | `library/django_utils/filter_rules.py` |

---

## Stage A — give `DatatableConfig` the hooks (`snpdb/views/datatable_view.py`)

Purely additive; every new option defaults to today's behaviour, so existing tables are untouched.

1. **`RichColumn`** gains:
   - `width: Optional[int]` → definition `"width": "150px"` (these tables lay out `table-layout: fixed`)
   - `header_title: Optional[str]` → definition `"headerTitle"`
   - `server_side_formatter: Optional[Callable[[dict, str], Any]]` — signature `(row, field)`, the
     same one the export path uses, so grid and CSV can never disagree
   - `data_type: Optional[str]` (`'int'`/`'float'`/`'date'`) — drives the filter builder's widget and
     `node_column_summary`'s quantitative/categorical split (was colmodel `sorttype`)
   - `filterable: bool` / `filter_choices: Optional[dict]` — what the filter builder offers
   - `queryset_field: bool = True` — False for columns that exist in the column list but are not
     selected from the queryset
   - `sort_annotation: Optional[Any]` — a Django expression; when ordering by this column the config
     annotates it under an alias and orders on that. This is what replaces the `index` string packing
     for genotype columns, and is a genuinely useful general addition.
2. **`DatatableConfig`** gains:
   - `get_known_count(self, qs) -> Optional[int]` (default `None` = count the queryset)
   - `approximate_records: Optional[str]` — set during counting; the view renders it as `~1.2M`
   - `get_extra(self) -> dict` (lifted from `DataFrameDatatableConfig`, which already emits `extra`)
   - `max_page_length: Optional[int] = 100` — the current hard cap becomes per-config; variant grids
     set `None`
   - `filter_builder` / `filter_builder_toolbar` / `table_class` / `defer_loading` /
     `cache_stable_params` / `compact_controls` / `ajax_type` flags, all defaulting to today's values
   - `filter_rules` (cached) — the `{groupOp, rules: [...]}` JSON off the request, and
     `filter_queryset` applies it via `filter_rules_to_q` when the config opts in
   - `row_values(self, qs) -> Iterable[dict]` (default `qs.values(*self.value_columns())`) — the one
     place row dicts are produced, so a config can run its server-side formatters over them
   - `csv_rows(self, qs)` (default raw `.values().iterator()`) — same idea for the CSV download
3. **Extract `datatable_definition(config)` and `datatable_data(config)`** as module-level functions
   over a config, with `DatabaseTableView` calling them. `NodeGridHandler` / `NodeGridConfig` are not
   `DatabaseTableView`s (they carry the per-user lock, `@cache_page` and the redirect logic), so they
   need to build both from a config directly — exactly as they call the adapter's functions today.
4. **`DatabaseTableView`** uses `config.get_known_count()` for *both* `recordsTotal` and
   `recordsFiltered` when it returns a value, honours `max_page_length`, emits the new definition
   keys, and echoes `draw` only when the client sent one (so a `@cache_page`d response is not
   discarded as stale).

Point 4 also fixes the item the original plan carried out of Stage 1: `VariantTag` (~412k rows at SA
Path) currently pays for two counts per draw. Any table can now supply one.

**Tests:** the existing datatable tests are the regression net. Add unit coverage for
`get_known_count` short-circuiting the `COUNT(*)`, `max_page_length = None`, and `sort_annotation`
ordering.

## Stage B — the filter-rule vocabulary comes out (`library/django_utils/filter_rules.py`)

The `{groupOp, rules: [{field, op, data}]}` JSON is not a jqGrid artefact — it is FilterNode's
**persistence format** (`FilterNodeItem` rows) and the wire format the DataTables filter builder
already emits. It just lives in the wrong file.

- `FILTER_OPERATIONS`, `FILTER_OPERATION_LABELS`, `format_operation` move over verbatim.
- `JqGrid.get_q` becomes `filter_rules_to_q(rules: dict, *, ignore_case: bool = True) -> Optional[Q]`,
  keeping the two hard-won details in place: `__in` splits on commas, and `isnull` passes the
  negation as the lookup value rather than inverting the `Q` (inverting generated a full table scan
  over 100M+ variant rows).
- `filter_rules_from_request(querydict) -> Optional[dict]` reads the `filters` param.
- `FilterNode._get_node_q` becomes `return filter_rules_to_q(self.get_filters())` — **`FakeFilterGrid`
  is deleted**, and with it the only reason `FilterNode` imported a grid class at all.
- `JqGrid.lookup_foreign_key_field` → `lookup_field_path()` in `library/django_utils/__init__.py`
  (`analysis/views/nodes/node_views.py` and the column builder use it).
- `KnownCountPaginator` → `library/django_utils/paginators.py`; `library/tests/test_jqgrid_paginator.py`
  → `library/tests/test_paginators.py`.

## Stage C — variant columns as `RichColumn`s (`snpdb/grid_columns/variant_columns.py`)

The column set is built per request from the user's `CustomColumnsCollection`, so it is authored as a
`{field_name: {**overrides}}` dict merged with `update_dict_of_dict_values`. That authoring format
stays (it is DB-driven and merge-shaped); what changes is that it now describes a `RichColumn`.

- `column_overrides(...)` (was `jqgrid_sql.get_overrides`) moves here.
- `build_variant_columns(model, fields, overrides, ...) -> list[RichColumn]` replaces
  `JqGrid.get_colmodels` / `field_to_colmodel` / `add_search_options` / `get_field_choices` /
  `get_datetime_fields` / `get_field_formatters`. Django field introspection still supplies label,
  `data_type`, choices (as both a display formatter and the filter builder's `filter_choices`),
  and the datetime local-time formatter.
- `format_rows(columns, items)` replaces `iter_format_items`, keeping its optimisation of probing the
  formatters against the first row once rather than catching a `KeyError` per absent field per row.
- **Override keys get native names** — this is the substance of "native", and confines the jqGrid
  vocabulary to nothing at all:

  | was | becomes |
  |---|---|
  | `formatter` | `client_renderer` |
  | `sortable` | `orderable` |
  | `hidden: True` | `visible: False` |
  | `classes` | `css_class` |
  | `headerTitle` | `header_title` |
  | `search` | `filterable` |
  | `sorttype` | `data_type` |
  | `index` (packed genotype) | `sort_annotation` |
  | `editable`, `editrules`, `searchrules`, `stype`, `fixed`, `hide_non_admin`, `formatter_kwargs` | **deleted** — nothing reads them (`hide_non_admin` has no setters; `formatter_kwargs`/`linkFormatter` has no authors) |

  `model_field` / `queryset_field` keep their names — they are `VariantGridColumn` model fields.
  `index` for ordinary fields was always just the field name, so it goes too.
- Authors updated: `snpdb/grid_columns/custom_columns.py`, `AbstractVariantGrid._get_standard_overrides`,
  `VariantGrid.get_grid_genotype_columns_and_overrides` / `get_source_columns_and_overrides`,
  `TaggedVariantGrid.TAG_COUNT_OVERRIDE`, `analysis/models/nodes/cohort_mixin.py`,
  `analysis/models/nodes/sources/cohort_node.py`, `FilterNode._get_inherited_colmodel_overrides`.
  The `AnalysisNode` hooks keep their names (`get_extra_colmodel_overrides` →
  `get_extra_column_overrides`, docstrings updated off "JQGrid").

## Stage D — `AbstractVariantGrid(DatatableConfig[Variant])`

`snpdb/grids.py`, then the five subclasses. Class names stay (`VariantGrid`, `AllVariantsGrid`,
`NearbyVariantsGrid`, `TaggedVariantGrid`, `GeneSymbolVariantsGrid`, `ExportVariantGrid`) — "grid" is
the domain word all through the UI, the URLs, `UserGridConfig` and the issue tracker, and renaming
six classes would bury the behavioural diff under churn.

- `__init__(self, request, ...)`: `DatatableConfig` gives `self.user`. The four page grids read their
  own params (`genome_build_name`, `variant_id`, `region_type`, `gene_symbol`, `extra_filters`) via
  `get_query_param` — which already reads URL kwargs off `resolver_match` — so they can be served by
  a plain `DatabaseTableView.as_view(column_class=...)`. `VariantGrid` keeps an explicit `node`
  argument: every caller already resolved it (with permission checks).
- `rich_columns` built in `__init__` from `build_variant_columns`.
- `grid_name = caption`, so `UserGridConfig` rows-per-page keeps working through the existing
  `DatatableConfig.grid_name` mechanism (this is what `JqGridUserRowConfig` did).
- `get_initial_queryset()` = today's `get_queryset()`: the build-restricted `variantallele` join
  (#1626), `annotate_global_germline_counts`, `_get_q()`, grid-only annotations, `.values(...)`.
- `filter_queryset()` applies the filter-builder rules once, via `filter_rules_to_q`.
- `row_values()` / `csv_rows()` run `format_rows`.
- `get_extra()` = today's `get_datatable_extra()`.
- Per subclass: `VariantGrid.ordering()` carries `sorting_disabled()` and the `sort_annotation` path;
  `AllVariantsGrid.ordering()` forces `DEFAULT_ORDER_BY`; `get_known_count()` moves across unchanged
  on both; `TaggedVariantGrid._get_q` keeps its `UserGridConfig.show_group_data` logic.
- `JqGridUserRowConfig.delete_row` is dropped — only `JQGridView` called it.

## Stage E — views, URLs and templates

- The four page grids move to `DatabaseTableView.as_view(column_class=...)` with
  `server_csv_download = True`. That **drops the `<slug:op>` URL segment**: the CSV comes off the same
  URL with `?dataTableCsv=1`, the mechanism every other server-side CSV table already uses.
  6 templates (`variants.html`, `nearby_variants.html` ×5 tables, `variant_tags.html`,
  `view_gene_symbol.html`), `variantopedia/urls.py`, `genes/urls.py`, and 3 test files.
- **Drop `_search=true`** — the presence of `filters` is the signal. One line of JS
  (`datatable_definition.js:160`), the `_NODE_GRID_ALLOWED_PARAMS` set, and the reader.
  Node export params-hashes change once, so in-flight `CachedGeneratedFile`s regenerate; harmless.
- `NodeGridHandler._get_data` / `NodeGridConfig._get_data` call the Stage A
  `datatable_data(config)` / `datatable_definition(config)`.
- **Deleted:** `library/django_utils/jqgrid_datatable_adapter.py`,
  `library/django_utils/jqgrid_view.py` (`JQGridViewOp`, `create_grid_from_request`,
  `export_grid_as_csv`, `grid_export_request`; `VARIANT_GRID_LABEL_OVERRIDES` moves to
  `analysis/grid_export.py`, its only real consumer).

## Stage F — the export path

- `ExportVariantGrid` keeps its PK-batched, contig-ordered walk (#1257) as
  `export_rows() -> Iterator[dict]`, replacing the `sort_items`/`paginate_items` overrides — with no
  paginator to subvert, the batching reads as what it is.
- `analysis/grid_export.py`: `grid.get_items(request)[2]` → `grid.export_rows()`;
  `grid.get_colmodels()` → `config.csv_columns()` (already the `{name, label}` shape
  `grid_export_csv` and `_get_colmodel_info_dict` want); `grid.get_field_names()` →
  `[rc.name for rc in config.enabled_columns]`; `grid.iter_format_items()` → `format_rows(...)`.
- `variantopedia/views.py`: `tagged_variant_export` becomes a `TaggedVariantGrid` subclass overriding
  `ordering()` to force genomic order, streamed through the shared CSV helper — the
  `_get_sidx_and_sord` override was reaching into the jqGrid protocol to say "sort by position".
- The Celery tasks are untouched: they already build a `FakeRequest`, which is all the config needs.

## Stage G — the remaining colmodel consumers, then deletion

Each of these reads two or three keys off a colmodel dict and moves to the `RichColumn` attribute:

| File | Reads |
|---|---|
| `analysis/forms/forms.py` `ColumnSummaryForm` | `label`, `name`, `index` |
| `analysis/views/views.py` `node_column_summary` | `name`, `label`, `sorttype` |
| `analysis/grids.py` `NodeColumnSummaryConfig` | `get_column_colmodel(...)["label"]`, `get_field_formatters()` |
| `analysis/views/nodes/node_view.py` | `name` (graph form + column summary form) |
| `analysis/management/commands/analysis_grid_find_float_nan_data.py` | `sorttype`, `index` |

Then: **delete `library/jqgrid/`**, remove the `LEGACY-JQGRID-ENGINE` markers, and verify
`grep -rn LEGACY-JQGRID` is clean. Update `variantgrid_to_datatables_plan.md` to record that its
deferred native rewrite landed.

## Tests

The named safety net, and where its assertions land:

- **`analysis/tests/test_jqgrid_datatable_adapter.py`** → `analysis/tests/test_variant_datatable.py`.
  The param-translation class goes away with the protocol it tested; everything else is ported onto
  the config: column shape and widths, initial order, the envelope, `recordsTotal` from the stored
  node count with no `COUNT(*)`, `draw` echoing, the definition's `gridName`/`pageLength`/`extra`/
  `downloadUrl`, the filter builder's fields and operations, and a rule that really narrows the grid.
  `test_packed_genotype_sort_index_round_trips` becomes a `sort_annotation` test — same guarantee
  (sorting by a sample's zygosity column produces a real sort), expressed against the new API.
- **`analysis/tests/test_grid_export.py`** — the CSV/VCF content, ordering, batching, per-row query
  and Celery-launch assertions are all behavioural and stay as they are; only the
  `{"_search": "true", "filters": ...}` request helper changes.
- **`analysis/tests/test_node_grid_sort_limit.py`** — `sorting_disabled()` keeps its tests;
  `_sort_items(qs, sidx, sord)` assertions become `ordering(qs)` with DataTables order params, and
  the colmodel-`sortable` ones become `orderable` on the built columns.
- **`variantopedia/tests/test_all_variants_grid.py`** — filter and approximate-count tests are
  untouched; the sort tests move onto `ordering()`, the colmodel allowlist test onto `orderable`.
- **`variantopedia/tests/test_tagged_variant_grid.py`**, **`variantopedia/tests/test_urls.py`**,
  **`genes/tests/test_urls.py`** — grid construction and URL kwargs (the `op` segment goes).
- **`library/tests/test_jqgrid_paginator.py`** → `library/tests/test_paginators.py`, unchanged.

New tests earn their keep where the rewrite creates a branch that did not exist before:
`get_known_count` short-circuiting both record counts, `max_page_length = None` serving more than 100
rows, `sort_annotation` ordering, and `filter_rules_to_q`'s `isnull`/`__in` handling (previously only
covered indirectly through a grid).

Full run: `python3 manage.py test --keepdb`, then `./scripts/linting/run_pylint.sh`.

## Sequencing

A-B-C are additive and mechanical and can land together; D-E-F are the move and want to be one
reviewable unit (the grids, their views and their export change together or not at all); G is
deletion. One PR against `jqgrid_to_datatables_1785`, commits following the stages, retargeted to
`master` once #1787 merges.

## Outcome — what actually landed

All seven stages, in one PR against `jqgrid_to_datatables_1785`. `library/jqgrid/`,
`library/django_utils/jqgrid_datatable_adapter.py` and `library/django_utils/jqgrid_view.py` are
deleted; `grep -rn LEGACY-JQGRID` and `grep -rn colmodel` are both clean outside these plan files
(the old plan keeps its own history).

Beyond the plan as written:

- **`KnownCountPaginator` deleted rather than moved.** It existed to stop `Paginator.count` running a
  `COUNT(*)`; `DatabaseTableView` pages by slicing and takes its count from `get_known_count()`, so
  there is no paginator left to subvert. Its test went with it - the behaviour it guarded is covered
  at the config level in `analysis/tests/test_grid_known_count.py`.
- **`RichColumn.sort_nulls_first`.** The jqGrid engine sorted ascending with `nulls_first=True` (nulls
  are the lowest value); `DatatableConfig` puts them last in both directions. Variant columns set the
  flag so the ordering users see is unchanged.
- **Strict JSON on every DataTables endpoint.** The adapter encoded with `allow_nan=False` so a NaN in
  a float annotation column raised server side rather than going down as a bare `NaN` token and
  breaking `JSON.parse`. `JSONResponseMixin` now does the same, so the four page grids keep that
  protection.
- **`DatatableConfig.count_unfiltered`.** `DatabaseTableView` counts the unfiltered queryset as well
  as the filtered one, so the pager can say "filtered from N". The adapter reported one number for
  both, and on the variant grids the whole-table total costs more than the page it's counting - so
  `AbstractVariantGrid` turns it off and reports the filtered count as both.
- **`VariantGrid.server_csv_download = False`.** A node export can take minutes and runs under Celery
  from its own toolbar button, so the node grid's definition carries no `downloadUrl`.
- **`get_node_sql`** (the node editor's debug tab) reached into the engine, setting `grid.fields =
  ['id']` and re-running `get_queryset`. It now shows the node's own queryset as the node SQL and the
  annotated grid queryset as the grid SQL - which is what the two labels always meant.
- **Dead client code removed with its Python authors:** `VariantGridFormat.link`, the `renderKwargs`
  plumbing it read, and `detailsLink`'s `ctx.kwargs` branch. Nothing had authored `formatter_kwargs`
  for some time.
- `grid_export_csv`'s `label_overrides` parameter went too - the one caller
  (`VARIANT_GRID_LABEL_OVERRIDES`, renaming `id` to `variant_id`) is now
  `AbstractVariantGrid.csv_columns()`, where that knowledge belongs. Its `colmodels` argument is now
  `columns`.
- **`tagged_variant_export` gained a test.** It was the one path in this rewrite with no coverage, and
  its shape changed. Writing it surfaced a latent isolation bug in that test class: a test client call
  leaves its request in the threadlocals, and `UserSettings` are cached against it, so the next test
  read settings computed before it changed them. `setUp` now clears the threadlocal request.

Net: 48 files, roughly -1200 lines. `python3 manage.py test --keepdb` green (2334 tests);
pylint 9.91/10, unchanged.

## Left for follow-up

- `DataFrameDatatableConfig` (`library/django_utils/datatable_dataframe.py`) is a second, parallel
  implementation of the definition/data/CSV protocol. It is already native and out of scope here, but
  once `DatatableConfig` has the Stage A hooks the two could share the definition builder.
- `VariantGridColumn` still stores `model_field` / `queryset_field` per column, which only exist
  because the column list is half model fields and half query annotations. Worth revisiting when the
  custom columns UI is next touched.
- `close_existing_grid()` (`snpdb/data/data.html`, `seqauto/sequencing_data.html`, called from
  `upload/upload_forms.html`) still calls `grid.jqGrid('GridUnload')` on a global that no longer
  exists. It is guarded by a `typeof` check so it is inert, and it belongs to #1787's client
  deletion rather than this change - left alone here.
