# uicore — agent notes
Owns: template tags (labelled/modal/tabs/menus/help), js_tags JSON-for-JS filters, ValidatedJson, AJAX form embedding
  (LazyRender), "other" choice widgets, crispy FormHelper builders, page chrome (uicore/page/base.html). No models, no urls.
Start with: templatetags/ui_utils.py (labelled, modal, preview, if_user_can_edit, embed), templatetags/js_tags.py
  (jsonify, format_value, timestamp, code_json), templatetags/ui_tabs_builder.py (ui_register_tab, ui_render_tabs),
  views/ajax_form_view.py (LazyRender, AjaxFormView), json/validated_json.py (ValidatedJson, JsonMessages),
  variantgrid/static_files/default_static/js/global.js (the processors list that wires behaviour to markup).
Patterns here:
- Behaviour is attached to markup, not called: global.js keeps a `processors` list of {test: selector, func} that
  runs on document ready and again for every node the MutationObserver sees added
  (variantgrid/static_files/default_static/js/global.js:checkNode). Render the attribute/class
  (`table[data-datatable-url]`, `.format-json`, `[data-preview-db]`, `[data-toggle="ajax-modal"]`, `[data-replace]`,
  `input[type=checkbox][data-cookie]`) and it works in a page, an AJAX tab and a modal alike; add a behaviour as
  a new processor entry.
- Put Python values into a `<script>` with `{{ value|jsonify }}` (uicore/templatetags/js_tags.py:jsonify_for_js) -
  it escapes `</script>` and marks safe; `{{ value|js_symbol }}` for identifiers. Never build JS literals with `{{ }}`.
- Label/value rows are `{% labelled label="..." %}...{% endlabelled %}` (uicore/templatetags/ui_utils.py:LabelledValueTag),
  which takes hint, help, admin_only, show_if - use it rather than a hand-built Bootstrap row. New block tags parse
  kwargs with uicore/templatetags/ui_utils.py:parse_tag.
- Tabs: `{% ui_register_tab tab_set="x" label="..." url="view_name" param=obj.pk %}` then `{% ui_render_tabs tab_set="x" %}`
  (uicore/templatetags/ui_tabs_builder.py:ui_register_tab) loads the view by AJAX on first click; `url_check=True`
  hides the tab when the URL is unregistered on this deployment; `{% ui_register_tab_embedded %}` inlines content.
- Deployment-visible URLs come from variantgrid/perm_path.py:get_visible_url_names - menus
  (uicore/templatetags/ui_menus.py:menu_item) and snpdb/grids.py:url_if_visible check it; guard any link to a view
  a site unregisters (patients on Shariant) the same way.
- An object viewed inline, as a card, or edited in a modal is one template rendered through
  uicore/views/ajax_form_view.py:LazyRender - `{% embed lazy_render %}` in the page (uicore/templatetags/ui_utils.py:_embed),
  `.render(request, saved=True)` from the POST view; its wrappers (`embed-wrapper`, `data-replace`,
  `auto-close-modal`) are what global.js swaps.
- Hover cards: implement library/preview_request.py:PreviewModelMixin on the model (preview_icon, preview_with) and
  render `{% preview obj %}` / `{% preview_icon obj %}`; other apps add rows with
  `@receiver(preview_extra_signal, sender=Model)` (library/preview_request.py:preview_extra_signal). Text mentions
  (PMID:, OMIM:) answer library/preview_request.py:preview_request_signal, served by
  library/preview_request.py:preview_view (url name `preview_data`, registered in snpdb/urls.py).
- ValidatedJson carries messages beside data: build with uicore/json/validated_json.py:JsonMessages.error/warning,
  read `all_messages`/`has_errors` (uicore/json/validated_json.py:ValidatedJson). `to_json()`/`pure_json()` drop
  the messages, `serialize()`/`deserialize()` keep them - store the serialized form, send the pure form.
Gotchas:
- `{% load %}` names matter: labelled/modal/preview/badge/boolean live in `ui_utils`;
  jsonify/format_value/timestamp/code_json/duration/get_item live in `js_tags` (uicore/templatetags/js_tags.py:format_value).
  The wrong library only fails when the page renders.
- The MutationObserver marks processed top-level nodes `data-p` and skips subtrees under `ignoreSelectors` (third
  party widgets) - a processor that never fires on injected HTML is usually inside one
  (variantgrid/static_files/default_static/js/global.js:checkNode).
- Preview handlers run under send_robust: an exception becomes a report_message and the card shows nothing
  (library/preview_request.py:PreviewRequest.preview_data).
Tests: no uicore/tests. Tag logic is tested by rendering a `Template("{% load x %}...")` with a Context
  (variantgrid/tests/test_tips.py) or `render_to_string` of the template (analysis/tests/test_node_display.py);
  ValidatedJson in classification/tests/utils/test_json_utils.py. Pages are covered by URL tests:
  library/django_utils/unittest_utils.py:URLTestCase - `_test_urls` for pages, `_test_datatable_urls` for grid
  endpoints (hits the URL and `?dataTableDefinition=1`). No JS tests - check processors in a browser (`run` skill).
Deep reference: __uicore_readme.md · claude/research/uicore.md

## Grids
- One engine: subclass snpdb/views/datatable_view.py:DatatableConfig, set `self.rich_columns` (list of
  snpdb/views/datatable_view.py:RichColumn) in `__init__`, implement `get_initial_queryset`, and put per-request
  narrowing in `filter_queryset` reading `self.get_query_param(...)` (URL kwargs first, then GET/POST). Worked
  example with tests: snpdb/tests/test_datatable_server_csv.py:CohortCsvColumns.
- Wire it with `path('.../datatable', DatabaseTableView.as_view(column_class=XColumns), name='x_datatable')`
  (snpdb/views/datatable_view.py:DatabaseTableView). That one endpoint answers `?dataTableDefinition=1` (columns),
  the DataTables row requests, and `?dataTableCsv=1` (server CSV).
- Mount it with markup only: `<table data-datatable-url="{% url 'x_datatable' %}" data-datatable-data="jsFunc"></table>`;
  the global.js processor builds a variantgrid/static_files/default_static/js/datatable_definition.js:DataTableDefinition
  from it. `data-datatable-data` names a JS function whose return object is merged into every request.
- RichColumn: `key` is a `.values()` path and the sort/search key; `name` is the JSON field (no dots); `renderer`
  is a server callable taking snpdb/views/datatable_view.py:CellData (list `extra_columns` it reads);
  `client_renderer` is a JS expression, usually `TableFormat.*`
  (variantgrid/static_files/default_static/js/datatable_definition.js:TableFormat) or a function in
  datatables_client_renderers.js; `visible=False` sends data without a column; `detail=True` moves it to the expand row.
- Sorting: orderable needs a `key` or `sort_keys`; `default_sort=SortOrder.ASC` picks the initial order, else the
  first enabled column sorts (snpdb/views/datatable_view.py:DatatableConfig.initial_order); a pk tiebreaker is
  always appended (`_get_sort_tiebreaker` - override on grouped querysets); `null_order` places NULLs.
- Search box, CSV button, filter builder (`filter_builder` + `column_filter=FilterField(...)`), row expansion
  (`expand_client_renderer`) and per-user page length (`grid_name` -> snpdb/models/models_user_settings.py:UserGridConfig)
  are class flags; the definition JSON from snpdb/views/datatable_view.py:datatable_definition is all the client reads.
- CSV follows the config: `server_csv_download = True` streams `csv_columns()` over `iter_export_rows()`
  (snpdb/views/datatable_view.py:DatatableConfig.export_columns) with the grid's filters and ordering; set `csv_rendered=True`
  where the server renderer changes the value, `include_in_csv=False` on action columns. The client CSV button only
  exports rows already fetched.
- COUNTs cost: `count_unfiltered = False` skips the unfiltered one; `known_count()` / `approximate_count_enabled`
  supply a stored or estimated total (snpdb/views/datatable_view.py:datatable_response).
- Variant grids differ: snpdb/grids.py:AbstractVariantGrid builds `rich_columns` per user from the UserSettings
  CustomColumnsCollection via snpdb/grid_columns/custom_columns.py:get_variant_grid_columns (VariantGridColumn
  catalogue, composite cells, columns dropped when this build's annotation version never annotates them) plus
  snpdb/grids.py:get_standard_overrides for renderers. They use `ajax_type = 'GET'` + `cache_stable_params` for
  `@cache_page`, `table-layout: fixed` so every column needs a `width` (`default_column_width` fills in), sort the
  Variant column genomically (`_genomic_order_by`) and filter `variantallele` to the grid's build so shared contigs
  (MT) do not duplicate rows.
- The analysis node grid skips DatabaseTableView and calls snpdb/views/datatable_view.py:datatable_response
  directly; grid tests do the same on a RequestFactory request (analysis/tests/test_grid_export.py:GridExportTestCase).
- Filter rules become Django lookups on the column key via library/django_utils/filter_rules.py:rules_to_q; they are
  applied once in `apply_filters` - repeating them in `get_initial_queryset` adds a second JOIN and multiplies rows.
