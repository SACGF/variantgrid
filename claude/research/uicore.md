# uicore — research notes

Verified against a96540a68 on 2026-09-25

uicore is the server-rendered UI kit every other app's templates lean on: the page chrome
(`uicore/templates/uicore/page/base.html`), the template tag libraries (label/value rows, modals, tabs, menus, help,
JSON-for-JS), a handful of form widgets, and the `LazyRender` pattern for editing an object in place over AJAX. It has
no models ([models map](../maps/models.md#uicore) says "(none)"), no URLs, no tasks and no signals; `commands.md` and
`signals.md` have nothing for it. The other half of the UI - behaviour wired to markup - lives in
`variantgrid/static_files/default_static/js/global.js`, and the grid engine lives in snpdb
(`snpdb/views/datatable_view.py:DatatableConfig`); `uicore/AGENTS.md` has the how-to for both, including its `## Grids`
section. This document is the why and the traps.

## Flows

### A page request

Views render a template that extends one of about ten `menu_*_base` variables rather than a fixed path:
`snpdb/processors.py:settings_context_processor` puts `menu_variants_base` etc. into every context (each resolving to
`snpdb/menu/<name>.html`, which extends the page base), so a deployment can swap the menu shell without
touching page templates. About 70 templates extend a menu base and about 100 extend the page base directly; a few menu-less pages use
`uicore/templates/uicore/page/base_external.html`, and the login pages have their own
`variantgrid/templates/default_templates/external_abstract_base.html`. The same context processor supplies
`url_name_visible` (`variantgrid/perm_path.py:get_visible_url_names`) and `form_helper`
(`uicore/utils/form_helpers.py:FORM_HELPER_HELPER`), so templates can write `{% crispy form form_helper.horizontal_nested %}`
and hide links without loading anything.

Menus are `uicore/templatetags/ui_menus.py:menu_top` and `menu_item`: each reverses its URL only when
`get_visible_url_names` says the name is registered on this deployment (`settings.URLS_NAME_REGISTER`), otherwise it
renders nothing. An unregistered name is not removed from the URLconf - `variantgrid/perm_path.py:_perm_path` wraps
the view in `require_superuser` - so the menu hides it from everyone while superusers can still reach it directly.
The per-area menu bars are inclusion tags in `uicore/templatetags/ui_menu_bars.py`.

### Markup that turns into behaviour

The base page loads jQuery, Bootstrap 4, DataTables and then `global.js`, whose `processors` list pairs a selector with
a function. It runs every processor over the document on ready, then a `MutationObserver` on `body` runs
`checkNode` over each added node recursively, marking matched nodes (and each added top-level node) `data-p` so they
are not processed twice. This is why an AJAX tab, a modal and a full page all behave the same: server code only ever
emits attributes (`data-toggle="ajax-modal"`, `table[data-datatable-url]`, `.format-json`, `[data-help]`,
`[data-preview-db]` ...). Processors with `func: null` (`[role="gridcell"]`, select2 bits) mark subtrees to leave alone
during observation.

### Tabs

`{% ui_register_tab %}` / `{% ui_register_tab_embedded %}` append a `uicore/templatetags/ui_tabs_builder.py:TabBuilderTab`
to a `TabBuilder` stored in the template context under `ui-tab-<tab_set>`; `{% ui_render_tabs %}` renders
`uicore/templates/uicore/tags/tabs.html` as tabs, an accordion or a carousel. A URL tab renders as
`a[data-href][data-toggle="tab"]` and `global.js` loads it when shown (immediately if it is the active one), so the tab
view returns a fragment, not a page. An embedded tab whose rendered content starts with `/` becomes a URL tab too
(`uicore/templatetags/ui_tabs_builder.py:LocalTabContent.render`). The active tab survives reloads through
`?activeTab=<tab_set>:<tab_id>`: global.js writes it with `history.replaceState` on `shown.bs.tab`, and
`uicore/templatetags/ui_tabs_builder.py:check_active_tab` reads it from GET, or from the Referer when the tab content
itself was fetched by AJAX, or from an `active_tab` context variable. Tab ids are compared after stripping a trailing
`_<number>` or `-tab`, so `activeTab=lab_view` still matches `lab_view_12`. `TabBuilder.tabs_required` renders a single
embedded tab straight onto the page without tab chrome.

### Editing an object in place (LazyRender)

`uicore/views/ajax_form_view.py:LazyRender` holds an object, the context name for it, one template and a
static/dynamic context. The page embeds it (`{% embed lazy %}` via `uicore/templatetags/ui_utils.py:_embed`, which
defaults to `EMBEDDED_CARD` and wraps the output in `.card.embed-wrapper`); the template shows a read-only card with an
edit link marked `data-toggle="embed-content"`, and global.js loads the edit form into the closest `.embed-wrapper` or
`.modal-content`. The form is `form[data-toggle="ajax-form"]`; global.js POSTs it and replaces the wrapper with the
response, which the view produces with `lazy_render(obj).render(request, saved=True)`. In MODAL mode a save returns a
hidden `.auto-close-modal` block with `data-replace="#<name>-<pk>"`: the `[data-replace]` processor moves it over the
inline copy on the page and `cardToModal` closes the dialog. The mode travels as `?mode=` or a `mode` context key and
reaches the template as `mode`, so one template branches between read-only card, inline and form. The only users are
the classification triage views (`classification/views/discordance_report_triage_view.py:DiscordanceReportTriageView`,
`classification/views/overlaps_view.py:TriageView`); `uicore/views/ajax_form_view.py:AjaxFormView` itself is only a
marker base with an empty `lazy_render` classmethod.

### Label/value rows, modals, help

`uicore/templatetags/ui_utils.py:LabelledValueTag` is the house replacement for a hand-written Bootstrap row: `hint`
picks the column split (default 3/9, `tiny`, `chunky`, `inline`, `large-label`), an empty body renders a grey `-`, a
bare zero gets `zero-value`, `admin_only` prefixes a key icon and hides the row from others, and
`visible_fields`/`show_if` drop the row when its key is not in a caller-supplied set. It emits a real
`<label for>` only when the body contains an input/select/textarea with an id (accessibility, b1733f35e); otherwise
the label is a `div.field-label`. `ModalTag` renders a Bootstrap modal plus its toggle link. Help comes from
`uicore/templatetags/ui_help.py:page_help_embedded` (help written inline in the template, preferred) or `page_help`,
which reads `page_help/<page_id>.html` via staticfiles finders and reports a missing file with `report_message`.
`uicore/templatetags/ui_utils.py:InstallInstructionsTag` is the superuser-only collapsible "how to install this
annotation source" block on the annotation page, open and pink only when the component is missing (e29f286b3).

### JSON on the page

Two paths. Data for JavaScript goes into a `<script>` with `{{ x|jsonify }}`
(`uicore/templatetags/js_tags.py:jsonify_for_js`: `json.dumps`, `</script>` escaped, marked safe). Data for humans goes
through `{% code_json x %}` (`uicore/templatetags/js_tags.py:code_json`), which drops the JSON into a `.format-json` div
that global.js parses and replaces with coloured, collapsible HTML (`_formatJson`). `code_json` serialises a
`uicore/json/validated_json.py:ValidatedJson` first; the `"*wrapper$": "VJ"` envelope written by
`ValidatedJson._serialize` is what `_formatJson` recognises to draw each message inline beside the value it is about.
That is the ClinVar export preview: `classification/models/clinvar_export_convertor.py` builds the submission as nested
ValidatedJson (a missing required value becomes `ValidatedJson.make_void(messages)`, dropped from the pure JSON but
kept with its error), `classification/models/clinvar_export_models.py` stores the serialised form, and the pure form
(`to_json` / `pure_json`) is what is sent.

### Widgets

`uicore/widgets/radio_other_widget.py:RadioOtherWidget` / `CheckboxOtherWidget` with `ChoiceFieldWithOther` /
`MultiChoiceFieldWithOther` are a radio/checkbox list plus a free-text "other", used by the review forms
(`review/widgets/multi_lab_selector.py`); ticking "other" with no text returns a `ValuesMissingOther`, which the
field's `validate` turns into an error. `uicore/widgets/describe_difference_widget.py:DescribeDifferenceField` is the
review app's "how do these differ" control. `uicore/widgets/date_widget.py:NativeDateInput` is a plain
`<input type="date">`, which replaced the jQuery UI datepicker (2278f950d). `uicore/templatetags/tips_tags.py:tip_box`
picks a random feature tip server-side (from `variantgrid/tips.py`, filtered by visible URLs) and hands the rest to
`tips.js`.

## Why it is shaped this way

- **Server-rendered HTML, behaviour by attribute.** The readme's stated position: VariantGrid is traditional
  Django templates, not a JS framework. Attaching behaviour by selector means a fragment returned by any AJAX view
  works without its caller knowing what is in it, and a new behaviour is one processor entry rather than an init call
  at every call site. The cost is that the wiring is invisible from the template - you find it by grepping `global.js`.
- **One template, many placements.** LazyRender exists (4dcb89f20, 2023, for discordance triage) so the same
  template shows an object read-only inline, as a card, as an edit form in place, or in a modal, with the save response
  re-rendering the same template; the alternative was a template per placement.
- **Deployment-visible URLs in one register.** Shariant, SA Pathology and variantgrid.com enable different apps. One
  `URLS_NAME_REGISTER` drives both the URL guard and every menu, tab (`url_check=True`) and grid link, so a site that
  turns off patients loses the links as well as the pages.
- **ValidatedJson** (4546736a1, 2021) keeps validation messages attached to the exact node that caused them, so the
  ClinVar preview can say "this condition has no MONDO term" next to the condition instead of in a list at the top.
- **Tags over includes.** `labelled`, `modal`, the tab tags and `if_user_can_edit` are block tags built on
  `uicore/templatetags/ui_utils.py:parse_tag` + `TagUtils` so they take keyword arguments and wrap arbitrary template
  content; the readme admits there are probably too many and that `english_tags`, `js_tags` and `ui_utils` overlap.

## History

- 2020-09-30 "blank slate": uicore arrives with the tag libraries, the tab builder and the page base.
- 2021: ValidatedJson with the ClinVar export rewrite. 2023: LazyRender/AjaxFormView (discordance triage) and the
  "other" / describe-difference widgets (review discussions).
- Aug 2026 front-end cleanup: jqGrid replaced by the DataTables engine in snpdb (a16d58a47 onwards, #1462), jQuery UI
  removed (89e53d513, #1790) with native date inputs, lodash/d3 retired, CDN files self-hosted (3ca8513ee, #1741),
  feature tips (382ef1e21, #233).
- Sep 2026: model icons move to `{% preview_icon %}` with an SVG sprite for pedigree/node shapes
  (`uicore/templates/uicore/tags/svg_icon_sprite.html`, themed via `style` not presentation attributes, 16ec67411);
  `URLS_NAME_REGISTER` gates DRF router URLs too (`variantgrid/perm_path.py:router_urls`, e3588762f); `global.scss`
  split into partials (17178d469); install-instruction ids made unique per render (e29f286b3).

## Traps

- **`jsonify` does not escape quotes in a plain string.** In `uicore/templatetags/js_tags.py:jsonify_for_js` the str
  branch does `replace('"', '\"')`, which in Python is a no-op: `{{ 'a"b'|jsonify }}` renders `"a"b"`, broken JS (and
  an injection vector for user text). Dicts and lists go through `json.dumps` and are fine; wrap a lone string in a
  dict or use `json.dumps` in the view until it is fixed.
- **`jsonify` output is only safe inside `<script>`.** `json.dumps` leaves `<` and `>` alone and the result is marked
  safe, so `uicore/templates/uicore/tags/code_block_json.html`, which puts it into a `div.format-json`, renders any
  HTML inside a JSON string as markup before global.js reads it back with `.text()`. Treat `code_json` of
  user-supplied values with suspicion.
- **`admin_only` on `ui_register_tab_embedded` crashes for non-superusers.** `LocalTabContent.render` returns `None`
  instead of `""` for a hidden tab, and Django's `NodeList.render` joins strings, so the page raises `TypeError`
  (reproduced on Django 6.1). `snpdb/templates/snpdb/labs_graph_detail.html` uses it. The body is also rendered
  before the admin check.
- **`tab_id` on `ui_register_tab` only affects active-tab matching.** The stored id is always `url + "_" + param`
  (`uicore/templatetags/ui_tabs_builder.py:ui_register_tab`), and `url` is required even though it defaults to `None`.
  Tabs registered inside a `{% for %}` or `{% with %}` vanish when the sub-context pops - call
  `{% ui_register_tabs tab_set="x" %}` at the outer level first (the error from `ui_render_tabs` says so).
- **`labelled` does not escape `label` or `help`.** Both are interpolated into an f-string (`help` only has `"`
  swapped for `'`); pass user text through `|escape` first.
- **`current_record` points at a template that does not exist** (`uicore/templatetags/ui_menus.py:current_record`
  renders a `current_record` template under uicore/menus that is not in the tree); nothing uses it yet.
- **`menu_item(href=...)` skips the URL register** - only the `url_name` path is checked. `menu_top` accepts
  `a|b` and uses the first visible name.
- **`ChoiceFieldWithOther.valid_value` is always `True`** (`uicore/widgets/radio_other_widget.py`): any posted string
  is accepted as a choice, since "other" text is legitimately anything. Validate the value downstream if it matters.
- **global.js's initial pass ignores its skip list.** The ready-time loop does `for (const badTest in
  badElementTests)`, iterating indices, so `func: null` selectors only protect nodes added later through the
  MutationObserver. The observer's `ignoreSelectors` (`.wiki-tag`) is tested against the mutation target's parents
  only, and only the top-level added node is marked `data-p`, so a child re-added elsewhere can be processed twice.
- **There is no `jsstring` filter** (commented out in `js_tags`); which library holds which tag is in `uicore/AGENTS.md`.
- **No tests in uicore.** Tag behaviour is exercised by rendering templates in other apps' tests
  (`variantgrid/tests/test_tips.py`, `analysis/tests/test_node_display.py`, ValidatedJson in
  `classification/tests/utils/test_json_utils.py`) and by URL tests; nothing tests global.js - check processors in a
  browser.
