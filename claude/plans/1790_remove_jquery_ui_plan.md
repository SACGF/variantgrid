# #1790 — retire jQuery UI in favour of Bootstrap 4

[#1790](https://github.com/SACGF/variantgrid/issues/1790) asks for a cut-down jQuery UI build
containing only the nine widgets we call. This plan goes further: it moves every one of those widgets
onto Bootstrap 4, native HTML inputs, or a small MIT-licensed library, and deletes jQuery UI.

The site already runs Bootstrap 4 everywhere. Every page currently loads both frameworks, so buttons,
tabs, dialogs and accordions each have two competing appearances and two APIs, and several templates
apply both at once (`cohort_hotspot.html:58` renders `class='btn btn-primary'` and then calls
`.button()` on it). Consolidating on Bootstrap removes 255KB of JS and 31KB of CSS (the CSS is
vendored three times, once per site theme), and leaves one set of UI idioms.

## End state

- `uicore/templates/uicore/page/base.html`, `base_external.html` and
  `analysis/templates/analysis/analysis_editor_and_grid.html` load Bootstrap and no jQuery UI.
- `variantgrid/static_files/{default,runx1,shariant}_static/js/lib/jquery-ui-1.13.2.custom/` is gone.
- Two new vendored libraries, both MIT (VariantGrid ships under BSL 1.1 with an AGPL v3 change
  licence, so vendored JS stays permissive): html5sortable and noUiSlider.
- No `ui-*` class names remain in our templates, JS, SCSS or CSS.

## Inventory

Every jQuery UI entry point in the codebase, and what replaces it.

| Widget | Call sites | Replacement | New dependency |
| --- | --- | --- | --- |
| `tabs` | 24 files (15 load panes over AJAX) | `nav-tabs` + `tab-pane` + `data-toggle="tab"`, reusing `loadAjaxTab()` in `global.js` for AJAX panes | — |
| `button` | 33 calls in 21 files | `class="btn btn-*"` in the markup | — |
| `dialog` | 6 files | `createModal()` / `createModalShell()` / `loadAjaxModal()`, already in `global.js` | — |
| `accordion` | 4 calls in 3 templates | `.collapse` + `data-parent` | — |
| `datepicker` | 6 sites | `<input type="date">` | — |
| `slider` (single handle) | 5 sites | `<input type="range">` | — |
| `slider` (two handles) | 1 site | noUiSlider | noUiSlider 15.8.1, 27.7KB + 4.2KB CSS |
| `sortable` + `disableSelection` | 6 calls in 4 templates | html5sortable | html5sortable 0.14.0, 16.7KB |
| `selectable` | 1 site | ~50 lines of marquee code in `analysis_nodes.js` | — |
| `selectmenu` | 1 custom widget subclass | Bootstrap dropdown driving the underlying `<select>` | — |

Net: −286KB vendored, +49KB vendored.

## Libraries to vendor

Both go in `variantgrid/static_files/default_static/js/lib/`, matching the existing per-library
directory layout, and load from the templates that need them rather than from `base.html`.

**[html5sortable](https://github.com/lukasoppermann/html5sortable) 0.14.0, MIT, 16.7KB minified.**
Written as a jQuery UI `sortable` replacement, and the option and event names line up with what we
already pass: `items`, `placeholder`/`placeholderClass`, `acceptFrom` in place of `connectWith`, and
it emits a `sortupdate` event — the exact event all three of our sortable pages bind. It uses the
native HTML5 drag-and-drop API; add the DragDropTouch polyfill alongside it if touch support for
these admin-style pages is wanted.

**[noUiSlider](https://refreshless.com/nouislider/) 15.8.1, MIT, 27.7KB + 4.2KB CSS.** Built
specifically as a jQuery UI slider replacement, with first-class two-handle support. Needed only by
the allele frequency filter; the other four sliders become native range inputs.

## Phase 1 — buttons — landed in f68444579

The largest count and the simplest change, and it makes the rest of the work visible as it lands.

For each of the 33 `.button()` / `.button('enable')` / `.button('disable')` calls
(`grep -rn "\.button(" --include=*.js --include=*.html`, excluding `js/lib/`):

- Give the element Bootstrap classes in the markup — `btn btn-primary` for the primary action of a
  form, `btn btn-outline-secondary` otherwise. Several already carry them and can simply lose the JS
  call.
- Replace `.button('disable')` / `.button('enable')` with `.prop('disabled', true)` / `.prop('disabled', false)`.
  Sites: `create_classification_from_hgvs.html`, `create_classification_for_variant.html`,
  `cohortnode_editor.html`, `view_vcf.html:96`.
- `analysis.js:780` passes `{icon: "ui-icon-arrowthick-1-n"}`; render a Font Awesome icon inside the
  button instead, matching the `<i class="fas fa-...">` usage elsewhere in that file.

`analysis_output_node_downloads.html:7` uses the `ui-icon` sprite class — give it the Font Awesome
icon that matches its neighbours.

## Phase 2 — tabs

### Markup

jQuery UI reads `<div id="x-tabs"><ul><li><a href="#pane">Label</a></li></ul><div id="pane">…</div></div>`.
Bootstrap wants:

```html
<div id="x-tabs">
  <ul class="nav nav-tabs" role="tablist">
    <li class="nav-item"><a class="nav-link active" data-toggle="tab" href="#pane">Label</a></li>
  </ul>
  <div class="tab-content">
    <div class="tab-pane fade show active" id="pane">…</div>
  </div>
</div>
```

### The Bootstrap tab convention already in the codebase

`uicore/templatetags/ui_tabs_builder.py` and its `uicore/templates/uicore/tags/tabs.html` already
render Bootstrap 4 tabs — `ul.nav.nav-tabs` / `a.nav-link[data-toggle="tab"]` / `div.tab-pane` — and
36 templates use it. It is the target markup for this phase; match it rather than inventing a second
shape. Note it also puts `data-tab-set` on both the link and the pane, which drives the `activeTab`
URL-state tracking in `global.js:238-250`.

The 24 files still calling `.tabs()` are a separate, older set and are what this phase converts.

### AJAX panes

`global.js` already has the shim this phase needs. `loadAjaxTab(tab)` (`global.js:708`) reads
`data-href` off the link, resolves it relative to the page when it has no leading slash, and loads it
into `$(tab.attr('href'))` via `loadAjaxBlock()` (which supplies the loading overlay and spinner). It
caches on `tab.data('loaded')`. Two entries in the `enhanceAndMonitor` autoloader list wire it up:

```js
// global.js:252-258
{test: '.nav-tabs a.active[data-href][data-toggle="tab"]', func: (node) => {loadAjaxTab(node);}},
{test: '.nav-tabs a[data-href][data-toggle="tab"]',
 func: (node) => {node.on('shown.bs.tab', function(e) {loadAjaxTab($(this));});}},
```

That covers first-show loading, the initially-active pane, and dynamically added tabs, so converting
an AJAX tab is a markup change only: move the Django URL to `data-href`, point `href` at an empty
`<div class="tab-pane" id="…">`, and give the link `data-toggle="tab"`.

15 templates give the tab link a Django URL instead of a fragment, and jQuery UI fetches the pane on
activation: `view_qc`, `view_sequencing_run`, `view_gene_list`, `view_genomic_intervals`, `view_vcf`,
`view_sample`, `view_trio`, `view_quad`, `view_cohort`, `view_custom_columns`,
`view_tag_colors_collection`, `view_pedigree`, `view_ped_file`, `view_candidate_search_run`,
`view_karyomapping_analysis`.

`view_vcf.html:264`, `view_custom_columns.html:103` and `view_tag_colors_collection.html:59` pass
`{cache: true}`, which is `loadAjaxTab`'s existing behaviour.

The one behaviour `loadAjaxTab` lacks is refetching on every activation. `data.html:16`, `maps.html:20`
and `upload.html:348` currently empty sibling panes inside `activate` to force a refetch — add a
`data-reload` opt-out to the `tab.data('loaded')` check in `loadAjaxTab`, give those three links
`data-reload`, and drop the `.ui-tabs-panel` emptying.

### Programmatic API

| Current | Bootstrap |
| --- | --- |
| `tabs('option', 'active')` | index of `.nav-link.active` within the nav |
| `tabs('option', 'active', i)` | `$('.nav-link', tabs).eq(i).tab('show')` |
| `tabs({active: 1})` | `active` class on that link and pane in the markup |
| `tabs({activate: fn})` | `.on('shown.bs.tab', fn)` |
| `tabs({disabled: [1]})` | `disabled` class on the link and drop its `data-toggle` |

Sites: `analysis.js:310-313` (also update the `:has(> ul.ui-tabs-nav)` guard to test for
`ul.nav-tabs`), `base_editor.html:160-190`, `grid_editor.html:20`, `tagnode_editor.html:45`,
`view_cohort.html:231-292`, `maps.html:9`.

### Vertical tabs

`base_editor.html:190` adds `ui-tabs-vertical ui-helper-clearfix` to `#node-summary-tabs`, styled by
`variantgrid/static_files/default_static/css/vertical_tabs.css`. Bootstrap does this with
`nav flex-column` on the nav and a flex row on the container; rewrite `vertical_tabs.css` against the
Bootstrap class names, and rewrite the `#node-editor-tabs > ul.ui-tabs-nav` rules in
`analysis_nodes.css:883-894` to match.

`global.css:1216` / `global.scss:1449` and `flags.css:346` / `flags.scss:350` also style
`.ui-tabs`; port those to `.nav-tabs`. Edit the `.scss` and hand-apply the same change to the `.css`,
per the SCSS convention in CLAUDE.md.

## Phase 3 — dialogs

`global.js` already has `createModal()`, `createModalShell()`, `cardToModal()` and `loadAjaxModal()`,
and 11 other files already use Bootstrap modals, so this phase is a port onto existing helpers.

- `global.js:1342` `showReloadPageErrorDialog` — already carries a `// FIXME turn into Bootstrap modal`
  comment. Build the shell with `createModalShell()`, add the Reload Page and (conditionally) Close
  buttons to the footer, and give it `data-backdrop="static" data-keyboard="false"` where the current
  code passes `dialogClass: "no-close"`. Drop `.no-close .ui-dialog-titlebar-close` from
  `analysis_nodes.css:733`.
- `grid.js:213` error dialog — same helper.
- `create_patient_dialog.html`, `view_vcf.html:152-196` (new project dialog) — these are form dialogs
  with an OK/Cancel button pair and `dialog("open")` / `dialog("close")` calls; `.modal('show')` and
  `.modal('hide')` on a `createModal()` shell. Their `ui-state-error` validation classes become
  `is-invalid`, and the `ui-widget-content ui-corner-all` inputs (`view_vcf.html:723-726`) become
  `form-control`.
- `quad_wizard.html:102`, `trio_wizard.html:124` — form-submitted confirmation, same helper.
- `map_matrix.html:41` loads a URL into a dialog — `loadAjaxModal()` is the existing pattern. Drop the
  `.ui-dialog` rules at `map_matrix.html:112-116`.

## Phase 4 — accordion

Bootstrap accordions are `.card` headers with `data-toggle="collapse"` and `data-parent` on the
panels.

- `base_editor.html:55-71` `accordionForm()` — keeps the same job of writing the open panel index to
  `#id_accordion_panel`; read it from `shown.bs.collapse` and set the initial `show` class from the
  stored value.
- `moinode_editor.html:72-100` — `accordion({active: 0})` / `({active: 1})` become
  `.collapse('show')` on the target panel. `ui-state-disabled` on the patient control header becomes
  Bootstrap's `disabled` class plus `pointer-events: none`.
- `allvariantsnode_editor.html:6` — plain markup conversion.
- `vc_form.js:1433` holds a commented-out `accordion('instance')` line; remove it with the rest.

## Phase 5 — datepicker

All six sites format as `yy-mm-dd`, which is what `<input type="date">` submits, so Django form
parsing is unaffected. `changeYear`/`yearRange: "-120:+0"` becomes `min`/`max` attributes on the
input.

- `global.js:301` — the `.date-picker` autoloader; set `type="date"` on matching inputs, or better,
  set the widget input type in the Django forms that render `.date-picker` and retire the autoloader
  entry.
- `vc_form.js:1933` — classification form date fields.
- `view_patient_specimens.html:15`, `view_specimen.html:8`, `view_patient_extractions.html:15`,
  `view_extraction.html:8`.

`global.js:231` excludes `.ui-datepicker-prev` / `.ui-datepicker-next` from the tooltip autoloader;
that selector simplifies once the widget is gone.

## Phase 6 — sliders — landed in 7ec48e88b

### Single handle

`analysis.js:1128` `setupSlider(inputSelector, sliderSelector)` already reads `min`, `max`, `step` and
`value` off a hidden `<input>` and writes the value back to it on change, so it converts to
`<input type="range">` almost directly: give the range input those attributes, bind `input` for the
live label and `change` for the value write-back, and keep the `.min-value` / `.max-value` /
`.slider-value` label updates as they are. A blank field means "not set", so park the handle at `min`
and let only the events write back, keeping the field blank until the user drags. That covers
`conservationnode_editor.html`, `damagenode_editor.html` (24 markup sites across the
`columns_version` branches) and `builtinfilternode_editor.html:76` — the COSMIC count slider, live
whenever the built-in filter is set to COSMIC. Their `slider("enable")` / `slider("disable")` /
`slider('value', x)` calls become `prop('disabled', …)` and `val(x)` on the range input.

Two call sites need more than a mechanical port:

- `conservationnode_editor.html:104` sets `slider('option', 'change', …)` after `setupSlider`, which
  *replaced* setupSlider's own change callback so the master's write-back came only from `slide`.
  With a range input both handlers coexist — bind `masterSliderChanged` with `.on('change', …)` and
  the individual sliders still follow the master on release.
- `damagenode_editor.html` `addCutoffBands()` prepends absolutely positioned band divs *inside* the
  `.ui-slider` div, which an `<input>` cannot hold. Wrap each damage slider in
  `<span class="slider-track">` and prepend the bands into that. A native range's track is opaque, so
  the bands move from behind it (`z-index: 0`, full height) to over it (`z-index: 1`,
  `pointer-events: none`, a centred 8px bar) to stay visible.

`hotspot_graph.html:319` builds its own single-handle slider with a `<label>` per step; keep the label
loop, driving it from `GNOMAD_PERCENT.length` rather than `data().uiSlider.options`. The labels move
from inside the slider element to the `.gnomad-af-slider` wrapper, which then needs
`position: relative` and an explicit `top` (the old `margin-top` was relative to a 12px jQuery UI div).

`analysis_nodes.css:546` `#node-size-slider` matches nothing in any template or JS — delete it.

### Two handles

`allele_frequency_tag.html:68` is the one two-handle slider. noUiSlider:

```js
noUiSlider.create(sliderElement, {start: [value1, value2], connect: true, range: {min: 0, max: 1}, step: 0.01});
```

`slider('values', 0)` / `slider('values', 1)` become `sliderElement.noUiSlider.get(true)`, which
returns numbers so the values post as floats, and the `slide` callback that drives
`afSetSliderMinMaxLabels` becomes noUiSlider's `update` event — which also fires on create, so it
covers the initial labels too. `getAlleleFrequencyValues()` iterates `$(".af-slider")` to collect
pairs — keep that shape and read each element's `noUiSlider.get(true)`. noUiSlider wants the element
in the DOM before `create()`, so `createSlider()` and `addSlider()` fold into one function.

### Loading noUiSlider

`analysis/templates/analysis/analysis_includes.html` is the include shared by `analysis.html`'s
`{% block head %}` and `analysis_editor_and_grid.html`'s `{% if stand_alone %}` block, so one
`<link>`/`<script>` pair there reaches both editor homes with no dependency on `globalSetup()`.

## Phase 7 — sortable

Six calls across four templates, all lists of `<li>` with a placeholder, three of them connected
pairs. html5sortable:

```js
sortable('#my_columns_sortable', {items: 'li.user-column', placeholderClass: 'sort_placeholder',
                                  acceptFrom: '#available_columns_sortable'});
```

The existing `.bind("sortupdate", …)` handlers stay as they are — html5sortable fires the same event
name. Sites: `view_custom_columns.html:126,141`, `analysis_settings_node_counts_tab.html:71,80`,
`cohort_sort.html:46`, `view_tag_colors_collection.html:61`.

`.disableSelection()` becomes `user-select: none` in CSS on those lists. `ui-state-disabled` on the
placeholder `<li>`s (`view_custom_columns.html:213,227`) and the `ul.ui-sortable` rule
(`analysis_settings_node_counts_tab.html:55`) get project-owned class names.

`view_custom_columns.html:93` calls `sortable("widget")` to reach the list element; that is the
element itself with html5sortable.

## Phase 8 — analysis canvas

Two pieces, both in the analysis UI, and the only two that are written rather than ported.

### Marquee selection

`analysis_nodes.js:819` `setupNodeModifications()` uses `selectable({filter: 'div.window', start: fn,
cancel: '.cancel, div.window'})` so a drag on empty canvas rubber-bands a set of nodes. The
integration surface is one constant — `analysis_nodes.js:2`, `const ACTIVE_CLASS = "ui-selected"` —
and one rule, `analysis_nodes.css:541` `.ui-selected.design-a-node`.

Write a marquee in `analysis_nodes.js`: on `mousedown` on `#analysis-container` where the target is
the container itself, draw an absolutely positioned `<div class="node-marquee">`, track it on
`mousemove`, and on `mouseup` toggle `ACTIVE_CLASS` on every `div.window` whose bounding box
intersects it. Keep the existing `start` behaviour — calling `replaceEditorWindow()` and
`loadGridAndEditorForNode()` when a new drag begins. Rename `ACTIVE_CLASS` to a project-owned value
such as `"node-selected"` and update `analysis_nodes.css:541`; `syncDragSelection()` reads the
constant, so it follows automatically.

### Node type picker

`analysis.js:640-672` `setupNodeTypeSelect()` subclasses `$.ui.selectmenu` via the jQuery UI widget
factory to render each option of `#id_node_types` as an icon plus a label. Replace it with a
Bootstrap dropdown built from the same `<select>`: keep `renderNodeTypeItem()` unchanged, render one
`.dropdown-item` per option and one for the button face, and on click write the value back to the
`<select>` and trigger `change` so existing listeners keep working. Style rules
`analysis_nodes.css:139` (`#id_node_types-button.ui-state-active`) and `:168`
(`.node-type-menu .ui-state-active`) get the Bootstrap `.active` / `.show` class names.

`analysis_nodes.css:567` styles `table.grid tr.ui-state-highlight`; give it a project-owned class.

## Phase 9 — removal

1. Drop the three `<link>`/`<script>` triples: `base.html:113-115`, `base_external.html:71-73`,
   `analysis_editor_and_grid.html:7-14`.
2. Delete `js/lib/jquery-ui-1.13.2.custom/` from `default_static`, `runx1_static` and
   `shariant_static`.
3. Add the two new libraries to whichever templates use them: `js/lib/nouislider-15.8.1/` is already
   loaded from `analysis_includes.html` for `allele_frequency_tag.html`; add `js/lib/html5sortable/`
   to the four sortable templates.
4. Confirm `grep -rn "ui-tabs\|ui-dialog\|ui-accordion\|ui-slider\|ui-widget\|ui-state\|ui-icon\|ui-helper\|ui-corner\|ui-sortable\|ui-selected\|ui-button\|ui-datepicker\|jquery-ui"` over templates, JS,
   SCSS and CSS (excluding `sitestatic/`) comes back clean.
5. `help_tags.html:53` (`ui-widget`) and `external_help.html:61` (`ui-widget-header`) take
   project-owned classes with the same visual result.

## Testing

Each phase is independently shippable, so test as they land rather than at the end.

- `python3 manage.py test --keepdb` covers the URL smoke tests, which will catch a template that fails
  to render after a markup rewrite.
- The AJAX tab shim is the piece with the most behavioural surface. Walk every tab of `view_sample`,
  `view_vcf`, `view_cohort`, `view_qc` and `view_sequencing_run` and confirm each pane populates,
  including the one active on page load.
- The analysis editor exercises tabs, accordion, sliders, marquee selection and the node type picker
  in one page — open a saved analysis, add a node of each type from the picker, rubber-band a group of
  nodes and drag them, and open the conservation and damage node editors.
- `view_custom_columns` and the node counts settings tab cover connected-list sorting and the
  `sortupdate` POST-back.
- Bootstrap 4 uses `data-toggle` / `data-target`, matching the convention already recorded in
  CLAUDE.md.

## Interaction with #1790 as written

The issue's cut-down jQuery UI download stands on its own and ships the bulk of the size win in a
single commit with no template changes. Taking it first is compatible with this plan: re-download
after Phase 8 with only the widgets still in use, then Phase 9 removes the file entirely. The Plotly
half of #1790 is independent of all of this.
