# Defer expensive grid loads (issue #51)

Goal: clicking a node should be **snappy** and must never auto-run an expensive
variant query. People click merge nodes fed by 2× VCF samples with ~6M variants
just to *edit* settings, and today that fires the full row query immediately,
hammering the DB.

The fix: **don't auto-load a node's variant rows when the node is large.** Show
the node's (already-cached) count plus a "Load variants" button and a CSV
download, and only run the row query when the user explicitly asks. Small nodes
behave exactly as they do today.

This subsumes the original "make the editor snappy" framing of #51: if the data
doesn't auto-load, the editor is trivially instant. This is mostly a UI-layer
change; the server-side scheduling/gating already in place is preserved.

## The central mechanism: the grid already loads in two phases

`setupGrid` (grid.js:699) makes **two** requests, chained:

1. **Config** → `node_grid_config` (`NodeGridConfig`, views_grid.py). Returns the
   jqGrid **colModel** + `postData` + caption. This is **cheap** — it builds
   column definitions only; no variant query.
2. **Data** → `node_grid_handler` (`NodeGridHandler.get` → `grid.get_data`). This
   is the **expensive row query** — the thing that melts down on a 6M-variant
   merge node.

Today, `grid.jqGrid(data)` (grid.js:751) initialises the grid from the config
and *immediately* triggers phase 2. **The whole feature is: stop auto-firing
phase 2 for large nodes.** Load config, render a placeholder, and fire the data
fetch only on demand.

### Why this is clean

- **The count is free.** `node.count` is an `IntegerField` (analysis_node.py:94)
  populated by the Celery node-update pass, entirely independent of the grid. So
  we can render *"6,000,000 variants — click to load"* without running anything.
  The headline number costs nothing.
- **FilterNode keeps working.** `filternode_editor.html:12-25` reparents
  jqGrid's `searchGrid` dialog into the editor form — but that dialog is built
  from the **colModel** (phase 1), not the row data (phase 2). As long as we
  still fetch config, the filter-builder UI works with zero expensive query.
  This also dissolves the race condition the issue itself describes (editor
  waiting on the grid): the editor only ever needed the cheap config.
- **CSV/VCF export already exists and only needs config.** `export_grid`
  (grid.js:5) builds its download URL from `grid.jqGrid('getGridParam',
  'postData')` — which comes from the config, not the data fetch — and hits
  `node_grid_export` (async via `CachedGeneratedFile`). So "a CSV you can
  download" is already built; it just needs surfacing in the placeholder.

## Decision: threshold auto-load

A node auto-loads its rows when `node.count` is below a configurable threshold;
otherwise it shows the placeholder. Small nodes feel exactly like today; only the
heavy ones require a click. The threshold applies **uniformly to all node types**
(including tag/selection nodes — nobody clicks 50k checkboxes, so they're
effectively always under the threshold anyway).

- `node.count < threshold` → behave as today (config → data, auto).
- `node.count >= threshold`, **or `node.count is None`** (not yet computed) →
  render placeholder, load **config only**, defer the data fetch to a button.
- **Already loaded this session** → auto-load regardless of count (see "Session
  re-display" below): the result is in the page cache, so re-showing is instant.

The decision is made **server-side** in `node_data_grid.html`, which is rendered
with `node` in context — so `node.count` is known at render time and no extra
round-trip is needed.

### Threshold = cascading setting (Global → Org → Lab → User)

The threshold is a **cascading user setting**, not just a global constant, so a
lab or individual user can tune it. This reuses the existing `SettingsOverride`
cascade (snpdb/models_user_settings.py:149) — Global → Org → Lab → User, where
later non-null levels override earlier ones.

1. Add a nullable field to `SettingsOverride`:

   ```python
   node_grid_auto_load_max_variants = models.IntegerField(
       null=True, blank=True,
       help_text="Analysis nodes with at least this many variants don't auto-load "
                 "their grid — the user clicks 'Load variants' to run the row query. "
                 "Blank inherits the next level up.")
   ```
   (DB migration; field is inherited by `GlobalSettings`,
   `OrganizationUserSettingsOverride`, `LabUserSettingsOverride`,
   `UserSettingsOverride` automatically.)

2. Add `node_grid_auto_load_max_variants: Optional[int]` to the `UserSettings`
   dataclass annotations (models_user_settings.py:357) so
   `dataclasses.fields(UserSettings)` picks it up — the `get_for` merge
   (line 447-459) then cascades it with no further code.

3. **Code default** still lives in `default_settings.py` as the ultimate
   fallback when the whole cascade is null:

   ```python
   # Fallback when no Global/Org/Lab/User override is set. None = always auto-load.
   ANALYSIS_NODE_GRID_AUTO_LOAD_MAX_VARIANTS = 50_000
   ```

   Resolution: `UserSettings.get_for_user(user).node_grid_auto_load_max_variants`
   `or settings.ANALYSIS_NODE_GRID_AUTO_LOAD_MAX_VARIANTS`.

4. The new field needs adding to whichever settings forms list fields explicitly
   (Global/Org/Lab/User settings forms in `snpdb/forms.py` + their templates) so
   it's editable in the UI.

### Session re-display (already-loaded nodes)

If the user already loaded a large node's grid this session, re-show it
automatically rather than re-prompting — the row query result is in the
`cache_page` server cache, so it returns instantly. Track loaded
`(node, version)` client-side in the analysis window and let it force auto-load:

- On a successful data load, `gridComplete` records the node-version:
  `getAnalysisWindow().loadedGridVersions[node_id] = node_version`.
- In `node_data_grid.html` docready, before deciding to defer, check that set —
  if this `(node_id, node_version)` is present, treat as auto-load even when
  `count >= threshold`.

Editing a node bumps its `version`, so an edited node is a genuinely new query
(new cache key) → not in the set → placeholder again, which is correct.

## Why we don't also need to cancel previous Celery jobs

The issue floats terminating previous jobs. Not needed — three layers already
dedup, and with deferral the expensive query mostly doesn't run at all:

1. **Celery scheduling** — `NodeTask` has `unique_together = (node, version)`
   and is created via `bulk_create(..., ignore_conflicts=True)`, so duplicate
   scheduling for the same `(node, version)` is silently dropped.
2. **`node_load`** (`@never_cache`) is idempotent: it inspects status and
   redirects to `node_data_grid` / `node_async_wait` / `node_errors`.
3. **`NodeGridHandler.get`** takes a per-`(user, URL+params)` `cache.add` lock
   (views_grid.py:54-82); concurrent callers sleep 2s and redirect to the same
   URL, which is `cache_page(WEEK_SECS) + vary_on_cookie` — so the redirect
   resolves from cache. It's also wrapped in `major_operation(user,
   "node_grid")` to cap concurrent expensive queries per user. Spam-clicks don't
   multiply DB work.

So "finish the existing task, cache it, next visit is free" already holds — and
deferral means we usually never start the heavy query in the first place.

## Proposed change

### 1. Setting

Add the cascading `node_grid_auto_load_max_variants` field + the
`ANALYSIS_NODE_GRID_AUTO_LOAD_MAX_VARIANTS` code default (see "Threshold =
cascading setting" above). In the `node_data_grid` view compute `grid_auto_load`
server-side:

```python
max_variants = (UserSettings.get_for_user(request.user).node_grid_auto_load_max_variants
                or settings.ANALYSIS_NODE_GRID_AUTO_LOAD_MAX_VARIANTS)
grid_auto_load = (max_variants is None) or (node.count is not None and node.count < max_variants)
```

(The session-re-display override is applied client-side in the template, since
that state lives in the analysis window.)

### 2. `node_data_grid.html` — branch on `grid_auto_load`

Today (node_data_grid.html:60-69) the docready unconditionally calls
`setupGrid(...)`, which loads config then data. Split into:

- **Always** load the editor (`load_node_editor`) and the grid **config** (so
  colModel exists → FilterNode + CSV/VCF work).
- **Conditionally** fetch data:
  - `grid_auto_load` true (or this node-version already loaded this session) →
    fetch immediately (as today).
  - otherwise → render the placeholder block; wire the **Show Grid** button to
    fire the data fetch on click.

```diff
 $(document).ready(function() {
     const node_view_url = "{% url 'node_view' ... %}";
-    load_node_editor(node_view_url);
+    const unique_code = "{{ node_id }}_{{ node_version }}";
+    load_node_editor(node_view_url, unique_code);

     const grid_url = "{% url 'node_grid_config' ... %}";
-    const unique_code = "{{ node_id }}_{{ node_version }}";
-    setupGrid(grid_url, {{ analysis_id }}, {{ node_id }}, {{ node_version }}, unique_code, gridComplete, gridLoadError, on_error_function);
+    // Re-show instantly if this exact node-version was already loaded this session (page cache hit).
+    const aWin = getAnalysisWindow();
+    aWin.loadedGridVersions = aWin.loadedGridVersions || {};
+    const alreadyLoaded = aWin.loadedGridVersions[{{ node_id }}] === {{ node_version }};
+    const autoLoad = {{ grid_auto_load|yesno:"true,false" }} || alreadyLoaded;
+    setupGrid(grid_url, {{ analysis_id }}, {{ node_id }}, {{ node_version }}, unique_code,
+              gridComplete, gridLoadError, on_error_function, autoLoad);
+    if (!autoLoad) {
+        $("#load-variants-{{ node_id }}").click(function() {
+            $("#grid-placeholder-{{ node_id }}").hide();
+            loadNodeGridData({{ node_id }}, unique_code);  // fires deferred phase 2
+        });
+    }
 });
```

And `gridComplete` records the load so a later revisit re-shows automatically:

```diff
 function gridComplete() {
     const unique_code = "{{ node_id }}_{{ node_version }}";
     if ($("#" + unique_code, "#node-data-container").length === 0) { return; }
     ...
     registerComponent(unique_code, GRID);
+    const aWin = getAnalysisWindow();
+    aWin.loadedGridVersions = aWin.loadedGridVersions || {};
+    aWin.loadedGridVersions[{{ node_id }}] = {{ node_version }};
 }
```

**Placeholder UI — reuse the three existing icon buttons.** The patient
phenotype toolbar (`patients/templates/patients/tags/phenotype_entry_tag.html:357-363`)
already has exactly this trio — a "Show Grid" icon plus CSV/VCF download icons —
backed by CSS in `global.css` (`.show-grid-icon` :2371, `.csv-icon` :2261,
`.vcf-icon` :2265). Reuse those classes rather than new buttons; CSV/VCF reuse
the node grid's existing `export_grid()` (which only needs `postData` from the
config, so it works with no rows loaded):

```html
{% if not grid_auto_load %}
<div id="grid-placeholder-{{ node_id }}" class="node-grid-placeholder">
    <p>{{ node.count|default:"?" }} variants — not loaded.</p>
    <a id="load-variants-{{ node_id }}" title="Show Grid" href="javascript:void(0)">
        <div class="show-grid-icon icon32"></div> Show grid
    </a>
    <a title="Download as CSV"
       href="javascript:export_grid({{ analysis_id }}, {{ node_id }}, '{{ node_id }}_{{ node_version }}', 'csv')">
        <div class="csv-icon icon32"></div> CSV
    </a>
    <a title="Download as VCF"
       href="javascript:export_grid({{ analysis_id }}, {{ node_id }}, '{{ node_id }}_{{ node_version }}', 'vcf')">
        <div class="vcf-icon icon32"></div> VCF
    </a>
</div>
{% endif %}
```

(The `.icon32 .has-phenotypes-icon hidden` combo in the phenotype toolbar is its
own show/hide logic; here we just want the static icon, so use the icon class +
`icon32` without `hidden`.)

### 3. `grid.js` — split config-load from data-fetch

`setupGrid` gains an `autoLoad` arg (default `true` to preserve other callers).
Initialise the grid from config either way; defer the row fetch via jqGrid's
`datatype: 'local'`, then flip to the server datatype and `reloadGrid` to fire
phase 2 on demand.

```diff
-function setupGrid(config_url, analysisId, nodeId, versionId, unique_code, gridComplete, gridLoadError, on_error_function) {
+function setupGrid(config_url, analysisId, nodeId, versionId, unique_code, gridComplete, gridLoadError, on_error_function, autoLoad) {
+    if (typeof autoLoad === "undefined") { autoLoad = true; }
     $(function () {
         $.getJSON(config_url, function(data) {
             ...
+            // Remember the server datatype so a deferred load can flip back to it.
+            window.nodeGridServerDatatype = window.nodeGridServerDatatype || {};
+            window.nodeGridServerDatatype[nodeId] = data["datatype"] || "json";
+            if (!autoLoad) {
+                data["datatype"] = "local";  // build colModel/pager/nav, fetch no rows
+            }
             const grid = getGrid(nodeId, unique_code);
             grid.jqGrid(data).navGrid(...);
             ...  // CSV / VCF / canonical-transcript nav buttons still added (need only postData)
         });
     });
 }
+
+// Fire the deferred row query (phase 2) for a node whose grid was config-loaded only.
+function loadNodeGridData(nodeId, unique_code) {
+    const grid = getGrid(nodeId, unique_code);
+    const datatype = (window.nodeGridServerDatatype || {})[nodeId] || "json";
+    grid.jqGrid('setGridParam', {datatype: datatype}).trigger('reloadGrid');
+}
```

Notes:
- `gridComplete` already fires after the data load and registers the GRID
  component / restores selected-variant checkboxes — that still runs on the
  deferred load. With `datatype: 'local'` and no local data, jqGrid does not run
  `gridComplete` for the empty config-only init, so the GRID component isn't
  registered until a real load. FilterNode's `connectFilterNodeEditorToGrid`
  uses `searchGrid` (colModel only) and does **not** depend on `gridComplete`,
  so it still wires up — confirm during implementation.
- `gridLoadError` already calls `hideLoadingOverlay()` (grid.js) and clears the
  container; unchanged.

### 4. Client-side stale-response guard (carried over from prior plan)

When the user clicks B while A's editor/data AJAX is in flight, A's late
response can clobber B. The server-side dedup covers DB cost; this is purely a
DOM concern. Drop stale responses by checking the unique-code marker still
exists:

```diff
 // base_node_data.html
-function load_node_editor(node_edit_url) {
+function load_node_editor(node_edit_url, unique_code) {
     $("#node-editor-container").empty();
     $.ajax(node_edit_url, {
-        success: function(html) { $("#node-editor-container").html(html); },
+        success: function(html) {
+            if ($("#" + unique_code, "#node-data-container").length === 0) { return; }
+            $("#node-editor-container").html(html);
+        },
         ...
```

```diff
 // node_data_grid.html — gridComplete()
 function gridComplete() {
     const unique_code = "{{ node_id }}_{{ node_version }}";
+    if ($("#" + unique_code, "#node-data-container").length === 0) {
+        return;  // user navigated away; ignore stale callback
+    }
     ...
     registerComponent(unique_code, GRID);
 }
```

## Diff summary by file

- `snpdb/models/models_user_settings.py`
  - Add `node_grid_auto_load_max_variants` to `SettingsOverride`; add the
    annotation to the `UserSettings` dataclass (cascade is automatic).
- `snpdb/migrations/` — new migration for the field.
- `snpdb/forms.py` (+ settings templates) — add the field to the Global/Org/Lab/
  User settings forms that list fields explicitly.
- `variantgrid/settings/components/default_settings.py`
  - Add `ANALYSIS_NODE_GRID_AUTO_LOAD_MAX_VARIANTS = 50_000` (cascade fallback).
- `analysis/views/views.py` (`node_data_grid` view)
  - Compute `grid_auto_load` from the cascaded setting + `node.count`; add to
    context.
- `analysis/templates/analysis/node_data/node_data_grid.html`
  - Branch docready on `grid_auto_load` (+ session `loadedGridVersions`); render
    placeholder with the reused Show-Grid / CSV / VCF icons; pass `autoLoad` to
    `setupGrid`; pass `unique_code` to `load_node_editor`; record the load and
    add the unique-code stale guard in `gridComplete`.
- `variantgrid/static_files/default_static/js/grid.js`
  - `setupGrid(..., autoLoad)`: init from config; `datatype:'local'` when
    deferred. Add `loadNodeGridData(nodeId, unique_code)` to fire phase 2.
- `variantgrid/static_files/default_static/css/global.css` (only if the
  placeholder needs layout polish — the icon classes already exist).
- `analysis/templates/analysis/node_data/base_node_data.html`
  - `load_node_editor(url, unique_code)`: drop stale response.

One DB migration (the cascading-settings field). No change to `NodeGridHandler` /
`NodeGridConfig` / `node_grid_export`.

## Edge cases / regressions to test

- **Small node** (`count < threshold`): loads rows automatically, identical to
  today.
- **Large node** (`count >= threshold`): placeholder shows the count; editor is
  fully interactive immediately; no row query runs until "Load variants".
- **Load variants click**: fires the deferred query once; pager, sorting,
  paging, selected-variant checkboxes all work afterwards.
- **FilterNode** (any size): the `searchGrid` filter UI appears from colModel
  with no row fetch; saving a filter still `reloadNodeAndData`s.
- **CSV / VCF / canonical-transcript download** from the placeholder (no rows
  loaded): `export_grid` builds the URL from `postData` and downloads. Verify
  the export iterator runs the query server-side independent of the browser grid
  state.
- **`count is None`** (node never computed, or just edited → dirty): treated as
  "large" → placeholder. After an edit + node-update cycle, `count` repopulates.
- **Tag / selection node** (`ANALYSIS_TAGS_NODE_ID` via `viewTags()`): threshold
  applies uniformly; in practice always under it (nobody tags 50k variants), so
  it auto-loads as today.
- **Session re-display**: load a large node (placeholder → click), switch to
  another node, switch back → grid re-shows automatically (no placeholder),
  served from the page cache. Edit the node → version bumps → placeholder again.
- **Click A then B while A in flight**: only B's editor/grid render (stale
  guard).
- **Dual-screen mode** (`secondWindow`): placeholder + Load button live in the
  grid window; verify the button handler resolves the right grid.
- **Read-only analysis** (`set_form_read_only`): editor fields disabled; Load
  variants + CSV still work (they're reads).
- **`extra_filters`** (e.g. show damage predictions) re-enters
  `node_data_grid`: re-evaluates `grid_auto_load`; same behaviour.
- **Node errors path** (`node_errors.html`): no grid; unaffected.

## Resolved decisions

1. **Threshold** = cascading setting (Global → Org → Lab → User), default
   `50_000`, reusing the `SettingsOverride` mechanism.
2. **Tag/selection nodes** — threshold applies uniformly (no special-case).
3. **Placeholder UI** — reuse the existing Show-Grid / CSV / VCF icon buttons
   from the patient phenotype toolbar (`.show-grid-icon` / `.csv-icon` /
   `.vcf-icon`).
4. **Session re-display** — if a node-version was already loaded this session,
   auto-load it again (page-cache hit, instant); editing bumps the version so it
   re-defers.

## Tests

Worth adding (the rest is interactive JS, system-tested by hand):

- **`grid_auto_load` computation** (`analysis` view/unit test): a node with
  `count >= threshold` → `grid_auto_load` False and the rendered
  `node_data_grid` contains the placeholder (`load-variants-<id>`) and no
  auto data fetch; a node with `count < threshold` → True, no placeholder. Use a
  fake node/analysis from the existing analysis test helpers; assert on the
  template context and/or rendered HTML.
- **Cascading setting resolution** (`snpdb` settings test): set the field at
  Global vs Lab vs User levels and assert
  `UserSettings.get_for_user(user).node_grid_auto_load_max_variants` returns the
  most-specific non-null value, and falls back to
  `settings.ANALYSIS_NODE_GRID_AUTO_LOAD_MAX_VARIANTS` when all are null.
