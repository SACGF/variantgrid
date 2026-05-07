# Always-snappy node editor (issue #51)

Goal: when the user clicks a node, show the editor form **immediately** so they
can change settings (cohort, sample, filters, etc.) without waiting for the
grid to load. The grid loads independently underneath.

This is a UI-layer change; the server-side scheduling/gating that already
exists is preserved.

## Today's flow

```
click node
  └─> loadNodeData(nodeId) ─ analysis.js:27
        └─> loadGridAndEditorForNode(nodeId) ─ analysis.js:459
              ├─ dataContainer.attr("node_url", load_node_url)  ─── new URL stamp
              ├─ #node-editor-container.empty()                  ─── editor cleared
              ├─ showLoadingOverlay()  ─ analysis.js:515         ─── *** overlay covers whole right-panel ***
              └─ dataContainer.load(Urls.node_load(...))         ─── AJAX

server: node_load (views.py:675)
  ├─ NodeStatus.is_error  → renders node_errors.html
  ├─ NodeStatus.is_ready  → renders node_data_grid.html  (grid <table> + script)
  └─ otherwise            → renders node_async_wait.html (polls via messagePoller, re-fires
                                                          loadNodeData when status flips to ready)

node_data_grid.html docready:
  ├─ load_node_editor(node_view_url)   ── AJAX into #node-editor-container
  └─ setupGrid(...)                    ── jqgrid AJAX

editor finishes:  finishedLoadingEditor() → registerComponent(unique_code, EDITOR, everythingLoaded)
grid finishes:    gridComplete()        → registerComponent(unique_code, GRID)
                                          ↓
                  registerComponent (analysis_editor_and_grid.html:34) requires BOTH
                  EDITOR + GRID, then runs all queued callbacks
                                          ↓
                  everythingLoaded()  → hideLoadingOverlay()
```

### Why the editor "feels slow"

- `showLoadingOverlay()` appends `#overlay-container` to `#right-panel`, which
  contains both `#node-editor-container` and `#node-grid-container`.
- The overlay isn't dismissed until `everythingLoaded`, which fires only after
  *both* EDITOR and GRID have registered.
- The editor HTML usually comes back much faster than the grid (the grid view
  triggers Celery-driven node updates and waits on them via `node_async_wait`),
  so the editor is sitting fully rendered behind the overlay.

### Server-side "gate" (the dedup story you remembered)

Three layers cooperate so duplicate clicks don't run duplicate work:

1. **Celery task scheduling** — `NodeTask` (analysis_node.py:1110) has
   `unique_together = (node, version)`. `analysis_update_tasks.py:71` does a
   `bulk_create(..., ignore_conflicts=True)`, so a second scheduling attempt
   for the same (node, version) is silently dropped.

2. **`node_load`** itself is `@never_cache` and idempotent: it inspects current
   status and redirects to one of three child templates. `node_async_wait` then
   registers a one-shot `messagePoller.observe_node(id, "ready")`; when the
   node flips to ready it calls `loadNodeData(id)` again, which the server
   now redirects to `node_data_grid`. The poll callback also guards re-entry:
   only re-fires if `node_id == getLoadedNodeId()`.

3. **`NodeGridHandler.get`** (the heavy grid-data endpoint) — *this* is the
   per-request gate the user remembers:

```python
@method_decorator([cache_page(WEEK_SECS), vary_on_cookie], name='get')
class NodeGridHandler(NodeJSONViewMixin):
    def get(self, request, *args, **kwargs):
        """ This can be a really expensive operation (ie a few mins)
            And users can sometimes click multiple times, causing the DB to get
            slow running the same query multiple times, interfering with itself
            — so make a per-user lock, and redirect any further calls which
            should hopefully hit the cache next time
        """
        LOCK_EXPIRE = 60 * 10
        node = self._get_node(request)
        url = reverse("node_grid_handler", kwargs={"analysis_id": node.analysis_id})
        url = _add_allowed_node_grid_params(url, request.GET.dict())
        lock_id = sha256sum_str(f"{url}_{request.user}")
        if cache.add(lock_id, "true", LOCK_EXPIRE):
            try:
                response = self.get_response(request, *args, **kwargs)
            finally:
                cache.delete(lock_id)
        else:
            time.sleep(2)
            response = HttpResponseRedirect(url)
        return response
```

Per-(user, URL+params) `cache.add` lock. The first caller computes; concurrent
callers sleep 2s then redirect to the same URL, which is `cache_page(WEEK_SECS)`
+ `vary_on_cookie` — so by then the response is in the page cache and the
redirect resolves instantly without re-running the query. Hence: the grid AJAX
"only runs once" even if you spam clicks; further clicks land on the cache.

So "make one request, redirect when ready" is really:
**Celery dedup at the scheduling layer + status polling for the page render +
per-user lock-then-cache for the grid AJAX**. Cancelling celery tasks for
nodes the user clicked away from is unnecessary — the existing task finishes,
its result populates `cache_page`, and the next visit is free.

## Proposed change

Single principle: **the editor is its own loading lane, independent of the grid.**

### 1. Scope the loading overlay to the grid area only

`analysis.js`, `showLoadingOverlay`:

```diff
 function showLoadingOverlay() {
     const oc = $("#overlay-container");
     if (!oc.is(":visible")) {
-        // Move to right-panel (with top z-order), then things can load underneath.
-        oc.appendTo("#right-panel");
+        // Cover only the grid — editor renders independently and is
+        // interactive while the grid is still loading.
+        oc.appendTo("#node-grid-container");
         oc.show();
         …
```

`#node-grid-container` already exists (analysis_editor_and_grid.html:94) and
sits to the right/below of `#node-editor-container`, so positioning the
overlay's CSS (`top:0; bottom:0; left:0; right:0;`) within it covers the right
panel grid only.

### 2. Drop "wait for grid" from editor activation

Currently `everythingLoaded` (analysis.js:546) is what dismisses the overlay,
and it requires both EDITOR and GRID. We change the *meaning* of the registry:

- Editor render is no longer gated on the grid.
- The registry is kept **only** for components that genuinely need cross-wiring
  (FilterNode connecting its editor to the grid).

```diff
 function finishedLoadingEditor(node_id, version_id) {
-    const everythingLoaded = function () {
-        hideLoadingOverlay();
-    };
-
     const unique_code = node_id + "_" + version_id; // make sure only attach editor to grid that requested
-    registerComponent(unique_code, EDITOR, everythingLoaded);
+    registerComponent(unique_code, EDITOR);
 }
```

And the overlay is dismissed when the **grid** finishes loading:

`node_data_grid.html`:

```diff
 function gridComplete() {
     …
     const unique_code = "{{ node_id }}_{{ node_version }}";
-    registerComponent(unique_code, GRID);
+    registerComponent(unique_code, GRID);
+    hideLoadingOverlay();
 }
```

`grid.js` also calls `gridLoadError` on failure — push `hideLoadingOverlay()`
into that callback too, so a grid error doesn't leave a permanent overlay.

### 3. Keep registry for FilterNode wiring

`filternode_editor.html:64` uses `registerComponent(EDITOR, connectFilterNodeEditorToGrid)`
to defer DOM wiring until the grid table exists. That continues to work as-is
because `registerComponent` still fires queued callbacks when both components
register.

We do need to make sure the FilterNode editor body shows a placeholder while
waiting for the grid, so the user sees something:

```diff
 <div id="FilterNode-editor"></div>
+<div id="FilterNode-editor-loading">Loading filters…</div>
 <div id="filter-messages"></div>
```

And inside `connectFilterNodeEditorToGrid`:

```diff
 if (grid && grid.length > 0) {
+    $("#FilterNode-editor-loading").hide();
     // Make search always appear
     …
 }
```

The rest of the FilterNode editor (group operation toggle, save state) is
already interactive immediately because it's part of the editor template —
the only thing waiting is the search dialog injection.

### 4. Form-save + dirty-state interaction

When the user changes settings and submits the editor form:

1. `ajaxForm` POSTs to `node_view` (UpdateView).
2. `NodeView.form_valid` (node_view.py:65) sets `queryset_dirty = True` and
   calls `update_analysis(...)` which schedules a new `NodeTask`.
3. On success, the JS calls `reloadNodeAndData(node.pk)` which calls
   `loadNodeData(node.pk)` → `loadGridAndEditorForNode(node.pk)` again.
4. Server returns `node_async_wait` (status is now LOADING) → poll until ready
   → re-fire load.

If the user submits while a previous grid load is still in flight, the
re-entry into `loadGridAndEditorForNode` already empties the editor container
and replaces the data container, so old in-flight responses populating into
the new container is the only race we need to worry about.

### 5. Stale in-flight responses on node switch

Given the server-side gate above (per-user lock + `cache_page`), spamming
clicks doesn't multiply DB work — the second click sleeps then redirects into
the cached response. The remaining concern is purely *client-side*: when the
user clicks B while A's grid AJAX is still in flight, A's response will
arrive after B's containers have been swapped in.

What actually happens today:

- `load_node_editor` (`base_node_data.html:47`) injects HTML into
  `#node-editor-container`. If A's response arrives *after* B has loaded, it
  overwrites B's editor. ⚠️
- `setupGrid` populates `#grid-{node_id}` *inside* `#{node_id}_{node_version}`
  (a uniquely-scoped subtree). When A's data container subtree is replaced
  by B's, the orphan `<table>` is gone from the DOM — jqgrid operating on a
  detached element is harmless (no visible UI effect).
- `gridComplete` calls `registerComponent(unique_code, GRID)` and (after this
  change) `hideLoadingOverlay()`. If A fires after B started, A's callback
  could prematurely hide the overlay for B. ⚠️

So the correctness fix is small and scoped:

```diff
 // base_node_data.html
-function load_node_editor(node_edit_url) {
-    $("#node-editor-container").empty();
-    $.ajax(node_edit_url, {
-        success: function(html) { $("#node-editor-container").html(html); },
+function load_node_editor(node_edit_url, unique_code) {
+    $("#node-editor-container").empty();
+    $.ajax(node_edit_url, {
+        success: function(html) {
+            // Drop response if user navigated to a different node.
+            if ($("#" + unique_code, "#node-data-container").length === 0) {
+                return;
+            }
+            $("#node-editor-container").html(html);
+        },
         …
```

```diff
 // node_data_grid.html — gridComplete()
 function gridComplete() {
     const unique_code = "{{ node_id }}_{{ node_version }}";
+    if ($("#" + unique_code, "#node-data-container").length === 0) {
+        return;  // user navigated away; ignore stale callback
+    }
     …
     registerComponent(unique_code, GRID);
+    hideLoadingOverlay();
 }
```

Both checks rely on the unique-code DOM marker that `node_data_grid.html`
already emits (`<div id="{{ node_id }}_{{ node_version }}">…</div>`).

## Diff summary by file

- `variantgrid/static_files/default_static/js/analysis.js`
  - `showLoadingOverlay()` → append to `#node-grid-container` instead of
    `#right-panel`.
  - `finishedLoadingEditor()` → drop `everythingLoaded` gating of overlay.
- `analysis/templates/analysis/node_data/node_data_grid.html`
  - `gridComplete()` → call `hideLoadingOverlay()` directly; pass `unique_code`
    to `load_node_editor`.
- `analysis/templates/analysis/node_data/base_node_data.html`
  - `load_node_editor(url, unique_code)` → drop response if user navigated.
  - In `gridLoadError` path (defined in calling templates), also call
    `hideLoadingOverlay()`.
- `analysis/templates/analysis/node_editors/filternode_editor.html`
  - Add `Loading filters…` placeholder; hide it inside
    `connectFilterNodeEditorToGrid`.

No server-side changes.

## Edge cases / regressions to test

- Click slow node → editor visible and interactive within ~100 ms.
- Click slow node A → click slow node B while A loading → only B's editor +
  grid render; no flash of A.
- Click node, change cohort, submit → `reloadNodeAndData` triggers new wait
  cycle; overlay reappears over grid only.
- FilterNode click: editor headers/group-op visible immediately; search
  dialog appears once grid loads; "Loading filters…" placeholder hidden then.
- `cancelNodeLoad` button (lives in `#overlay-container`): still appears with
  the overlay over the grid area; clicking still revokes the celery task.
- Dual-screen mode (`secondWindow`): overlay lives in the grid+editor window;
  `#node-grid-container` exists there too — verify path.
- `node_errors` view (`base_node_data.html` parent): on error there is no
  grid; `on_error_function` already calls `hideLoadingOverlay()` (line 10 of
  base_node_data.html). Confirm the editor is still rendered (it is — error
  path goes through `load_node_editor`).
- Tag node (`ANALYSIS_TAGS_NODE_ID` via `viewTags()`): same path as a normal
  node load; nothing special.
- Read-only analysis (`set_form_read_only` in NodeView.get_form): editor
  fields disabled; nothing about loading changes.

## Open questions

1. **FilterNode body while grid loads.** The proposal shows a "Loading
   filters…" placeholder. Alternative: render the search dialog skeleton
   server-side from the colmodels we already have and have it be inert
   until the grid AJAX returns. Worth the extra code, or is the placeholder
   fine?
2. **Cancel-load button placement.** With the overlay scoped to
   `#node-grid-container`, the cancel button is anchored over the grid area.
   Acceptable, or do you want it pinned to a fixed corner of the right panel
   so it's reachable while the user is editing the form?
3. **Stale-response gating.** §5 is now small (two `unique_code` checks)
   given the server-side dedup already covers the expensive part. Land it
   alongside §1–4, or only if observed?
4. **Anything else relying on `everythingLoaded`?** I only see overlay
   dismissal. Worth a final grep before I implement — e.g. confirm no plugin
   or external skill hooks the registry's "both" callback path.
5. **Test coverage.** This is JS UI behaviour — no Python tests would catch a
   regression. Do you want a manual regression checklist baked into the PR
   description, or would you accept a lightweight Selenium/Playwright smoke
   test for "editor visible within X ms after click"?
