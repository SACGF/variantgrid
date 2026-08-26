# Analysis layout rework: horizontal mode with Swap + Drawer editors

**Goal:** a per-analysis "horizontal mode": nodes flow left→right on a full-width canvas, the
variant grid docks along the bottom at full width, and the node editor is available two ways —
docked as a bottom-pane tab beside the grid ("Swap") or sliding in over the canvas ("Drawer") —
with a live toggle between them so both can be used and compared. Mockups and rationale:
`claude/design/analysis_grid_rework.html`.

The variant grid stays jqGrid throughout — the DataTables conversion is a separate plan.
The existing vertical layout keeps working unchanged and remains the default; everything new is
gated on `analysis.analysis_horizontal_mode`.

**Prerequisite:** `claude/plans/jsplumb_upgrade_plan.md` has landed (anchor work below is written
once, against the new API).

## 1. Settings wiring (fields already exist)

- `Analysis.analysis_horizontal_mode` (`analysis/models/models_analysis.py:65`): add to
  `AnalysisForm.Meta.fields` (`analysis/forms/forms.py`) so it appears in analysis settings.
- `UserSettings.analysis_horizontal_mode` (`snpdb/models/models_user_settings.py:216`): expose the
  form field by making `get_settings_form_features().analysis_horizontal_mode` return
  `self.analysis_enabled` (`snpdb/forms.py:483` — currently hard-coded off) and pointing the
  `_hide_unused_fields` entry (`snpdb/forms.py:606`) at that feature flag.
- New analyses inherit the user setting in `set_defaults_and_save`
  (`analysis/models/models_analysis.py:~222`, alongside `custom_columns_collection` etc.).
  Analyses created from a template take the **template's** mode — node positions copy from the
  template, and positions are orientation-specific.
- **Mode switch transposes positions:** when `AnalysisForm` save changes the mode, swap x/y on all
  the analysis's nodes in the same transaction. Swap is self-inverse, so toggling back restores the
  original layout exactly. No version bump — positions don't affect queries.

## 2. Page layout (templates + Split.js)

- `view_analysis` (`analysis/views/views.py:198`) passes the mode; `analysis.html` renders the
  horizontal structure when set: toolbar, then `#canvas-row` (full width), then a horizontal
  gutter, then `#bottom-panel` (full width) holding `analysis_editor_and_grid.html`.
- `layoutAnalysisPanels` (`analysis.js:88`) grows a horizontal variant:
  `Split(['#canvas-row', '#bottom-panel'], {direction: 'vertical', ...})`.
  `analysis_panel_fraction` keeps its meaning as "fraction taken by the node-canvas side" in both
  modes, and `analysis_set_panel_size` (`analysis/views/views_json.py:431`) is reused as-is.
- `resizeGrid` (`analysis.js:32`) sizes to the bottom panel (full window width) and also sets the
  grid height to fill the pane; call it on splitter drag-end and on grid-tab show (a jqGrid sized
  while hidden measures zero width).
- The template/analysis-variables top bar and its vertical split stay exactly as they are.

## 3. Swap: tabbed bottom pane

`analysis_editor_and_grid.html` gains a tab strip — **Editor | Grid** — over the existing
`#node-editor-container` and `#node-grid-container`. Container ids and the
`registerComponent(EDITOR/GRID)` pairing stay untouched; the tabs only control visibility.

- **Tab behaviour:** the active tab is sticky while browsing — clicking nodes with the Grid tab up
  loads each node's grid (results-browsing mode); clicking nodes with the Editor tab up shows each
  node's editor. Saving an editor switches to the Grid tab, which shows the existing async-wait /
  loading states until the node re-computes.
- **Toolbar-driven editor content** (`analysisSettings()`, `inputSamples()`, `viewTags()` via
  `replaceEditorWindow`) switches to the Editor tab when it loads.
- **Lazy grid:** while the Grid tab is hidden, defer the grid data fetch until the tab is first
  shown — integrate with the existing deferred machinery in
  `analysis/templates/analysis/node_data/node_data_grid.html` (`grid_auto_load` placeholder,
  `loadedGridVersions` page-cache, and the `datatype:'local'` empty init that lets the editor's
  `everythingLoaded` proceed without a row fetch). Editing a chain of nodes with the Editor tab up
  therefore issues zero grid queries; the `grid_auto_load_max_variants` placeholder then applies as
  normal when the tab is shown.

## 4. Drawer: the same editor, undocked

- A drawer element on the right edge of `#canvas-row`, hidden by default, ~450px wide, resizable by
  dragging its left edge, closable via a button and Escape.
- **One editor, two homes:** toggling reparents the existing `#node-editor-container` between the
  Editor tab pane and the drawer body — all load logic targets that id and keeps working. In drawer
  mode the bottom panel is grid-only (tab strip hidden) and node selection opens the drawer; the
  grid stays visible below, so after a save the re-computed results appear without any tab flip.
- **Toggle placement:** an "undock" button on the Editor tab strip, and a "dock as tab" button in
  the drawer header. Persist the choice in `localStorage` for the comparison period; once a winner
  (or a lasting preference) emerges, promote it to a `UserSettings` field in a follow-up change.

## 5. Node canvas: horizontal flow

- Centralise orientation in one config object consumed by `setupConnections` and the instance
  defaults in `analysis_nodes.js` — vertical: input `TopCenter` (Venn: `[0.25,0]`/`[0.75,0]`),
  output `BottomCenter`; horizontal: input `LeftMiddle` (Venn: left edge at 0.25/0.75 height),
  output `RightMiddle`. Everything else about endpoints (uuids, styles, MergeNode multi-input,
  invalid-target highlighting) is orientation-independent and stays shared.
- New-node placement: horizontal mode places an added/copied node to the right of the active node
  (vertical mode keeps placing below).
- Node cards, chips, badges and the counts strip are unchanged — the card design works in both
  orientations.

## 6. Interactions to keep working in both modes

Dual-screen mode, Delete-key node deletion, marquee multi-select + multi-drag, node count
click-through (`extra_filters`), grid export CSV/VCF (including from the placeholder), the
tags view, and read-only analyses (no editing, no dragging, editor read-only).

## Verification (manual, `python3 manage.py runserver`)

1. A vertical analysis renders and behaves exactly as before (regression pass).
2. Toggle an analysis to horizontal in analysis settings: nodes transpose sensibly; toggle back
   restores the original layout exactly.
3. New analysis inherits the user setting; analysis-from-template inherits the template's mode.
4. Swap: sticky tabs while clicking around the graph; save flips to Grid with the loading state;
   with the Editor tab up, no grid queries fire (check the network tab); Grid tab first-show loads
   or shows the auto-load placeholder per the row-count threshold.
5. Drawer: undock/redock round-trips, preference survives reload, Escape closes, drawer-edited
   save shows refreshed grid below.
6. Splitter drag persists `analysis_panel_fraction`; grid resizes on drag-end and tab show.
7. Items in §6 all pass in horizontal mode.
