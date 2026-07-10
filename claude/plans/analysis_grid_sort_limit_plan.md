# Plan: Disable analysis grid sorting on large nodes

## Problem

Analysis node grids (`/analysis/<id>/node_grid/handler/`) let the user sort by any column,
including deeply-joined, unindexed text columns such as
`variantannotation__gene__geneannotation__omim_terms` (a plain `TextField` on `GeneAnnotation`,
reached via `variant → variantannotation → gene → geneannotation`).

Sorting on such a column forces Postgres to fully sort the entire joined result set before it can
return even the first page (the `LIMIT` applies after the sort, and no index can satisfy the
ordering across the join). On large nodes this runs past the 300s `major_operation`
`statement_timeout` and Postgres kills the query.

### Evidence (from the production event log, July 2026)

- The failing requests die with `psycopg2.errors.QueryCanceled: canceling statement due to
  statement timeout` — genuine 300s DB timeouts, not the client or the concurrency cap.
- **210** `node_grid` statement-timeouts are logged across dozens of analyses; it is a long-running
  systemic issue, not a one-off.
- **Zero** `major_operation` concurrency-cap events exist — contention was not the cause.
- Sampled node sizes (`node.count`) for the reported URLs: the timeouts line up with the large
  nodes (e.g. 4,697,584 rows and 23,607 rows sorted by `omim_terms`), while the small nodes in the
  same analyses (23–412 rows) sort in milliseconds and are not the ones dying.

Medical scientists reach for `omim_terms` sorting to surface variants that have OMIM records. The
correct workflow is to **filter** the grid down first; sorting a multi-million-row grid is the
wrong tool. This plan removes the footgun and points users at filtering.

## Behaviour change

When a node grid has **`count >= ANALYSIS_GRID_SORT_MAX_ROWS`** (default 10,000), or its count is
unknown (`None`):

- All columns become non-sortable in the grid UI (no sort arrows, no `sidx` sent).
- The server ignores any incoming `sidx`/`sortname` and orders the page by `-pk` only — the variant
  primary key, an index scan + `LIMIT`, which is fast and gives stable pagination.
- A short banner tells the user why, and that filtering under the limit re-enables sorting.

Below the limit, sorting behaves exactly as today.

This intentionally disables *all* sorting above the limit (including otherwise-cheap columns like
position), trading that for simplicity. Per-column exceptions can be added later if a specific cheap
column proves worth re-enabling.

## Changes

### 1. Setting

`variantgrid/settings/components/default_settings.py`, alongside the `MAJOR_OPERATION_*` block
(~line 107):

```python
# Analysis node grids sort in-DB via ORDER BY. Sorting by joined/unindexed columns (e.g. OMIM) on a
# large result set forces a full sort that blows the statement_timeout. Above this row count we
# disable sorting entirely and fall back to ORDER BY -pk (indexed). Users filter down to re-enable.
ANALYSIS_GRID_SORT_MAX_ROWS = 10_000
```

### 2. Server-side gate — `analysis/grids.py`

Add a helper on `VariantGrid` that reports the current view's row count (mirroring `get_count`'s
source of truth) and whether sorting is allowed:

```python
def _grid_row_count(self) -> Optional[int]:
    if self.node_count:
        return self.node_count.count
    return self.node.count

def sorting_disabled(self) -> bool:
    count = self._grid_row_count()
    return count is None or count >= settings.ANALYSIS_GRID_SORT_MAX_ROWS
```

In `VariantGrid._sort_items` (grids.py:281), short-circuit at the top so the packed-field branch and
the incoming `sidx` are skipped when over the limit — `super()._sort_items` with `sidx=None` already
appends `-pk`:

```python
def _sort_items(self, items, sidx, sord):
    if self.sorting_disabled():
        return super()._sort_items(items, None, sord)  # ORDER BY -pk only
    ...  # existing CohortGenotype packed-field handling
```

When over the limit, also leave `extra_config['sortname']` unset (grids.py:63-65) so the initial
grid load does not request a sort on `default_sort_by_column`.

### 3. Colmodels — `analysis/grids.py`

In `VariantGrid.get_colmodels` (grids.py:101), when `self.sorting_disabled()`, set
`sortable: False` on every returned colmodel so the header sort arrows disappear and the client
never sends a `sidx`.

### 4. User instruction — `analysis/templates/analysis/node_data/node_data_grid.html`

The template already shows `node.count` in the placeholder (line 117). Pass a context flag from the
node grid data view (the view rendering this template) — `grid_sorting_disabled` and
`grid_sort_max_variants` (`settings.ANALYSIS_GRID_SORT_MAX_ROWS`) — and render an info banner above
the grid table (around line 133) shown only when `grid_sorting_disabled`:

> ⓘ This grid has {{ node.count|intcomma }} variants. Sorting is disabled above
> {{ grid_sort_max_variants|intcomma }} rows because it is too slow. **Add a filter to narrow your
> results, then sorting re-enables.** To find variants with OMIM records, filter for them rather
> than sorting.

Locate the view supplying this template's context (renders `node_data_grid.html`) and add the two
context variables there.

## Verification

- `python3 -m py_compile` the edited Python files.
- Add a unit test (run in a non-production environment) covering `VariantGrid._sort_items`:
  a node with `count >= ANALYSIS_GRID_SORT_MAX_ROWS` produces `order_by("-pk")` regardless of the
  requested `sidx`; a node below the limit keeps the requested ordering.
- Add a test that `get_colmodels` marks all columns `sortable: False` when over the limit.

## Out of scope (deliberately deferred)

- Per-column "cheap vs expensive" sort classification, index introspection, and a `sortable` flag on
  `VariantGridColumn`.
- Sweeping/clearing saved `default_sort_by_column` values or notifying users — unnecessary, because
  the gate is evaluated live at render time and a stored sort column is simply ignored above the
  limit.
- Adding new filter nodes/columns (e.g. a broad "has OMIM / gene" filter). Filtering already exists;
  the banner steers users toward it.
