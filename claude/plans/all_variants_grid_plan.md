# All Variants Grid Timeout — Fix Plan

**Issue:** SACGF/variantgrid#1279

## Problem Summary

The All Variants Grid (`AllVariantsGrid` in `variantopedia/grids.py`) times out at the PostgreSQL statement timeout. The variant table can reach hundreds of millions of rows; `VariantZygosityCount` is the same scale (one row per variant per collection).

There are two distinct timeout paths, both hit on every grid load:

1. **COUNT timeout** — `get_data` (`library/jqgrid/jqgrid.py:346`) accesses `paginator.num_pages` and `paginator.count`, both of which hit `COUNT(*)` over the full filtered+joined dataset.
2. **Sort timeout** — `get_items` runs `ORDER BY (het_count + hom_count) DESC LIMIT 25` over the join result. This is a computed expression on a joined table; no index can apply, so PostgreSQL must sort the entire join result to find the top 25 rows.

---

## Root Cause Analysis

### Double join (diagnosed previously)

`get_annotation_kwargs` (`snpdb/models/models_zygosity_counts.py:52`) annotates with a `FilteredRelation` → LEFT OUTER JOIN (alias `global_variant_zygosity`). Then `AllVariantsGrid._get_q()` (`variantopedia/grids.py:82`) adds:

```python
# benchmarking - I found it much faster to do both of these queries (seems redundant)
hom_nonzero = Q(**{f"{self.vzcc.hom_alias}__gt": 0})   # Q(global_variant_zygosity__hom_count__gt=0)
het_nonzero = Q(**{f"{self.vzcc.het_alias}__gt": 0})    # Q(global_variant_zygosity__het_count__gt=0)
filter_list.append(hom_nonzero | het_nonzero)
filter_list.append(Q(**{f"{self.vzcc.non_ref_call_alias}__gt": 0}))
```

Filtering on a FilteredRelation path (`global_variant_zygosity__hom_count__gt`) forces Django to create a **new** INNER JOIN (T7) rather than reusing the existing LEFT OUTER JOIN. The result is two joins to the same large table. The `hom | het` condition is logically redundant with `non_ref_call > 0` (since `non_ref_call = het + hom`).

The resulting SQL:

```sql
SELECT COUNT(*)
FROM snpdb_variant
LEFT OUTER JOIN snpdb_variantzygositycount global_variant_zygosity
    ON (snpdb_variant.id = global_variant_zygosity.variant_id AND global_variant_zygosity.collection_id = 1)
INNER JOIN snpdb_variantzygositycount T7
    ON (snpdb_variant.id = T7.variant_id AND T7.collection_id = 1)   -- redundant
WHERE (T7.hom_count > 0 OR T7.het_count > 0)
  AND (global_variant_zygosity.het_count + global_variant_zygosity.hom_count) > 0
  AND [contig filter]
ORDER BY (global_variant_zygosity.het_count + global_variant_zygosity.hom_count) DESC NULLS LAST
LIMIT 25
```

### Why fixing the double join alone is not enough at scale

Even with one join, at hundreds of millions of rows:

- `COUNT(*)` still requires scanning and joining the entire filtered dataset — there is no shortcut.
- `ORDER BY (het_count + hom_count) DESC` is a computed expression on a joined column. PostgreSQL cannot use any index for this; it must materialise and sort the full join result before returning page 1.

The benchmarking result that motivated the double join was likely measured on a much smaller dataset.

---

## Proposed Fix — Three Parts

### ~~Part 1 — Remove the double join~~ (deferred — see note)

The double join was added based on benchmarking ("I found it much faster to do both of these queries"). Parts 2 and 3 already fix both timeout paths (COUNT and sort), so removing the double join now would risk a performance regression without measured benefit. Defer until Option C (stored `non_ref_count` column) is implemented, at which point the computed-expression filter is replaced entirely and the double join becomes moot.

---

### Part 2 — Replace precise COUNT with an approximate estimate, displayed honestly

The `COUNT(*)` is only used for jQGrid's pagination display (`records` and `total` fields in the JSON response). An approximate count is perfectly acceptable here — users will never manually page through millions of records. But the UI must not present a planner estimate as if it were an exact figure.

#### 2a — Server: EXPLAIN-based estimate + formatted label

Override `get_count` in `AllVariantsGrid` to use PostgreSQL's query planner estimate, and override `get_data` to add an `approximate_records` label to the JSON response.

**File:** `variantopedia/grids.py`

```python
def get_count(self, request) -> int:
    import re
    from django.db import connection
    qs = self.get_queryset(request)
    sql, params = qs.query.sql_with_params()
    with connection.cursor() as cursor:
        cursor.execute(f"EXPLAIN {sql}", params)
        first_line = cursor.fetchone()[0]
    match = re.search(r'rows=(\d+)', first_line)
    return int(match.group(1)) if match else 10_000_000

def get_data(self, request) -> dict:
    data = super().get_data(request)
    data['approximate_records'] = _format_approx_count(data['records'])
    return data
```

Add a module-level helper in `variantopedia/grids.py`:

```python
def _format_approx_count(n: int) -> str:
    """Format a large approximate count as '~100M', '~1.2B', etc."""
    for threshold, suffix in ((1_000_000_000, 'B'), (1_000_000, 'M'), (1_000, 'K')):
        if n >= threshold:
            rounded = n / threshold
            fmt = f"{rounded:.0f}" if rounded >= 10 else f"{rounded:.1f}"
            return f"~{fmt}{suffix}"
    return f"~{n}"
```

`get_count` returns the raw integer (used by the paginator for page-count math). `approximate_records` carries the display string — the two are kept separate so pagination still works correctly.

Note: `get_count` in `AbstractVariantGrid` (`snpdb/grids.py:489`) calls `self._get_queryset(request)`. The override above calls `self.get_queryset(request)` instead — verify that both resolve to the same queryset, or use whichever is consistent with `get_items`.

#### 2b — Client: display the approximate label in the jQGrid pager

jQGrid's `loadComplete` callback receives the full JSON response (including any extra fields). Add handling in `snpdb/templates/jqgrid/jqgrid.html` so that when `approximate_records` is present, the pager's record count is replaced with the formatted label.

jQGrid renders the record count in a `.ui-paging-info` element with text like `"1 - 10 of 102012921"`. Replace the trailing number with the approximate label:

**File:** `snpdb/templates/jqgrid/jqgrid.html`, inside `setupJQGrid{{ unique_code }}`, after `data["pager"] = pagerId`:

```javascript
const existingLoadComplete = data.loadComplete;
data.loadComplete = function(responseData) {
    if (responseData.approximate_records) {
        $(pagerId).find('.ui-paging-info').each(function () {
            // Replace trailing integer (the records count) with the approximate label
            $(this).text($(this).text().replace(/\d[\d,]*\s*$/, responseData.approximate_records));
        });
    }
    if (existingLoadComplete) {
        existingLoadComplete.call(this, responseData);
    }
};
```

This means the pagination display shows `"1 - 10 of ~100M"` rather than `"1 - 10 of 102012921"`, making clear the figure is an estimate. The integer `records` field is still used for page-count math; only the rendered text changes.

The `loadComplete` wrapper pattern chains correctly with any existing `loadComplete` set by the `{% if grid_complete %}` block.

---

### Part 3 — Change default sort to `pk`

The default sort `sortname: non_ref_call_counts_N` forces sorting a hundreds-of-millions-row join result by a computed expression on every page load. Change the default sort to `pk` (the Variant primary key, a B-tree index).

**File:** `variantopedia/grids.py`, `AllVariantsGrid.__init__`

**Before** (line 75):
```python
self.extra_config.update({'sortname': self.vzcc.non_ref_call_alias,
                          'sortorder': "desc",
                          'shrinkToFit': False})
```

**After:**
```python
self.extra_config.update({'sortname': 'id',
                          'sortorder': "desc",
                          'shrinkToFit': False})
```

With `ORDER BY id DESC LIMIT 25`, PostgreSQL can:
1. Scan the Variant table via the PK index (cheapest possible)
2. Probe `VariantZygosityCount` via nested-loop for each candidate
3. Stop after finding 25 rows passing the `non_ref_count > 0` filter

This is orders of magnitude faster than a full sort of the join result.

Users can still click any column header to sort by that column — the change only affects the initial grid load. The `non_ref_call` column remains sortable; it will just be slow to sort by it until Option C (stored column) is addressed.

---


## Files to Change

| File | Change |
|------|--------|
| `variantopedia/grids.py` | Part 2a: override `get_count` + `get_data` with EXPLAIN estimate and `approximate_records` label |
| `snpdb/templates/jqgrid/jqgrid.html` | Part 2b: `loadComplete` wrapper to display approximate label in pager |
| `variantopedia/grids.py` | Part 3: change `sortname` default to `'id'` |

---

## Out of Scope — Option C (separate issue)

**Stored `non_ref_count` column on `VariantZygosityCount`**: Adding a real `non_ref_count = het_count + hom_count` integer field with a plain B-tree index on `(collection_id, non_ref_count DESC)` would make the `ORDER BY non_ref_count DESC LIMIT 25` fast at any scale — PostgreSQL could use an index scan and return the top 25 rows directly without sorting. This is the proper long-term fix for sortability by zygosity count.

Requires: model field addition, migration, updating all write paths that modify `het_count`/`hom_count`.

---

## Testing

- Load the All Variants Grid page and confirm no timeout on initial load
- Confirm `EXPLAIN ANALYZE` shows one join to `snpdb_variantzygositycount`, not two
- Confirm `ORDER BY id DESC` in the query plan (default sort)
- Confirm the pager shows `"~NM"` / `"~NB"` format rather than an exact integer
- Sorting by the non_ref_count column still works (will be slow but should not timeout due to Part 2)
- Existing URL tests (`URLTestCase`) cover basic 200 response
