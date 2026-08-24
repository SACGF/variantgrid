# Classification Explorer — evidence key and reclassification filters (#1395)

## Goal

Let a curator ask questions of the classification catalogue that today need a database session:

- *How have people actually filled in `splicing_assertion`? How often in a P/LP?* (#1395 as written)
- *Which records went down in severity while PM2 came off?* (the question the reclassification analytics
  page raises but cannot answer)

Both are the same shape: narrow the catalogue with a filter, get back **classifications**, then open one and
read its change history. The unit on screen is always a classification — a reclassification is something a
classification *did*, expressed as a filter over `ReclassificationEvent`, so one record appears once however
many times it has moved.

The page also answers #1395's second half — value counts and distributions for the selected keys, over
whatever the filters currently select.

## Design

### One page, one grid, two filter blocks

A new superuser page, `/classification/explorer`, whose grid is a `ClassificationColumns` subclass —
`ClassificationExplorerColumns`, one row per `ClassificationModification` (latest published). Reusing that
config brings the identifier / c_hgvs / classification / condition / date / flags columns already tuned for
this data, and `get_initial_queryset()` already routes through
`ClassificationModification.latest_for_user(published=True, exclude_withdrawn=True)`, so share level and
withdrawal are respected without new permission code.

Two filter blocks sit above the grid, both feeding `filter_queryset()`:

**A. Evidence values** — pick evidence keys (multi-select autocomplete over `EvidenceKeyMap`), then per key
say what you want: one or more option values for `SELECT` / `CRITERIA` / `BOOLEAN`, a numeric range for
`INTEGER` / `FLOAT` / `UNIT` / `AGE`, contains-text for `FREE_ENTRY` / `TEXT_AREA`, has-a-value for any key.

**B. Reclassification** — from significance (multi), to significance (multi), direction (towards benign /
towards pathogenic), event type (Initial / Reclassification / Re-evaluation), date range, points delta range,
and evidence movement: a key plus one of applied / unapplied / strengthened / weakened / changed.

The existing lab / gene symbol / allele origin / user filters from the classifications page carry over, since
they already live on `ClassificationColumns.filter_queryset()`.

### Reclassification filters are one subquery

```python
def reclassification_q(self) -> Optional[Q]:
    events = ReclassificationEvent.objects.all()
    # ... every reclassification filter narrows THIS queryset ...
    return Q(classification__in=events.values('classification_id'))
```

Two properties matter and both follow from building a single subquery:

- **one row per classification** — a record with three qualifying events still appears once
- **the conditions describe the same event** — "went down in severity *and* lost PM2" means one move that did
  both, rather than a record that dropped in 2022 and lost PM2 in 2024

`ReclassificationEvent` already indexes `(from_clinical_significance, to_clinical_significance)`,
`(event_type, reclassified_date)` and `(lab, reclassified_date)`, so the subquery is served by those.

`ReclassificationEventBuilder.tracked_classifications_qs()` covers records curated on the germline axis, so
this block selects from the germline catalogue — the same population the analytics page reports on.

### Evidence value filters read through the shared accessor

`ReclassificationTimeline._evidence_value(key)` already builds the `COALESCE(pe -> key -> 'value',
pe -> key)` expression that copes with both the dict-wrapped and bare shapes found in `published_evidence`.
Promote it to a module-level helper (`classification/models/evidence_json.py` or alongside `EvidenceMixin`)
and use it from three places: the timeline builder, the explorer filters, and the stats panel.

Each selected key is annotated once and filtered on the annotation, so a key used by both a filter and the
stats panel is extracted once per row.

Criteria-only questions have a cheaper path: `Classification.summary.criteria_labels` is already a jsonb list
of met criteria with strength (`"PM2_Supporting"`), already rendered as a grid column. Give it a GIN index
and "has PM2 met" is a containment lookup with no `published_evidence` access at all. Route criteria filters
through it, and everything else through the annotation above.

### Stats panel

Under the grid, a card per selected evidence key, computed over the filtered queryset:

| Key type | What is shown |
|---|---|
| `SELECT`, `CRITERIA`, `BOOLEAN` | value counts, in the key's own option order, with a bar per value |
| `MULTISELECT` | counts per element (`jsonb_array_elements_text`) |
| `INTEGER`, `FLOAT`, `UNIT`, `AGE` | count / min / median / max plus a histogram |
| `FREE_ENTRY`, `TEXT_AREA`, `DATE` | has-a-value against blank, plus the most common values |

Aggregated in one query per key, shaped like
`ReclassificationEventBuilder.criteria_points_by_modification()` — grouped in the database rather than
pulling evidence documents back to count in Python. Cached against the filter query string, the way
`_evidence_movement_counts` is.

This is the interactive version of what `ClassificationExportFormatterKeys` writes to JSON today, scoped to a
filter rather than the whole catalogue.

## Materialising the evidence diff

Evidence movement is the one filter with no column behind it. Today
`ReclassificationAnalytics._count_evidence_movements()` diffs both `published_evidence` blobs of every event
in SQL, aggregate-only, cached an hour — right for a bar chart, and worth promoting to a stored column now
that rows are being selected by it.

Add to `ReclassificationEvent`:

```python
evidence_changes = models.JSONField(null=True, blank=True)
"""
Which evidence keys moved between from_modification and to_modification, bucketed the way the analytics
page reads them. Written by ReclassificationEventBuilder, so a movement filter is a containment lookup.
{
  "applied": ["pp3"], "unapplied": ["pm2"], "strengthened": [], "weakened": ["ps1"],
  "changed": ["literature", "condition"],
  "strengths": {"pp3": [null, "PP"], "pm2": ["PM", null], "ps1": ["PS", "PM"]}
}
"""
```

with `GinIndex(fields=["evidence_changes"])` in `Meta.indexes`.

The bucket lists are what filters hit — "PM2 came off" is
`Q(evidence_changes__contains={"unapplied": ["pm2"]})`. `strengths` carries the from/to pair for display and
for the aggregate chart. Non-criteria keys land in `changed`, matching how
`ReclassificationAnalytics._evidence_movements()` classifies them today, so the existing chart can read the
column and the hourly diff-and-cache step retires with it.

Criteria keys are lab-namespaced. Store the key as it appears on the record; where the page is showing more
than one lab, expand a chosen group into its namespaced variants when building the `Q`, reusing
`ReclassificationAnalytics.evidence_group_by_key`.

**Writing it:** `ReclassificationEventBuilder.apply_criteria_points()` already gathers every from/to
modification pair per batch and runs one SQL pass over them. Add a second pass in the same place — the
`_count_evidence_movements` SQL with the `GROUP BY` replaced by per-event rows — and assign the result onto
each event before `bulk_create`. One diff at build time replaces every diff at read time.

**Filling it in:** `rebuild()` writes timelines from scratch, so the migration adds the column and index and
sets `ReclassificationEventBuildState.built_to = None`. The next visit to the analytics page finds more than
`RECLASSIFICATION_PAGE_BUILD_LIMIT` outstanding and hands the batch to `reclassification_events_update`,
which is the path that already exists for a large rebuild.

## Linking from the analytics page

Every chart and table on `/classification/reclassification_analytics` becomes clickable, following the
`tag_stats.html` → `variant_tags` pattern: build a URL from the page's current filters plus the clicked
element's own, and navigate.

`ReclassificationAnalytics.query_string` already serialises lab / organisation / dates / lab_external, and
the explorer reads the same parameter names, so a link is that string plus:

| Clicked | Adds |
|---|---|
| Sankey band, matrix cell | `from=V&to=LB` |
| Direction totals ("towards benign") | `direction=benign` |
| Time-to-reclassification series | `from=V` |
| Yearly activity bar | `date_from=&date_to=` (that year), `event_type=` |
| Lab league table row | `lab=<pk>` |
| VUS burden gene bar | `gene_symbol=X` |
| Evidence movement bar segment | `evidence_key=pm2&movement=unapplied&direction=benign` |
| Points transition box | `from=V&to=LB` |

`plotly_click` is already wired in that template for the evidence chart's expand toggle, so the handler has a
home. Tables get plain anchors.

The explorer reads these on page load and paints the filter widgets with them, so an arriving user sees which
question they are looking at and can widen it — the way `variant_tags` picks up `initial_tags` from GET.

## Files

| File | Work |
|---|---|
| `classification/models/classification_reclassification_models.py` | `evidence_changes` field + GIN index; builder writes it |
| `classification/migrations/0175_*.py` | add column + index, reset the build watermark |
| `classification/views/classification_explorer_view.py` | new page view, filter parsing, stats panel |
| `classification/views/classification_explorer_datatables.py` | `ClassificationExplorerColumns(ClassificationColumns)` |
| `classification/templates/classification/classification_explorer.html` | filter blocks, grid, stats cards |
| `classification/views/classification_reclassification_view.py` | evidence chart reads `evidence_changes`; link params |
| `classification/templates/classification/classification_reclassification_analytics.html` | click-through links |
| `classification/urls.py` | page + `DatabaseTableView` datatable route |
| `uicore/templates/uicore/menus/menu_bar_classifications.html` | `{% menu_item %}`, `admin_only=True` |
| `classification/models/classification.py` (or a new `evidence_json.py`) | shared evidence value accessor |
| `snpdb/migrations/…` or `classification/migrations/…` | GIN index on `Classification.summary` criteria labels |

## Steps

1. **Shared accessor** — lift `_evidence_value` out of `ReclassificationTimeline` into a helper both the
   timeline and the explorer import.
2. **`evidence_changes`** — field, index, builder pass, migration with the watermark reset. Point the
   analytics evidence chart at the column and retire the read-time diff.
3. **Explorer page, reclassification block first** — view, `ClassificationExplorerColumns`, template, menu
   item, url. This block is plain field lookups on an indexed table, so it lands and proves the page shape.
4. **Evidence value block** — key picker, per-type predicates, criteria fast path via `summary`.
5. **Stats panel** — per-key aggregates over the filtered queryset, cached on the filter string.
6. **Click-through from analytics** — link params on both sides, `plotly_click` handlers, GET-populated
   filter widgets.

## Sizing to check while building

Measure on a production-shaped database before settling the index strategy in steps 4 and 5:

- rows in `ClassificationModification` with `is_last_published=True` — sets the cost of an evidence value
  filter that has no criteria fast path, and says whether `published_evidence` wants its own
  `jsonb_path_ops` GIN index
- rows in `ReclassificationEvent`, and the typical size of one `evidence_changes` document

## Testing

- `ReclassificationEventBuilder` writes the buckets a known from/to modification pair implies: a criterion
  coming on, coming off, changing strength in each direction, and a non-criteria key changing.
- The reclassification subquery returns a classification once when it has several qualifying events, and
  matches only when one event satisfies every condition together.
- Namespaced criteria: a group chosen across labs matches each lab's own key; a single-lab view matches that
  lab's key.
- Evidence value predicates against both `published_evidence` shapes — `{"value": x}` and bare `x`.
- Stats aggregates against a small fixture where the expected counts are countable by hand, including a
  `MULTISELECT` key and a record with the key absent.
- `URLTestCase` entry for the page and its datatable route.

## Open questions

- **Sharing the page more widely.** It starts superuser-only, matching the analytics page and the
  #1395 discussion. Lab users would find the evidence value block useful over their own records;
  `latest_for_user` already scopes rows correctly, so widening it is a permission decision rather than a code
  one — worth revisiting once the page has been used.
- **Graduating the reclassification block onto the main classifications grid.** It is cheap enough to sit on
  `/classification/classifications`, where `ClassificationGroupingColumns.classification_filter_to_grouping()`
  already exists to lift a classification-level `Q` onto a grouping. Worth doing once the vocabulary settles.
- **Somatic records.** Timelines cover the germline axis today; the reclassification block picks up somatic
  records for free when `ReclassificationEvent` starts tracking `somatic:clinical_significance`.
