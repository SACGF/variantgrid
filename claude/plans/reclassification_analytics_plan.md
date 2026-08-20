# Reclassification Analytics Page

Design for [#1523](https://github.com/SACGF/variantgrid/issues/1523) — an admin page mining
`ClassificationModification` history to show how classifications change over time. Answers "how often do
we reclassify?", "what is the median VUS half-life?", "is the database trending benign?"

Scope for this plan: **germline records only**. Somatic curation is captured on a separate axis
(`somatic:clinical_significance`, tier_1..tier_4) and currently holds 8 change events across 1,486
records, so the somatic view earns its own issue once that data matures. The schema keeps
`allele_origin_bucket` so the somatic axis slots in later.

---

## What the data says

Measured on SA Path prod, 2026-08-20, over 62,905 classifications / 561,015 published modifications.

**Initial classification and reclassification are different events.** Walking consecutive published
modifications yields 10,323 transitions of `clinical_significance`, and 7,541 of them start from no
value at all. Those complete same-day (`median=0d p25=0d p75=0d p90=1d`) — they are the curator filling
in the form on a record whose first published version had no significance yet. They also account for
Part 0's 8,721 NULL-significance classifications: 86.5% of that population went on to gain a value.

That leaves **2,782 genuine reclassifications**, and their shape is the product:

| from → to | events |
|---|---|
| VUS → LB | 1,158 |
| VUS → LP | 554 |
| LP → P | 518 |
| VUS → B | 129 |
| LP → VUS | 132 |
| P → LP | 120 |
| VUS → P | 93 |

**The headline metric holds up.** 1,558 of 21,161 VUSes (7.4%) have moved at least once, at a median of
709 days to first change (p25 412d, p75 1,278d, max 2,501d). VUS resolves benign over pathogenic at
roughly 2:1 (1,287 vs 647), which is the trend the page exists to show.

**Provenance is a first-class filter.** 15 of 28 labs are `external=True` — synced in from Shariant.
Their 49,022 classifications carry 1,887 change events whose modification timestamp is the *sync* date,
not the curation date, arriving in lumps (344 events on 2025-06-18, 322 on 2026-08-05). By source,
change events are 8,391 `form`, 1,908 `api`, 24 `variantgrid`; the `api` share tracks the `admin_bot`
service account at 1,915. On Shariant every lab is local and the timestamps are curation dates, which
is why the page is validated there first.

The ~29-published-version cluster (16,264 records) is Shariant sync rewriting external records: 8.5%
carry a significance change, and the cluster is absent from the internal-lab-only view.

---

## Data model

One materialised table in `classification/models/`, holding a **timeline** of germline significance
states rather than changes alone. Charts 3 and 5 need to know what every record's significance was at an
arbitrary past date, and a timeline answers that with the same rows that answer chart 1.

```
ReclassificationEvent
    classification            FK Classification (CASCADE)
    lab                       FK Lab              # denormalised for filtering
    gene_symbol               FK GeneSymbol, null # denormalised from allele_info, for chart 5
    allele_origin_bucket      CharField(1)        # "G" for everything this phase writes

    from_clinical_significance  CharField(1), null # null marks the initial classification
    to_clinical_significance    CharField(1)
    from_modification         FK ClassificationModification, null
    to_modification           FK ClassificationModification
    reclassified_date         DateField, db_index # date of to_modification
    step                      PositiveSmallInteger # 1-based position in this record's timeline
    significance_delta        SmallInteger, null   # signed ClinicalSignificance.distance, drives chart 2
```

A row with `from_clinical_significance` null is the record's initial classification. Reclassification
charts filter `from_clinical_significance__isnull=False`; denominator queries ("how many VUSes existed
on 2024-01-01") take the latest row per classification with `reclassified_date <= D`. Expected volume on
SA Path is ~11.5k rows.

Index on `(lab, reclassified_date)` and `(from_clinical_significance, to_clinical_significance)`.

Store the significance codes, and render labels through `ClinicalSignificance.SHORT_LABELS` at display
time so the page follows any future relabelling.

### Deriving a timeline

For one classification: published modifications ordered by `created`, reading significance from
`ClassificationModification.clinical_significance`, falling back to the `clinical_significance` value in
`published_evidence` normalised through
`EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).matched_options()`. The fallback covers
79 rows where the column lags the evidence — `Classification.patch_value()` sets the column during patch
(`classification/models/classification.py:1809`, which carries a TODO about moving it to publish), so
the evidence JSON is the authority.

Emit a row wherever the significance differs from the previous published version, plus one row for the
first version that carries a significance.

---

## Population

**Backfill:** management command `reclassification_events_backfill`, streaming
`ClassificationModification.objects.filter(published=True).order_by('classification_id', 'created')` with
`.iterator()` and grouping by classification. `--rebuild` recomputes from scratch. The probe walks the
same 561k rows in a couple of minutes, so a single pass is fine.

Register it as a `ManualOperation` migration per `snpdb/migrations/0188_one_off_migrate_common_filter_gnomad_versions.py`,
with a `test=` callable that reports the task when `Classification.objects.exists()` and the event table
is empty.

**Incremental:** a `classification_post_publish_signal` receiver in a new
`classification/signals/classification_hooks_reclassification.py`. The signal already delivers
`previously_published` and `newly_published`, which is exactly the pair a row needs —
`clinical_significance_change_check` in `classification_hooks_significant_change.py` compares the same
two modifications to raise the significance-change flag, and `_SignificanceChange` there is the shape to
follow. The new receiver fires for every share level, since the admin page reports on all records.

**Safety net:** a daily Celery beat entry in `variantgrid/celery.py` alongside `notify-server-status`,
calling a task on the `db_workers` queue that reconciles classifications whose latest published
modification is newer than their newest event row.

---

## Page

`/classification/reclassification_analytics`, superuser-only via `require_superuser` from
`library.django_utils`, following `classification/views/classification_view_metrics.py`. Charts are
plotly, loaded as in `classification/templates/classification/classification_dashboard.html`
(`js/lib/plotly-latest.min.js` + `js/plotly_helper.js`).

Filters, applied to every chart: organisation → lab, date range, and a lab-provenance toggle
(local labs / synced labs / all) defaulting to local, so the default view reflects curation dates.

### 1. Sankey — how significance flows

Two ranked B / LB / VUS / LP / P columns, band width proportional to event count, over rows with
`from_clinical_significance` set. ~2,782 events populate this well.

### 2. Trend over time

One mark per reclassification, coloured by `significance_delta` sign — benign-direction one way,
pathogenic-direction the other — with a rolling average. At ~400 events/year, monthly aggregation reads
better than individual dots; keep the dots available as a hover detail.

### 3. VUS reclassification rate per year

Numerator: events with `from_clinical_significance = VUS` in the period. Denominator: VUS count at
period start, from the timeline as described above. Yearly by default — the measured 7.4% lifetime rate
across 8 years means a year is the granularity that carries signal.

### 4. Time to reclassification

Histogram of `reclassified_date` minus the date of that classification's initial-classification row,
one series per starting significance. This is the strongest result in the data (VUS n=1,558, median
709d) and the fastest to build, so it lands first.

### 5. VUS burden by gene

Sortable table plus bar chart over `gene_symbol`, counting records whose current significance is VUS,
normalised by total classifications for that gene. Uses the timeline's latest row per classification;
`gene_symbol` comes from `Classification.allele_info` → `ResolvedVariantInfo.gene_symbol`
(`classification/models/classification_variant_info_models.py:145`), captured at write time.

### 6. Evidence keys driving reclassification

Diff `from_modification.published_evidence` against `to_modification.published_evidence` per event, and
bar-chart the most frequently changed keys — "is gnomAD or functional data resolving our VUSes?" 2,782
diffs is comfortable to compute on demand with caching.

---

## Phasing

1. **Model + population** — table, backfill command, `ManualOperation` migration, publish receiver, beat
   task. Verify the backfill reproduces 2,782 reclassification rows and 8,720 initial-classification
   rows on an SA Path copy.
2. **Page + charts 4, 1, 3, 2** — in that order: chart 4 is the headline number, chart 1 the clearest
   picture, chart 3 the quality metric, chart 2 the long view.
3. **Charts 5 and 6** — gene burden, then evidence-key diffs.

Validate on Shariant, where every lab is local, before deciding whether labs see their own view.

---

## Open items

- **Per-lab access** stays superuser-only until the Shariant team has used it. When a lab view arrives,
  scope it to classifications the lab *owns* so unpublished movement elsewhere stays private.
- **Curation date for synced records** — `sync/` may carry the originating lab's timestamp in the
  payload; if it does, storing it alongside `reclassified_date` lets synced labs join the time-series
  charts. The provenance filter covers the page until then.
- **Deferred to their own issues:** lab follow-on convergence (when lab A reclassifies, do others
  follow?), and a stale-VUS list of records unreviewed in over N years.
