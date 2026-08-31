# #1273 — Liftover page: retry failed liftovers

Written by Claude Fable 5 (claude-fable-5), 2026-08-31

[#1273](https://github.com/SACGF/variantgrid/issues/1273): after fixing a bug or config problem we want to
re-run liftover for the alleles that failed, broken down by tool (e.g. "relaunch all failed bcftools
+liftover jobs").

---

## 1. Why there is currently no way to retry

Every path that creates liftover pipelines goes through `_get_build_liftover_dicts()`
(`snpdb/liftover.py:183`), which fetches `AlleleLiftover.get_failed_conversion_tools()`
(`snpdb/models/models_variant.py:1223`) and hands each allele's set of previously-failed tools to
`_liftover_using_dest_variant_coordinate()` / `_liftover_using_source_variant_coordinate()`. Those
skip any tool in the set. The failed set is "any `AlleleLiftover` with `status=ERROR` for this allele and
destination build, ever" — there is no age, no distinction between a transient failure (bcftools binary
missing, chain file misconfigured, `LiftoverRun.error()` marking a whole batch `"Liftover(N) failed"`)
and a genuine one (ClinGen has no coordinates for that build).

So once an allele has an `ERROR` row for a tool, that tool is permanently off the table for it:

- **Liftover page "Liftover variants" button** (`snpdb/views/views.py:2048` → `liftover_alleles` task →
  `Allele.missing_variants_for_build`) selects the right alleles, but every tool is skipped, so the
  `create_liftover_pipelines` call produces empty dicts and silently creates nothing. The "Missing
  Count" column stays the same and the user has no feedback.
- **Allele page "Create Variant" button** is hidden, because `VariantCard` (`classification/variant_card.py:43`)
  gates it on `allele_can_attempt_liftover()` (`snpdb/liftover.py:417`), which honours the same failed set.
- **Admin action "Liftover"** on `Allele` (`snpdb/admin.py:87`) → `liftover_alleles()` → same skip.

The only current workaround is deleting `AlleleLiftover` rows by hand in the admin, which also throws away
the error message history that the allele page's "Liftover Details" modal and the "Failed Liftover to
…" grid show.

The failed-tool memory is the right default — it's what stops a ClinGen-less allele being re-queried
on every import. The fix is a way to *deliberately override it for one tool*, keeping the earlier
failures as history. That is the first option in the issue; it's preferable to the delete-and-rerun
alternative because the error messages are the main thing you look at when diagnosing why a tool is
failing, and because a retry that fails again lands as a new `AlleleLiftover` under a new `LiftoverRun`
(`unique_together = ('liftover', 'allele')`), so history stays linear and `get_last_failed_liftover_run`
keeps pointing at the latest attempt.

---

## 2. Core: `retry_conversion_tools` through the pipeline builder

`snpdb/liftover.py`

Add a keyword argument, threaded through the three layers:

```python
def create_liftover_pipelines(user, alleles, import_source, inserted_genome_build,
                              destination_genome_builds=None,
                              retry_conversion_tools: Iterable[AlleleConversionTool] = ()):
```

`_create_liftover_pipelines_for_batch` passes it on unchanged; `_get_build_liftover_dicts` applies it
where the per-allele failed set is looked up:

```python
retry_tools = set(retry_conversion_tools)
...
failed_tools = build_failed_tools[genome_build].get(allele.pk, set()) - retry_tools
```

That is the whole behavioural change. With bcftools in `retry_tools`:

- ClinGen stays skipped for alleles it already failed on (usually "no ClinGen allele"), so retrying
  bcftools doesn't spray a fresh ClinGen `ERROR` row per allele.
- bcftools is attempted again exactly as on first run, including the `Variant.validate` and
  `bcftools_pre_liftover_error_check` pre-checks, so an allele that is still invalid gets a new `ERROR`
  row with the current reason rather than the stale one.
- Success writes a `VariantAllele` via the normal `BulkAlleleLinkingVCFProcessor` path, and
  `liftover_run_complete_signal` fires so `ImportedAlleleInfo.relink_variants` picks the variant up for
  classifications.

Alleles that have since acquired a variant in the destination build are already excluded by the
`existing_builds` check at the top of the per-allele loop, so a retry queryset can be generous.

### Selecting the alleles to retry

`snpdb/models/models_variant.py` — alongside `Allele.missing_variants_for_build`:

```python
@staticmethod
def failed_liftover_for_build(genome_build, conversion_tool) -> QuerySet['Allele']:
    """ Alleles still missing a variant in genome_build, where conversion_tool has an ERROR attempt """
    return Allele.missing_variants_for_build(genome_build).filter(
        alleleliftover__status=ProcessingStatus.ERROR,
        alleleliftover__liftover__genome_build=genome_build,
        alleleliftover__liftover__conversion_tool=conversion_tool).distinct()
```

This mirrors `get_failed_conversion_tools` (any `ERROR` ever for that tool + build), so the count shown on
the page is exactly the set of alleles the retry will re-attempt.

---

## 3. Celery task

`snpdb/tasks/liftover_tasks.py`

Both tasks gain an optional `retry_conversion_tool: str | None = None` (the enum value, since it goes
through the Celery message). Pull the queryset choice into one helper so the top-level task and the batch
task (which deliberately re-queries) agree:

```python
def _alleles_to_liftover(genome_build, retry_conversion_tool=None) -> QuerySet[Allele]:
    if retry_conversion_tool:
        return Allele.failed_liftover_for_build(genome_build, retry_conversion_tool)
    return Allele.missing_variants_for_build(genome_build)
```

`liftover_alleles` uses it for batching and passes `retry_conversion_tool` into each
`liftover_allele_batch.si(...)`. `liftover_allele_batch` uses it for the re-query and calls
`create_liftover_pipelines(..., retry_conversion_tools=[retry_conversion_tool] if retry_conversion_tool else ())`.
Everything else (batching by pk range, `log_traceback`, one task per batch) stays as is.

---

## 4. Liftover page

### View — `snpdb/views/views.py:liftover_runs`

**POST**: alongside the existing `liftover_to_{build}` buttons, accept `retry_{build}_{tool}` where
`tool` is the `AlleleConversionTool` value. Parse both in one loop over
`genome_builds × AlleleConversionTool`, queue `liftover_alleles.si(username, build.name, tool)` and add
an INFO message `"Retrying {tool label} liftover to {build} for {n} alleles"`. Keep the existing
`ValueError` for an unrecognised POST.

**GET context**: per destination build, the retryable counts by tool — one query per build:

```python
retry_counts = {}
for genome_build in genome_builds:
    missing_qs = Allele.missing_variants_for_build(genome_build)
    qs = AlleleLiftover.objects.filter(liftover__genome_build=genome_build,
                                       status=ProcessingStatus.ERROR,
                                       allele__in=missing_qs)
    counts = qs.values_list("liftover__conversion_tool").annotate(n=Count("allele", distinct=True))
    retry_counts[genome_build.name] = [(AlleleConversionTool(ct), n) for ct, n in counts if n]
```

Only tools with a non-zero count are listed, in `AlleleConversionTool` order.

### Template — `snpdb/templates/snpdb/liftover/liftover_runs.html`

In each existing "Failed Liftover to {{ genome_build }}" block, above the failures datatable, add a
small form (the page already has one POST form; make this a second `<form method="post">` scoped to the
block, or extend the existing one — either is fine, the button `name` carries the intent):

```
Tool                  | Alleles still missing | Action
BCFtools/liftover     | 1,234                 | [Retry failed liftovers]
ClinGen Allele Registry | 56                  | [Retry failed liftovers]
```

Add a one-line `text-info` note under the heading: a retry re-attempts that tool for every allele it has
failed on and still lacks a variant for; earlier attempts are kept and shown in the failures table
below. Bootstrap 4 (`btn btn-secondary`, `table`), matching the "Alleles Missing Variants" table.

---

## 5. Allele page and admin — same override, single-allele scale

These are the two other places a human explicitly asks for a liftover, and both currently go quiet after
every tool has failed. Same kwarg, so each is a one-liner:

- `variantopedia/views.py:create_variant_for_allele` — pass
  `retry_conversion_tools=list(AlleleConversionTool)`: the user clicked "Create Variant", so try every
  tool again.
- `classification/variant_card.py` / `allele_can_attempt_liftover()` — for the button to *appear* on an
  all-failed allele, `allele_can_attempt_liftover` needs the same override. Give it a
  `retry_conversion_tools` parameter, forward it into the two `_liftover_using_*` calls (they already
  accept `failed_tools`; pass `failed_tools=<fetched set> - retry_tools`), and have `VariantCard` call it
  with all tools. The button then shows whenever a tool *could* produce a coordinate, which is the
  question the user is actually asking.
- `snpdb/admin.py:AlleleAdmin.liftover` → `snpdb/liftover.py:liftover_alleles()` — add
  `retry_conversion_tools` to `liftover_alleles()` and have the admin action pass all tools. An admin
  selecting specific alleles and clicking "Liftover" wants them retried.
- `snpdb/management/commands/liftover_alleles.py` — add `--retry-tool <value>` (choices from
  `AlleleConversionTool`) forwarding to the task, so ops can kick a retry from a shell after a deploy.

---

## 6. Tests — `snpdb/tests/test_liftover.py`

`TestLiftover.setUpTestData` already has an allele with a GRCh37 variant and a ClinGen allele covering
both builds. Add:

1. **`test_retry_conversion_tools_overrides_failed_set`** — create a `LiftoverRun(conversion_tool=BCFTOOLS_LIFTOVER,
   genome_build=grch38)` plus an `ERROR` `AlleleLiftover` for the allele, and a ClinGen `ERROR` too so
   every tool is failed. Assert `_get_build_liftover_dicts([allele], grch37, [grch38])` returns empty
   dicts; then with `retry_conversion_tools={BCFTOOLS_LIFTOVER}` assert the second dict has a
   `BCFTOOLS_LIFTOVER` entry for grch38 carrying a coordinate (and no `CLINGEN_ALLELE_REGISTRY` entry —
   ClinGen stays skipped).
2. **`test_failed_liftover_for_build`** — two alleles with bcftools `ERROR` rows to grch38; give one a
   grch38 `VariantAllele`. Assert the queryset returns only the other, and returns nothing for
   `CLINGEN_ALLELE_REGISTRY`.
3. **`test_allele_can_attempt_liftover_with_retry`** — all tools failed → `False`; with
   `retry_conversion_tools=list(AlleleConversionTool)` → `True`.

A `URLTestCase` POST to `liftover_runs` with `retry_GRCh38_BL` under `CELERY_TASK_ALWAYS_EAGER` isn't
needed — the view logic is button-name parsing; the pipeline behaviour is covered by (1).

---

## 7. Files touched

| File | Change |
|---|---|
| `snpdb/liftover.py` | `retry_conversion_tools` kwarg on `create_liftover_pipelines`, `_create_liftover_pipelines_for_batch`, `_get_build_liftover_dicts`, `liftover_alleles`, `allele_can_attempt_liftover` |
| `snpdb/models/models_variant.py` | `Allele.failed_liftover_for_build()` |
| `snpdb/tasks/liftover_tasks.py` | `retry_conversion_tool` on both tasks, `_alleles_to_liftover()` helper |
| `snpdb/views/views.py` | `liftover_runs`: parse `retry_{build}_{tool}` POST, `retry_counts` context |
| `snpdb/templates/snpdb/liftover/liftover_runs.html` | per-build retry table + note |
| `variantopedia/views.py` | `create_variant_for_allele` retries all tools |
| `classification/variant_card.py` | `allele_can_attempt_liftover(..., retry_conversion_tools=list(AlleleConversionTool))` |
| `snpdb/admin.py` | admin action retries all tools |
| `snpdb/management/commands/liftover_alleles.py` | `--retry-tool` option |
| `snpdb/tests/test_liftover.py` | three tests above |
| `variantgrid/templates/default_templates/changelog.html` | entry for #1273 |

No migration: no model fields change.
