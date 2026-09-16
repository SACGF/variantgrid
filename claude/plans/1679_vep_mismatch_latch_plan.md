# #1679 — VEPVersionMismatchError: latch it on the VAV, gate the dispatcher, unblock explicitly

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-14
Status: draft

[#1679](https://github.com/SACGF/variantgrid/issues/1679). Related: #1642 (the vgaws `gencode_subset`
instance, closed), #1658 / #1660 (lease and reclaim, same task file), #1665 (dispatcher may not be on the
VEP box, still open - the reason the check stays worker-observed).

## The problem

A VEP version mismatch is a property of the deployment - the installed VEP or one of its data files no
longer matches the `VariantAnnotationVersion` pins - so once it is true it is true for every run against
that VAV. It is only discovered per run: `annotation/pipelines/vep.py:VEPRunner.annotate` calls
`annotation/vep_annotation.py:vep_check_command_line_version_match` immediately before dumping, the
`VEPVersionMismatchError` unwinds into the generic handler in
`annotation/tasks/annotate_variants.py:annotate_variants`, which marks the run ERROR, sends
`Annotation pipeline run <id>` to Rollbar, and re-raises; `variantgrid/celery.py:on_task_failure` then sends
the raw exception to Rollbar as well. The scheduler keeps making range locks, the dispatcher keeps leasing
and launching, and every run fails the same way - ~950 Rollbar events for two configuration problems in
Jul-Aug 2026, plus wasted `annotation_workers` slots re-dumping variants for runs that cannot succeed.
`annotation/models/models.py:VariantAnnotationVersion.get_annotation_run_blocker` already expresses "runs
against this VAV will fail" (for a missing `GeneAnnotationRelease`) but knows nothing about VEP mismatches
and has no effect on scheduling: it is consulted inside `annotate_variants` after the run has been created,
leased and launched.

## Where the issue has drifted from the code

Verified against master `ab90dadcc`:

- The check is no longer called from `dump_and_annotate_variants` in `annotate_variants.py` (the file is
  446 lines now); it moved to `annotation/pipelines/vep.py:VEPRunner.annotate` (line 85). The generic
  handler is `annotation/tasks/annotate_variants.py` line 258 (not 280), the in-run blocker check line 236
  (not 245).
- `annotation/vep_annotation.py`: `vep_check_command_line_version_match` is at 619 (not 567),
  `_vep_check_version_match` at 593 (not 555), `get_vep_version` at 323 (not 277), and the
  `gencode_subset` / `distance` snapshot-field exclusion at 594-601 (not 540-546).
- `get_annotation_run_blocker` is at `annotation/models/models.py` line 920 (not 794); the runs page
  renders it at `annotation/templates/annotation/variant_annotation_runs.html` line 191 (not 146).
- `_dispatch_for_vav` no longer exists. Dispatch is `annotation/tasks/annotation_scheduler_task.py:_dispatch_sweep`
  (line 348) over the list from `_dispatchable_variant_annotation_versions` (325), leasing through
  `_lease_across_vavs` (396) in two lanes: CREATED runs to `annotation_workers` (VEP), ANNOTATION_COMPLETED
  runs to `db_workers` (import).
- The second Rollbar event per failure is `variantgrid/celery.py:on_task_failure` (the celery
  `task_failure` signal reports any exception that is not a `RollbarIgnoreException`), not a second
  `report_message`. That shapes the fix: staying silent after the first failure means the re-raised
  exception must subclass `library/django_utils/rollbar_middleware.py:RollbarIgnoreException`.
- Item 7346 (`gencode_subset`): the exclusion landed in `244397ca7` on 2026-07-03 and the item fired until
  2026-07-07, which is consistent with vgaws deploying a few days later. Its shape is not reachable on
  current code; nothing to do here beyond confirming vgaws is past that sha.

## Data

```python
class VariantAnnotationVersion(DataArchiveMixin, SubVersionPartition):
    ...
    # #1679: latched by annotate_variants when the installed VEP no longer matches this version's pins.
    # Set once (first failure reports to Rollbar), cleared only by the re-check task. While set the
    # scheduler creates no range locks and the dispatcher launches no VEP runs against this version.
    annotation_blocked_reason = models.TextField(null=True, blank=True)
    annotation_blocked_date = models.DateTimeField(null=True, blank=True)
```

A new migration under annotation/migrations (0184_variantannotationversion_annotation_blocked): two `AddField`s,
nothing to backfill. `SubVersionPartition.save` only bumps `AnnotationVersion` on insert, so writes to
these fields are ordinary `update()` calls and never create a new `AnnotationVersion`.

## Code changes

### `annotation/models/models.py`

- `VariantAnnotationVersion.block_annotation_runs(reason) -> bool`: a conditional update,
  `objects.filter(pk=self.pk, annotation_blocked_date__isnull=True).update(annotation_blocked_reason=reason,
  annotation_blocked_date=timezone.now())`, returning whether this call made the transition. Several
  already-leased runs fail in the same wave on parallel workers; the row count is what makes exactly one
  of them the reporter.
- `VariantAnnotationVersion.unblock_annotation_runs()`: nulls both fields with `update()`.
- `get_annotation_run_blocker()`: after the existing missing-`GeneAnnotationRelease` reason, return
  `f"VEP mismatch since {annotation_blocked_date:%Y-%m-%d %H:%M}: {annotation_blocked_reason}"` when
  latched. Its docstring gains the scheduling effect (below) - it is no longer just a page warning.

### `annotation/tasks/annotate_variants.py`

- New `AnnotationRunBlockedError(RollbarIgnoreException)`: the run still fails (ERROR, visible on the runs
  page, retryable) but `on_task_failure` stays quiet. Used for both the in-run blocker check at line 236
  (currently `InvalidAnnotationVersionError`) and the mismatch below - a run leased before the latch was
  set is the wave the issue accepts; there is nothing new to say about it.
- In the `except Exception` handler of `annotate_variants`, after the `save_if_owner` ownership check and
  before the `report_message`: when `e` is a `VEPVersionMismatchError`, call
  `annotation_run.variant_annotation_version.block_annotation_runs(str(e))`. If that returns True, send one
  `report_message(level='error')` named `Annotation blocked <vav>` with the traceback and the message in
  `extra_data` and a `create_event` at ERROR; otherwise `logging.warning` only. Then
  `raise AnnotationRunBlockedError(msg) from e`. When `e` is already an `AnnotationRunBlockedError` (the
  in-run check), skip `report_message` and re-raise. Every other exception keeps today's path unchanged.
  The `finally` (lease release, dispatcher kick) is untouched.
- New task `recheck_vep_version(variant_annotation_version_id)` on `annotation_workers` (route in
  `variantgrid/settings/components/celery_settings.py` next to `annotate_variants`): calls
  `get_vep_version.cache_clear()` so this process re-runs VEP on the fake VCF rather than returning the
  value memoised before the operator's fix, then `vep_check_command_line_version_match(vav)`. On success:
  `unblock_annotation_runs()`, `create_event` at INFO, `_trigger_dispatch(vav.pk)`. On
  `VEPVersionMismatchError`: rewrite `annotation_blocked_reason` with the current message (the date stays),
  `create_event` at WARNING, return - the operator reads the result on the runs page, nothing reaches
  Rollbar. Other worker processes keep their memoised value: after a real VEP upgrade the fix is a new VAV
  (`create_new_variant_annotation_version`), and a worker restart is documented (below) for the case where
  the same VAV is meant to run again.

### `annotation/tasks/annotation_scheduler_task.py`

- `_dispatch_sweep`: the VEP lane call becomes
  `_lease_across_vavs([vav for vav in vavs if vav.get_annotation_run_blocker() is None], AnnotationStatus.CREATED, ...)`,
  logging each skipped version and its reason at warning level. The import lane, `reclaim_stalled_annotation_runs`
  and `_dispatch_counts` keep running for a blocked version: an annotated VCF produced before the mismatch
  matches the VAV (the import lane checks its header with `vep_check_annotated_file_version_match`
  anyway), importing it frees disk, and expired leases still need reclaiming. This is the same shape as
  the `has_free_disk_for_annotation` gate that already sits on that call.
- `annotation_scheduler`: skip `_handle_variant_annotation_version` for a VAV whose blocker is set (log
  the reason) and still kick `dispatch_annotation_runs` so the import lane drains. No new range locks, so
  nothing accumulates while blocked.
- The in-run check in `annotate_variants` stays as the last line of defence for a run leased in the same
  sweep the latch was set.

### `annotation/views_annotation_runs.py` and `annotation/templates/annotation/variant_annotation_runs.html`

- The blocker alert (template line 191) drops the hard-coded "Assign a GeneAnnotationRelease" sentence in
  favour of the blocker text plus, when `vav.annotation_blocked_date` is set, a
  `recheck-vep-version-{{ vav.pk }}` button and one line: the check runs on an annotation worker; after
  a VEP upgrade create a new version instead, and restart `annotation_workers` if this version is meant to
  run again. The alert also says failed runs stay in Error and are retried with the existing
  "Retry all failed runs" button once unblocked.
- The POST branch of `variant_annotation_runs`, inside the existing per-VAV loop: on
  `recheck-vep-version-<pk>`, `recheck_vep_version.si(vav.pk).apply_async()` and an INFO message that the
  re-check is queued and the page will show the result.

### Docs

- `annotation/CLAUDE.md`, under the existing mismatch gotcha: the mismatch latches
  `annotation_blocked_reason` on the VAV, the dispatcher's VEP lane and the scheduler skip a blocked VAV,
  `get_vep_version` is memoised per worker process so the re-check clears only its own cache.
- `claude/research/annotation.md` Traps paragraph: one sentence on the latch and the re-check button.
- `scripts/vg map` after the model, task and route change.

## Decisions the issue left open

- **Only the VEP lane is gated.** The issue says "skip a blocked VAV entirely"; this plan keeps the import
  lane, reclaim and the count lane running for a blocked version (reasons above). The missing-
  `GeneAnnotationRelease` blocker is gated the same way; today no such run reaches the import lane anyway
  because `annotate_variants` fails it first.
- **Failing runs raise `AnnotationRunBlockedError`, a `RollbarIgnoreException`, rather than re-raising
  `VEPVersionMismatchError`.** Re-raising the original would still produce one `on_task_failure` Rollbar
  event per run; the issue's "silent after the first" needs the ignore marker. The run still records the
  full traceback in `error_exception`.
- **The re-check never clears a latch on its own judgement, and neither does a status change.** Only a
  passing `vep_check_command_line_version_match` on a worker clears it, as the issue asks; promoting or
  demoting a VAV leaves the fields alone. A new VAV starts unlatched, which is the normal path after a VEP
  upgrade.
- **The latch is not reported again on each dispatcher sweep.** The skip is logged at warning level and
  shown on the runs page; the one Rollbar event at the transition is the alert.
- **The missing-`GeneAnnotationRelease` case gets no Rollbar event at all** once gated (today it produces
  one per run). A NEW VAV without a release is the normal state between
  `create_new_variant_annotation_version` and `link_gene_annotation_release`, so silence there is right.

## Tests (`claude/guides/testing.md`)

Fixtures: `annotation/fake_annotation.py:get_fake_vep_version` via the `_make_vav` helper in
`annotation/tests/test_variant_annotation_version_status.py`; runs from `_make_lock` / `_lease` in
`annotation/tests/test_annotation_dispatch.py`; eager `annotate_variants` with the VEP stage patched as in
`annotation/tests/test_annotation_run_lease_abort.py`.

- `test_variant_annotation_version_status.py`: `block_annotation_runs` returns True on the first call and
  False on the second, and the reason and date from the first call survive; `get_annotation_run_blocker`
  names the latched reason; `unblock_annotation_runs` clears it.
- `test_annotation_dispatch.py`: a blocked VAV's CREATED runs are not leased while another VAV's are, and
  the blocked VAV's ANNOTATION_COMPLETED run still launches on the import lane; the scheduler creates no
  range lock for a blocked VAV (patch `_handle_variant_annotation_version`, assert it is not called for
  that version and the dispatcher kick still fires).
- `test_annotation_dispatch.py`: with `vep_pipeline.vep_check_command_line_version_match` raising
  `VEPVersionMismatchError` and `report_message` patched, running `annotate_variants` eagerly for two runs
  latches the VAV once, calls `report_message` once, leaves both runs in ERROR, and the raised exception
  is a `RollbarIgnoreException`. A third run fails with the in-run blocker without any `report_message`.
- `recheck_vep_version`: the check passing clears both fields and kicks the dispatcher; the check raising
  keeps the date, rewrites the reason and reports nothing.

Not kept: a test that the migration adds two nullable fields, or that the template renders the button
(`annotation/tests/test_urls.py` already covers the page at 200 for a superuser).

## Manual verification

On a box with VEP (vg-test2), scheduler and workers running, Rollbar or its log visible:

1. `manage.py shell`: take the ACTIVE GRCh38 VAV, note `dbsnp`, set it to `dbsnp + 1` with
   `save(update_fields=["dbsnp"])`. Import a small VCF or run the scheduler from the runs page.
2. Expect: the first run to reach VEP fails, one `Annotation blocked` event, `annotation_blocked_*` set,
   the runs page warning with the re-check button, `annotation_workers` idle, the scheduler log naming the
   skipped version, no further Rollbar events. `vg status` shows the run in Error.
3. Click re-check: the latch stays, the reason is refreshed, nothing reaches Rollbar.
4. Restore `dbsnp`, click re-check: the latch clears and the dispatcher is kicked. Click "Retry all failed
   runs"; the runs finish.
5. Runs page: `vg page /annotation/variant_annotation_runs --queries` before and after the template change -
   the blocker is one call per VAV either way.

## Definition of done

- `scripts/vg tests --explain` names `annotation.tests` (migration) and the modules above pass.
- `annotation/CLAUDE.md` gotcha and `claude/research/annotation.md` trap updated; `scripts/vg docs check`
  passes; `scripts/vg map` refreshed.
- This plan's `Status:` records the outcome.
