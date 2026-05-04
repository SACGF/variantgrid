# Plan: SACGF/variantgrid #1536 — `DataArchiveMixin` + audit + VCF/GeneCoverage re-import

Step 2 of the umbrella ticket #1539. #1537 (partition pg_dump pipeline) is a separate ticket that consumes the contracts defined here.

## Scope summary

- Add `DataArchiveMixin` (abstract) with four fields to: **VCF**, **VariantAnnotationVersion**, **GeneCoverageCollection**. (There is no `VariantTranscriptAnnotationVersion`; `VariantTranscriptAnnotation` and `VariantGeneOverlap` are partitioned under VAV — see `RECORDS_BASE_TABLE_NAMES` at `annotation/models/models.py:617`.)
- Wire archive/restore actions for **VCF** (detail-page, write-access users) and **GeneCoverageCollection** (Django admin only — GCC has no Guardian permissions of its own and is accessed via Sample permissions in some places; if access surfaces as a problem it gets fixed under a separate issue).
- VAV gets the mixin only — no archive action in this ticket; #1537 wires that. Tests for VAV-archive guards stamp `data_archived_date` directly.
- Single chokepoint for archive-aware queryset behaviour (`Cohort.cohort_genotype_collection` raises `DataArchivedError` — loud failure over silent empty result, so callers surface "data archived" instead of returning zero rows).
- Detail-page banner (read-only) on VCF; listing badges on existing grids (VCF, GeneCoverageCollection, Cohort, Analyses).
- Eventlog events on archive/restore.

**Out of scope:** per-sample archive, custom (non-VCF) cohort archive, Sample/Cohort/Trio mixin, #1537 partition dump pipeline, #1534 declarative partitions.

**Sub-cohort note (Q8 follow-up):** sub-cohorts inherit `cohort_genotype_collection` from their parent. Archiving the parent VCF makes the sub-cohort's CGC throw `DataArchivedError` automatically — no separate cascade needed. Cross-VCF custom cohorts copy CGC data into their own CGC at build time, so they are insulated from later parent-VCF archive (per spec — per-sample archive is deferred).

---

## 1. Mixin + migrations

### `library/django_utils/data_archive_mixin.py` (new)

```python
class DataArchiveMixin(models.Model):
    data_archived_date = models.DateTimeField(null=True, blank=True)
    data_archived_by = models.ForeignKey(
        User, null=True, blank=True, on_delete=SET_NULL, related_name="+"
    )
    data_archive_reason = models.TextField(null=True, blank=True)
    data_restorable_from = models.TextField(null=True, blank=True)

    class Meta:
        abstract = True

    @property
    def data_archived(self) -> bool:
        return self.data_archived_date is not None
```

`related_name="+"` is required: the same User would otherwise back-reference four colliding `*_set` accessors.

Top-level imports clean — no app coupling.

### Apply mixin to:

- `snpdb/models/models_vcf.py:61` — VCF
- `annotation/models/models.py:613` — VariantAnnotationVersion
- `genes/models.py:2184` — GeneCoverageCollection

### Migrations (additive, three files)

- `snpdb/migrations/0XXX_vcf_data_archive.py` — `AddField` × 4 on `snpdb_vcf`
- `annotation/migrations/0XXX_data_archive.py` — `AddField` × 4 on `annotation_variantannotationversion`
- `genes/migrations/0XXX_genecoveragecollection_data_archive.py` — `AddField` × 4

All `AddField` only, no data migration. Auto-generate with `makemigrations` then verify reverse-runs cleanly on a `--keepdb` test DB.

---

## 2. Archive-path settings + helper

### Settings

In `variantgrid/settings/components/annotation_settings.py`:

```python
ANNOTATION_ARCHIVE_DIR = get_secret("ANNOTATION.ARCHIVE_DIR", "/data/annotation_archive/")
```

### Helper

`annotation/archive_paths.py` (new):

```python
def get_annotation_archive_path(partition_filename: str) -> str:
    db_name = settings.DATABASES["default"]["NAME"]
    base = os.path.join(settings.ANNOTATION_ARCHIVE_DIR, db_name)
    os.makedirs(base, exist_ok=True)
    return os.path.join(base, partition_filename)
```

The DB-name namespace prevents test/prod/demo deploys on shared NFS/EFS from clobbering each other. #1536 owns the helper + setting; #1537 consumes it for the actual dump/restore I/O.

---

## 3. Audit — single chokepoint design

The original "edit every source node" design is overkill. Source nodes (Sample/Trio/Cohort/Pedigree/Quad) all funnel through `Cohort.cohort_genotype_collection`. Intercept there.

### Primary chokepoint: `Cohort.cohort_genotype_collection` (`snpdb/models/models_cohort.py:182`)

Raise a new `DataArchivedError` (caught at upstream entry points) if `self.vcf and self.vcf.data_archived`. The intent is **loud failure over silent empty result** — every consumer (detail pages, source-node helpers, `VCF.get_variant_qs`, `Sample.get_variant_qs`, `Sample.get_genotype`, `genotype_vcf_processor_factory`, etc.) sees the exception and surfaces "data archived" rather than returning zero rows.

`cohort_genotype_collection` is currently a `cached_property` — leave that, but the underlying lookup raises before the result is cached so a later restore (which clears `data_archived_date`) starts working again on next access.

Custom cross-VCF cohorts have `cohort.vcf is None` and copy CGC data at build time, so they are not gated. Sub-cohorts inherit from their parent and the chokepoint catches them automatically through the parent's `vcf`.

### Source-node shared helper

In `analysis/models/nodes/sources/` — helper added to `AbstractCohortBasedNode` (covers `CohortNode`, `TrioNode`, `PedigreeNode`, `QuadNode`) and to `SampleMixin` (covers `SampleNode`):

- `_get_node_arg_q_dict` wrapper: catch `DataArchivedError` → return `{None: {str(q_none): q_none}}` (mirrors the existing "no sample" path)
- `_get_configuration_errors` wrapper: surface `f"Source data archived {date} by {user}: {reason}"` → drives `NodeStatus.ERROR_CONFIGURATION` via existing `get_status_from_errors` machinery (`analysis_node.py:723`)

`AllVariantsNode` is unaffected — not pinned to VCF/Sample/Cohort/Trio.

### Annotation queryset builders

| Code path | File:line | Planned |
|---|---|---|
| `get_variant_queryset_for_annotation_version` | `annotation/annotation_version_querysets.py:31` | Guard on `annotation_version.variant_annotation_version.data_archived` → `qs.none()` + `report_message` (from `library.log_utils`) |
| `get_variants_qs_for_annotation` | `:52` | Same guard |
| `VariantTranscriptAnnotation.get_overlapping_genes_q` (queries `VariantGeneOverlap`) | `annotation/models/models.py:1668–1676` | Filters `VariantGeneOverlap.objects.filter(version=variant_annotation_version, ...)`. If `variant_annotation_version.data_archived` → return `Q(pk__in=[])` + `report_message`. Used by gene-list node filtering. |

(`VariantTranscriptAnnotation` and `VariantGeneOverlap` are partitioned under VAV, so a VAV archive guard is sufficient — there is no separate VTAV class.)

### `Analysis.get_errors` (consumed by `_get_analysis_errors` at `analysis_node.py:694`)

Add archive check: if pinned `AnnotationVersion`'s VAV is archived → analysis-level error. Maps to `NodeStatus.ERROR_WITH_PARENT` (`analysis_node.py:725`). One change propagates to all node types.

`Analysis.get_warnings` (`analysis/models/models_analysis.py:168`): also surface archived VAV as warning so listing/detail UI shows status without re-running.

### `GeneCoverageCollection` chokepoint

Identify the read entry point during implementation (likely `GeneCoverageCollection.get_uncovered_gene_symbols()` at `genes/models.py:2315` or related grid/queryset builders). Add guard: `qs.none()` + `report_message` if `collection.data_archived`. Single intercept, same pattern as VCF.

### Other independent paths (not under the chokepoint)

| Code path | File:line | Planned |
|---|---|---|
| `variant_sample_information.show_variant_sample_information` | `snpdb/variant_sample_information.py:85–120` | Filter out rows where joined VCF archived. Joins `variant→cohortgenotype→...→vcf` directly, doesn't go through `cohort_genotype_collection`. **Critical** — only place not covered by the chokepoint. |
| `VCFAlleleSource.get_variant_qs` | `snpdb/models/models_vcf.py:587–598` | `qs.none()` when `vcf.data_archived` — prevents LiftoverRun on archived data |
| `VCF.filter_for_user` | `snpdb/models/models_vcf.py:120` | Add `include_archived: bool = True`. Default keeps detail pages working; autocompletes/source-node selectors pass `False` so archived rows aren't pickable for new analyses. Pinned existing analyses still resolve. |
| Trio/Sample/Cohort autocompletes | `snpdb/views/views_autocomplete.py` | Thread `include_archived=False` via the inner VCF lookup |
| `auto_run_analyses_for_vcf` / `_for_sample` | `analysis/tasks/auto_analysis_tasks.py:9,18` | Defensive early return when archived. Also consumed by the restore "no auto-launch" requirement (Q8) — see §4. |

### Cascades that don't need touching (called out explicitly)

- `vcf_pre_delete_handler` (`snpdb/models/models_vcf.py:223`) only deletes the cohort. Archive does not delete the VCF row, so the signal does not fire and the cohort is preserved (correct — the cohort is needed for sub-cohort relationships and sub-cohort archive cascades).
- Guardian permissions unchanged.

---

## 4. Archive + restore actions

### Archive

The whole point of archive is **keep the metadata, drop the big underlying data**. Calling `vcf.delete()` would destroy the VCF row and its mixin stamp along with it — wrong.

**Reuse existing cleanup, no new extraction.** The data-cleanup we want already exists, split between two pieces:
- `VCF.delete_internal_data()` (`snpdb/models/models_vcf.py:194`) — drops VCF filters, sample stats, `uploadedvcfpendingannotation`, `variantzygositycountforvcf_set`, and partitions. **Important:** it currently drops then **recreates** the partitions (so reload-in-place works). Archive does not want recreated empty partitions — split the recreate out into a small helper (e.g. `_recreate_cohort_genotype_partitions()`) and have `delete_internal_data()` call it. `archive_vcf()` skips the recreate step.
- `update_all_variant_zygosity_counts_for_vcf(vcf, '-')` (`snpdb/variant_zygosity_count.py`) — already used by `reload_vcf_task` and `remove_soft_deleted_vcfs_task` before they touch internal data.

`snpdb/archive.py` (new):

```python
class ArchivePreconditionError(Exception):
    pass


def archive_vcf(vcf: VCF, user: User, reason: str = "") -> None:
    if vcf.data_archived:
        return
    try:
        uf = vcf.uploadedvcf.uploaded_file
        upload_path = uf.get_filename()
    except (UploadedVCF.DoesNotExist, AttributeError):
        raise ArchivePreconditionError(
            "Uploaded file doesn't exist — cannot archive. "
            "If you want to free up space you can permanently delete."
        )
    if not os.path.exists(upload_path):
        raise ArchivePreconditionError(
            "Uploaded file doesn't exist — cannot archive. "
            "If you want to free up space you can permanently delete."
        )
    try:
        update_all_variant_zygosity_counts_for_vcf(vcf, '-')
    except Exception:
        log_traceback()  # mirrors existing reload/soft-delete tolerance
    vcf.delete_internal_data(recreate_partitions=False)  # drop, do not recreate
    vcf.data_archived_date = timezone.now()
    vcf.data_archived_by = user
    vcf.data_archive_reason = reason
    vcf.data_restorable_from = upload_path
    vcf.save()  # row stays; only the heavy data is gone
    eventlog.create_event(user, "vcf_archived", details={
        "vcf_id": vcf.pk, "vcf_name": vcf.name,
        "reason": reason, "restorable_from": upload_path,
    })
```

The view handler catches `ArchivePreconditionError` and re-renders the detail page with the message. The template's pre-check duplicates the same logic so the button is hidden in the common case — server-side raise is the safety net for race conditions (file deleted between page render and POST).

`VCF.delete()` is **not** overridden — Django's cascade plus the existing `vcf_pre_delete_handler` already drops everything when a delete actually happens. Reload-in-place callers (`reload_vcf_task`) keep using `delete_internal_data()` with the default `recreate_partitions=True`.

For GeneCoverageCollection: GCC's existing `gene_coverage_collection_pre_delete_handler` (`genes/models.py:2326`) calls `delete_related_objects()` (drops the partition tables). Archive calls `delete_related_objects()` directly without deleting the GCC row — same shape, no extraction needed.

Parallel `archive_gene_coverage_collection(gcc, user, reason)` — same pattern. `data_restorable_from` records `gcc.path` (the existing source path field — backend or uploaded). No sha256 capture: GCC restore re-imports like a VCF rather than restoring from a dump (per Q9), so the integrity model matches VCF re-import.

VAV: **no archive helper in this ticket.** Mixin lands; #1537 wires the action. Tests for the queryset / Analysis-error guards stamp `data_archived_date` directly to exercise the paths end-to-end.

### Restore

Per Q4, restore is purely VCF data re-import. Variants and annotations are untouched by VCF archive (only `CohortGenotype*` partitions go away). No re-annotation step.

```python
def restore_vcf(vcf: VCF, user: User) -> AsyncResult:
    if not vcf.data_archived:
        raise ValueError("VCF is not archived")
    if not vcf.data_restorable_from:
        raise ValueError("No restore path recorded")
    if not os.path.exists(vcf.data_restorable_from):
        raise ValueError(f"Restore source missing: {vcf.data_restorable_from}")
    expected = vcf.uploadedvcf.uploaded_file.sha256_hash
    if expected:
        actual = file_sha256sum(vcf.data_restorable_from)
        if actual != expected:
            raise ValueError("sha256 mismatch: file changed since archive")
    set_vcf_and_samples_import_status(vcf, ImportStatus.IMPORTING)
    upload_pipeline = vcf.uploadedvcf.uploaded_file.uploadpipeline
    return retry_upload_pipeline(upload_pipeline)
```

Reuses `retry_upload_pipeline` (`upload/uploaded_file_type.py:90`) → routes to `reload_vcf_task` (`upload/tasks/vcf/genotype_vcf_tasks.py:194`). The "import into existing VCF.id" mode the spec wants **already exists** at `genotype_vcf_tasks.py:42` — `ImportCreateVCFModelForGenotypeVCFTask.process_items` reuses `upload_pipeline.uploadedvcf.vcf` when set, calling `configure_vcf_from_header` to refresh format fields.

Clear `data_archived_date` inside the existing success path (`ImportGenotypeVCFSuccessTask.process_items`), not in the entry helper. Emit `vcf_restored` eventlog event there.

**Auto-launch behaviour on restore (Q6):** the existing reload-side hook `reload_auto_analyses_for_vcf` (`analysis/tasks/auto_analysis_tasks.py:28`) refreshes node counts on existing analyses without auto-launching new ones, and `auto_run_analyses_for_vcf` already takes `skip_already_analysed`. Verify during implementation that `vcf_import_success_signal` consumers do the right thing for a VCF that already has analyses — if so, no `is_restore` flag is needed. If a consumer turns out to auto-launch unconditionally, surface that and decide on a flag at that point rather than threading one preemptively.

Parallel `restore_gene_coverage_collection(gcc, user)` — GCC restores via the existing gene-coverage import pipeline (no partition dump/restore — Q9). Locate the reload entry point during implementation and reuse.

---

## 5. Permissions and entry points

### VCF — write-access users via the GUI

Archive/restore for **VCF** is available to any user with `can_write` on the object (Guardian permission via the existing `VCF.can_write`). Buttons live on the detail page, not in Django admin.

**Detail-page archive/restore buttons** on `view_vcf.html`:
- Show "Archive" button when `obj.can_write(user) and not obj.data_archived` and the precondition check below passes.
- Show "Restore" button when `obj.can_write(user) and obj.data_archived` and the recorded `data_restorable_from` exists on disk.
- Confirmation modal (Bootstrap 4 `data-toggle="modal"`) prompts for `data_archive_reason` on archive.

**Pre-archive precondition check.** Before showing the archive button (and re-checked server-side in the view handler):
- Resolve the uploaded source path: `vcf.uploadedvcf.uploaded_file.get_filename()`.
- If `os.path.exists(path)` is False → button is disabled / replaced with a message:
  > **Uploaded file doesn't exist — cannot archive.** If you want to free up space you can permanently delete this VCF.
  Render a "Delete permanently" button next to the message (links to the existing delete flow that already exists for write-access users).
- If `vcf.uploadedvcf` doesn't exist at all (legacy VCFs imported without an uploaded-file record) → same treatment: archive disabled, permanent-delete offered.

This is the gate that keeps "archive" honest — archive is only a meaningful operation if restore is possible. When the source has gone missing, the right answer is permanent delete, not archive-and-pray.

**Views** (`snpdb/views/views.py`):
- `archive_vcf_view(request, vcf_id)` — POST, `obj.check_can_write(request.user)`, runs the precondition check, calls `archive_vcf(...)`, redirects to detail page.
- `restore_vcf_view(request, vcf_id)` — POST, `obj.check_can_write(request.user)`, calls `restore_vcf(...)`, redirects to detail page.

URL routes added to `snpdb/urls/urls.py`.

### GeneCoverageCollection — Django admin only

GCC has no Guardian permissions of its own and is accessed via Sample permissions in some places (Q2). No detail-page UI work in this ticket. Archive/restore are wired as Django admin actions on `GeneCoverageCollectionAdmin`. If access surfaces as a problem, it gets fixed under a separate issue.

### VAV — admin only (no action wired in this ticket)

VAV is a global system-level object, not user-owned. Archive action lands with #1537 — when it does, it goes through Django admin (superuser-only), since regular write users have no business archiving annotation versions.

For this ticket, the mixin lands and admin list views show the four mixin fields. No action wired.

### Django admin (oversight + GCC actions)

- `snpdb/admin.py` — `VCFAdmin` with `list_display` and `list_filter` covering the mixin fields. No archive action — sysadmins use the GUI like everyone else.
- `annotation/admin.py` — `VariantAnnotationVersionAdmin` (mixin fields visible, no actions until #1537).
- `genes/admin.py` — `GeneCoverageCollectionAdmin` with `list_display` covering mixin fields **and** `actions=["archive_selected", "restore_selected"]` wiring the helpers.

---

## 6. Eventlog integration

Use the project's existing `eventlog` app (per `claude/research/eventlog.md`). Locate the event-creation pattern during implementation (`Event.objects.create(...)` or a `create_event(...)` helper).

Events emitted from inside the archive/restore helpers (post-success, not via model `save()` so emits tie to deliberate operator actions):

- `vcf_archived` — `{vcf_id, vcf_name, reason, restorable_from}`, actor = archive user
- `vcf_restored` — `{vcf_id, restored_from}`, actor = restore user
- `gene_coverage_collection_archived` / `_restored` — same shape

VAV: no events in this ticket; #1537 wires them when the archive action lands.

**Not using auditlog** — auditlog scope is analysis-specific in this codebase, not for system-wide model change tracking.

---

## 7. UI

### Detail-page banner + archive/restore controls

`snpdb/templates/snpdb/data/_data_archived_banner.html` (new):

```html
{% if obj.data_archived %}
<div class="alert alert-warning" role="alert">
    <strong>Data archived</strong> on {{ obj.data_archived_date|date:"Y-m-d" }}
    by {{ obj.data_archived_by|default:"unknown" }}.
    {% if obj.data_archive_reason %}<p>{{ obj.data_archive_reason }}</p>{% endif %}
    {% if obj.data_restorable_from %}
        <p>Restorable from: <code>{{ obj.data_restorable_from }}</code></p>
    {% endif %}
    {% if can_write and restore_source_exists %}
        <form method="post" action="{% url 'restore_vcf' obj.pk %}">{% csrf_token %}
            <button class="btn btn-primary" type="submit">Restore</button>
        </form>
    {% endif %}
</div>
{% endif %}
```

Archive control (only when not archived, on the detail page action area):

```html
{% if can_write and not obj.data_archived %}
    {% if uploaded_file_exists %}
        <button class="btn btn-warning" data-toggle="modal" data-target="#archiveModal">Archive data</button>
        {# modal POSTs reason to archive_vcf_view #}
    {% else %}
        <div class="alert alert-info">
            Uploaded file doesn't exist — cannot archive.
            If you want to free up space you can permanently delete this VCF.
        </div>
        {# permanent-delete button reuses the existing delete flow #}
    {% endif %}
{% endif %}
```

Bootstrap 4 — `data-toggle`, `data-target`, no `data-bs-*` anywhere.

View context contributes:
- `can_write = obj.can_write(request.user)`
- `uploaded_file_exists = <precondition check from §5>`
- `restore_source_exists = obj.data_restorable_from and os.path.exists(obj.data_restorable_from)`

Include banner from:
- `snpdb/templates/snpdb/data/view_vcf.html`
- The GeneCoverageCollection detail template (locate during implementation)

### Listing badges — keep existing tables

- VCF list (jqGrid in `snpdb/grids.py`): add an "Archived" badge via cell formatter on an existing column (or new read-only column). No replacement.
- GeneCoverageCollection list: same pattern — locate the existing grid and add a badge.
- Cohort listing: surface "Source data archived" when the cohort's `cohort_genotype_collection` raises `DataArchivedError` (sub-cohorts of archived VCFs — Q8). Cell formatter wraps the access in try/except, badges accordingly.
- Analyses listing (`analysis/grids.py:911` `AnalysesColumns`, already DataTables): add a "Source data archived" badge column annotated via `Exists(...)` subquery walking `analysisnode_set` for archived FK targets. Add a "Hide analyses with archived data" filter wired into `get_q_list` (`:942`). Detail-side analysis editor surfaces archive via `Analysis.get_warnings` for free.

---

## 8. Tests

All `django.test.TestCase`, run with `python3 manage.py test --keepdb`.

### `snpdb/tests/test_data_archive_mixin.py` (new)

- Mixin fields present on VCF, VAV, GCC; `data_archived` property reflects `data_archived_date`.
- `archive_vcf` stamps fields, populates `data_restorable_from` from upload path, decrements zygosity counts, drops partitions without recreate.
- After archive, partition table for the VCF's CGC is gone; VCF row + samples + cohort still exist.
- `archive_vcf` is idempotent — second call no-ops.
- Eventlog event emitted with correct payload.

### `snpdb/tests/test_filter_for_user_archive.py`

- `filter_for_user(user)` includes archived by default (detail page works).
- `filter_for_user(user, include_archived=False)` excludes archived.

### `analysis/tests/test_source_node_archive_tolerance.py`

- For each of `SampleNode`, `TrioNode`, `CohortNode`, `PedigreeNode`, `QuadNode`: archive the source's VCF, assert configuration error mentioning "archived", `_get_node_arg_q_dict` returns `q_none`, queryset returns 0 rows without raising at the node level.
- `Analysis.get_errors` surfaces VAV-archived as analysis-level error (stamp `vav.data_archived_date` directly — no archive helper for VAV in this ticket).

### `annotation/tests/test_annotation_queryset_archive.py`

- `get_variant_queryset_for_annotation_version` returns `qs.none()` and emits `report_message` when VAV archived (stamped directly).
- `get_variants_qs_for_annotation` same.
- `VariantTranscriptAnnotation.get_overlapping_genes_q` returns empty Q when VAV archived.

### `snpdb/tests/test_archive_vcf_view.py` (new)

- Archive button hidden when user lacks `can_write`.
- Archive button hidden + "uploaded file doesn't exist" message shown when source path missing on disk.
- POSTing archive without write access → 403.
- POSTing archive when source file missing → re-renders with `ArchivePreconditionError` message; archive does not run.
- POSTing archive happy path → `archive_vcf` runs, eventlog event emitted, redirect.

### `upload/tests/test_restore_vcf.py`

- Rejects when `data_restorable_from` missing on disk.
- Rejects on sha256 mismatch.
- Calls `retry_upload_pipeline`; success handler clears `data_archived_date` and emits `vcf_restored`.
- Restore button hidden when user lacks `can_write`.

### `genes/tests/test_gene_coverage_archive.py`

- `archive_gene_coverage_collection` stamps fields, drops the related partition tables via `delete_related_objects()`, eventlog event emitted.
- Restore validates path, calls the gene-coverage reload pipeline (located during implementation).
- Queryset chokepoint returns empty when GCC archived.
- Django admin actions `archive_selected` / `restore_selected` round-trip end-to-end.

---

## 9. File-by-file change list

### Mixin / migrations
- `library/django_utils/data_archive_mixin.py` (new)
- `snpdb/models/models_vcf.py:61` — VCF inherits mixin
- `annotation/models/models.py:613` — VAV inherits mixin
- `genes/models.py:2184` — GeneCoverageCollection inherits mixin
- `snpdb/migrations/0XXX_vcf_data_archive.py` (new)
- `annotation/migrations/0XXX_data_archive.py` (new — VAV only)
- `genes/migrations/0XXX_genecoveragecollection_data_archive.py` (new)

### Settings + path helper
- `variantgrid/settings/components/annotation_settings.py` — `ANNOTATION_ARCHIVE_DIR`
- `annotation/archive_paths.py` (new) — `get_annotation_archive_path()`

### Audit (chokepoint + independent paths)
- `snpdb/models/models_cohort.py:182` — `Cohort.cohort_genotype_collection` raises `DataArchivedError` when `self.vcf and self.vcf.data_archived`
- `analysis/models/nodes/sources/` — helper on `AbstractCohortBasedNode` (covers `CohortNode`/`TrioNode`/`PedigreeNode`/`QuadNode`) and on `SampleMixin` (covers `SampleNode`) to catch + surface error
- `annotation/annotation_version_querysets.py:31,52` — VAV archive guards
- `annotation/models/models.py:1668` — `VariantTranscriptAnnotation.get_overlapping_genes_q` VAV archive guard (covers `VariantGeneOverlap` queries)
- `analysis/models/models_analysis.py:168` — `Analysis.get_warnings` adds VAV warning
- `analysis/models/models_analysis.py` — `Analysis.get_errors` adds VAV error
- `snpdb/variant_sample_information.py:85–120` — filter archived VCFs
- `snpdb/models/models_vcf.py:587–598` — `VCFAlleleSource.get_variant_qs` archived-check
- `snpdb/models/models_vcf.py:120` — `VCF.filter_for_user(include_archived=True)`
- `snpdb/views/views_autocomplete.py` — pass `include_archived=False`
- `analysis/tasks/auto_analysis_tasks.py:9,18` — defensive guards
- `genes/models.py` — `GeneCoverageCollection` chokepoint guard (read entry point located during implementation, likely `get_uncovered_gene_symbols`)

### Archive/restore helpers + views
- `snpdb/archive.py` (new) — `archive_vcf`, `restore_vcf`, `ArchivePreconditionError`, `DataArchivedError`
- `genes/archive.py` (new) — `archive_gene_coverage_collection`, `restore_gene_coverage_collection`
- `snpdb/models/models_vcf.py` — split `_recreate_cohort_genotype_partitions()` out of `delete_internal_data()`; add `recreate_partitions: bool = True` kwarg so archive can skip the recreate step
- `snpdb/views/views.py` — `archive_vcf_view`, `restore_vcf_view` (`check_can_write` gated)
- `snpdb/urls/urls.py` — URL routes
- `upload/tasks/vcf/genotype_vcf_tasks.py` — `ImportGenotypeVCFSuccessTask` clears `data_archived_*` fields on success and emits `vcf_restored` if the VCF was archived

### Admin
- `snpdb/admin.py` — `VCFAdmin` with mixin fields in `list_display` / `list_filter`. No archive action.
- `annotation/admin.py` — `VariantAnnotationVersionAdmin` (mixin fields visible, no actions until #1537)
- `genes/admin.py` — `GeneCoverageCollectionAdmin` with mixin fields in `list_display` and `actions=["archive_selected", "restore_selected"]` calling the helpers

### Templates
- `snpdb/templates/snpdb/data/_data_archived_banner.html` (new — archived banner + restore button)
- `snpdb/templates/snpdb/data/_archive_action.html` (new — archive button + "uploaded file missing" message + permanent-delete fallback)
- `snpdb/templates/snpdb/data/view_vcf.html` — include both partials
- `snpdb/grids.py` — VCF jqGrid archived badge
- Cohort grid — archived badge (sub-cohort cascade — Q8)
- GeneCoverageCollection grid — archived badge (locate)
- `analysis/grids.py:911,942` — `AnalysesColumns` archive badge column + filter

### Eventlog
- `snpdb/archive.py` and `genes/archive.py` emit events (call existing eventlog helper)

### Tests
- `snpdb/tests/test_data_archive_mixin.py` (new)
- `snpdb/tests/test_filter_for_user_archive.py` (new)
- `analysis/tests/test_source_node_archive_tolerance.py` (new)
- `annotation/tests/test_annotation_queryset_archive.py` (new)
- `upload/tests/test_restore_vcf.py` (new)
- `genes/tests/test_gene_coverage_archive.py` (new)

---

## 10. Risks and implementation-time checks

1. **Reuse, don't extract.** `VCF.delete_internal_data()` already does the cleanup; the only refactor is splitting the partition-recreate step out so archive can skip it. `archive_vcf` then does: zygosity decrement → `delete_internal_data(recreate_partitions=False)` → stamp mixin → save. No `VCF.delete()` override needed (Django cascade + existing `vcf_pre_delete_handler` cover hard-delete).
2. **GCC cleanup is already encapsulated.** `GeneCoverageCollection.delete_related_objects()` (inherited from `RelatedModelsPartitionModel`) drops the partition tables. Archive calls it directly without deleting the GCC row. No extraction needed.
3. **Chokepoint blast radius (Q4).** Raising `DataArchivedError` from `Cohort.cohort_genotype_collection` propagates to many call sites — that's intentional. During implementation, walk every caller (`Sample.cohort_genotype_collection`, `Sample.get_genotype`, `Sample.get_variant_qs`, `VCF.get_variant_qs`, `genotype_vcf_processor_factory`, etc.) and confirm they either propagate cleanly or are already inside the source-node helper / detail-page banner that catches and surfaces the message. Anywhere a callee swallows the exception silently, surface "archived" through the existing UI channel instead.
4. **Sub-cohort badge UX (Q8).** Sub-cohorts inherit the parent's `cohort_genotype_collection`, so the chokepoint catches them automatically. Cohort listing/detail need to wrap the access in try/except to render the badge — no separate cascade flag needed.
5. **Node-version bump path.** Archive must invalidate the `AnalysisNode` Q-object cache (Redis) and `CompHet` `cache_memoize`. `delete_internal_data()` already runs as part of the existing reload flow without an explicit bump, so verify whether reload-side cache invalidation is implicit (e.g. via cohort version increment) before adding bump logic. If not, plumb an invalidation call into `archive_vcf` after the data drop.
6. **Cross-ticket boundary.** #1536 owns `DataArchiveMixin`, the path-resolution helper (`get_annotation_archive_path`), and the queryset-builder guards. #1537 consumes the helper to write/read partition dumps and wires VAV archive actions + eventlog events. Document this in both tickets.
7. **Restore path verification.** During implementation, test that `retry_upload_pipeline` → `reload_vcf_task` preserves `vcf_id` (via `ImportCreateVCFModelForGenotypeVCFTask.process_items:42`). Spec assumes this works; verify with a test before relying on it.
8. **No `is_restore` flag (Q6).** `reload_auto_analyses_for_vcf` (`analysis/tasks/auto_analysis_tasks.py:28`) is the existing reload-side hook that refreshes node counts without auto-launching, and `auto_run_analyses_for_vcf` already takes `skip_already_analysed`. Verify during implementation that `vcf_import_success_signal` consumers behave correctly for a VCF that already has analyses; if not, surface the issue and add a flag at that point rather than threading one preemptively.
9. **GCC `data_restorable_from` (Q9).** GCC restore re-imports through the existing pipeline (no partition dump). `data_restorable_from` records `gcc.path`. If the path can be either backend-mounted or uploaded-via-webapp, the existing import pipeline already handles both — no special branching needed in archive.
10. **Top-level imports.** Mixin module is safe (no app deps). `snpdb/archive.py` calls `eventlog.create_event(...)` (snpdb→eventlog — confirm this is an existing dependency direction or use lazy import). `analysis` imports `snpdb` already, so source-node helpers using the chokepoint exception are clean.
11. **Bootstrap 4 only** — `data-toggle`, never `data-bs-toggle`.
12. **DataTables vs jqGrid** — keep existing tables, add badge only.

---

## 11. Recommended deferrals

- Per-sample archive (separate axis from VCF; would require repacking CGC arrays).
- Custom (non-VCF) cohort archive.
- VAV archive action + admin wiring (lands with #1537).
- GCC detail-page UI / Guardian permissions overhaul (Q2 — separate issue if access surfaces as a problem).
- jqGrid → DataTables migration for VCF list pages.
- Auditlog scope expansion (eventlog covers this ticket's needs).

---

## 12. Cross-ticket contracts

**This ticket exposes for #1537:**
- `DataArchiveMixin` on VAV (fields available for stamping). There is no separate VTAV class — `VariantTranscriptAnnotation` and `VariantGeneOverlap` are partitioned under VAV.
- `get_annotation_archive_path(partition_filename)` helper.
- `ANNOTATION_ARCHIVE_DIR` setting.
- Annotation queryset builders already guard on `vav.data_archived` — #1537's archive action only needs to stamp + drop partition + write dump file.
- **Per-VAV partition dump set is three tables**: `VariantAnnotation`, `VariantTranscriptAnnotation`, `VariantGeneOverlap`. All three FK to VAV with `on_delete=CASCADE` and must be included in #1537's pg_dump and restore I/O.

**#1537 owes back:**
- Wiring the `archive_vav` action (admin action + eventlog event).
- pg_dump file format + sha256 storage convention.
- Restore-from-dump implementation reading from `data_restorable_from`.
