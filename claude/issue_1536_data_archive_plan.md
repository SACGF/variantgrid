# Plan: SACGF/variantgrid #1536 — `DataArchiveMixin` + audit + VCF/GeneCoverage re-import

Step 2 of the umbrella ticket #1539. #1537 (partition pg_dump pipeline) is a separate ticket that consumes the contracts defined here.

## Scope summary

- Add `DataArchiveMixin` (abstract) with four fields to: **VCF**, **VariantAnnotationVersion**, **VariantTranscriptAnnotationVersion**, **GeneCoverageCollection**.
- Wire archive/restore actions for **VCF** and **GeneCoverageCollection** (admin-only). VAV/VTAV get the mixin only — no archive action in this ticket; #1537 wires that.
- Single chokepoint for archive-aware queryset behaviour (cohort-genotype-collection level + annotation queryset builders).
- Detail-page banners (read-only) and listing badges on existing grids.
- Eventlog events on archive/restore.

**Out of scope:** per-sample archive, custom (non-VCF) cohort archive, Sample/Cohort/Trio mixin, #1537 partition dump pipeline, #1534 declarative partitions.

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
- `annotation/models/models.py:???` — VariantTranscriptAnnotationVersion (locate alongside VAV)
- `genes/models/...` — GeneCoverageCollection (locate exact module during implementation)

### Migrations (additive, three files)

- `snpdb/migrations/0XXX_vcf_data_archive.py` — `AddField` × 4 on `snpdb_vcf`
- `annotation/migrations/0XXX_data_archive.py` — `AddField` × 4 on both `annotation_variantannotationversion` and `annotation_varianttranscriptannotationversion`
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

The original "edit every source node" design is overkill. Source nodes (Sample/Trio/Cohort/Pedigree) all funnel through `Cohort.cohort_genotype_collection`. Intercept there.

### Primary chokepoint: `Cohort.cohort_genotype_collection` (`snpdb/models/models_cohort.py`)

Raise a new `DataArchivedError` if `self.vcf and self.vcf.data_archived`. (Custom cohorts have `cohort.vcf is None` → not archivable; no special handling per spec.)

### Source-node shared helper

In `analysis/models/nodes/sources/` — small mixin or helper used by `SampleNode`, `TrioNode`, `CohortNode`, `PedigreeNode`:

- `_get_node_arg_q_dict` wrapper: catch `DataArchivedError` → return `{None: {str(q_none): q_none}}` (mirrors the existing "no sample" path)
- `_get_configuration_errors` wrapper: surface `f"Source data archived {date} by {user}: {reason}"` → drives `NodeStatus.ERROR_CONFIGURATION` via existing `get_status_from_errors` machinery (`analysis_node.py:723`)

`AllVariantsNode` is unaffected — not pinned to VCF/Sample/Cohort/Trio.

### Annotation queryset builders

| Code path | File:line | Planned |
|---|---|---|
| `get_variant_queryset_for_annotation_version` | `annotation/annotation_version_querysets.py:31` | Guard on VAV archived → `qs.none()` + `report_message` (from `library.log_utils`) |
| `get_variants_qs_for_annotation` | `:52` | Same guard |
| VTAV equivalent queryset path | (locate during implementation) | Same `qs.none()` + `report_message` treatment |
| `VariantTranscriptAnnotation.get_overlapping_genes_q` (queries `VariantGeneOverlap`) | `annotation/models/models.py:1630–1638` | Filters `VariantGeneOverlap.objects.filter(version=variant_annotation_version, ...)` — VAV-partitioned. Returns empty Q after archive but silently. | If `variant_annotation_version.data_archived` → return `Q(pk__in=[])` + `report_message`. Used by gene-list node filtering. |

### `Analysis.get_errors` (consumed by `_get_analysis_errors` at `analysis_node.py:694`)

Add archive check: if pinned `AnnotationVersion`'s VAV **or** VTAV is archived → analysis-level error. Maps to `NodeStatus.ERROR_WITH_PARENT` (`analysis_node.py:725`). One change propagates to all node types.

`Analysis.get_warnings` (`analysis/models/models_analysis.py:168`): also surface archived VAV/VTAV as warning so listing/detail UI shows status without re-running.

### `GeneCoverageCollection` chokepoint

Identify the read entry point during implementation (likely `GeneCoverageCollection.get_gene_coverage()` or a queryset builder). Add guard: `qs.none()` + `report_message` if `collection.data_archived`. Single intercept, same pattern.

### Other independent paths (not under the chokepoint)

| Code path | File:line | Planned |
|---|---|---|
| `variant_sample_information.show_variant_sample_information` | `snpdb/variant_sample_information.py:85–120` | Filter out rows where joined VCF archived. Joins `variant→cohortgenotype→...→vcf` directly, doesn't go through `cohort_genotype_collection`. **Critical** — only place not covered by the chokepoint. |
| `VCFAlleleSource.get_variant_qs` | `snpdb/models/models_vcf.py:587–598` | `qs.none()` when `vcf.data_archived` — prevents LiftoverRun on archived data |
| `VCF.filter_for_user` | `snpdb/models/models_vcf.py:120` | Add `include_archived: bool = True`. Default keeps detail pages working; autocompletes/source-node selectors pass `False` so archived rows aren't pickable for new analyses. Pinned existing analyses still resolve. |
| Trio/Sample/Cohort autocompletes | `snpdb/views/views_autocomplete.py` | Thread `include_archived=False` via the inner VCF lookup |
| `auto_run_analyses_for_vcf` / `_for_sample` | `analysis/tasks/auto_analysis_tasks.py:9,18` | Defensive early return when archived. Also consumed by the restore "no auto-launch" requirement (Q8) — see §4. |

### Cascades that don't need touching (called out explicitly)

- `vcf_pre_delete_handler` and `pre_delete_cohort` only fire on row delete. Archive does **not** delete the row (see §4). Whatever data-cleanup these signals currently do must be moved into the extracted `_delete_underlying_data()` method so both `delete()` and `archive_vcf()` get the cleanup. Signals can stay as thin wrappers calling that method, or the signals can be removed if the override approach is cleaner.
- Guardian permissions unchanged.

---

## 4. Archive + restore actions

### Archive

The whole point of archive is **keep the metadata, drop the big underlying data**. Calling `vcf.delete()` would destroy the VCF row and its mixin stamp along with it — wrong.

**Refactoring step (precondition for archive):** extract the existing data-deletion logic out of `VCF.delete()` (or its `pre_delete` signal handler at `snpdb/models/models_vcf.py:223`) into a standalone method `VCF._delete_underlying_data()` that drops:
- `CohortGenotype*` partition rows for this VCF
- Decrements global zygosity counts (`update_all_variant_zygosity_counts_for_vcf(vcf, '-')`)
- Fires the analysis node-version-bump path (variantgrid_com#22 pattern)
- Any other cascade-cleanup the existing delete pre/post-signals do

Both `VCF.delete()` (existing behaviour: drop data **then** delete row) and `archive_vcf()` (new behaviour: drop data, **keep** row, stamp mixin) call the extracted method.

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
    vcf._delete_underlying_data()  # extracted from VCF.delete(): partitions, zygosity counts, node-version bumps
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

`VCF.delete()` becomes:

```python
def delete(self, *args, **kwargs):
    self._delete_underlying_data()
    super().delete(*args, **kwargs)
```

Equivalent extraction for GeneCoverageCollection: `GeneCoverageCollection._delete_underlying_data()` drops `GeneCoverage` rows (and any other heavy cascade data); both `delete()` and `archive_gene_coverage_collection()` call it.

Parallel `archive_gene_coverage_collection(gcc, user, reason)` — same pattern. Locate the GCC source link (`UploadedFile` FK or backend path field) during implementation; capture sha256 if not already recorded (compute at archive time and store via the mixin if needed — flag).

VAV/VTAV: **no archive helper in this ticket.** Mixin lands; #1537 wires the action.

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
    return retry_upload_pipeline(upload_pipeline, is_restore=True)
```

Reuses `retry_upload_pipeline` (`upload/uploaded_file_type.py:90`) → routes to `reload_vcf_task` (`upload/tasks/vcf/genotype_vcf_tasks.py:193`). The "import into existing VCF.id" mode the spec wants **already exists** at `genotype_vcf_tasks.py:42` — `ImportCreateVCFModelForGenotypeVCFTask.process_items` reuses `upload_pipeline.uploadedvcf.vcf` when set, calling `configure_vcf_from_header` to refresh format fields.

Clear `data_archived_date` inside the existing success path (`ImportGenotypeVCFSuccessTask.process_items`), not in the entry helper. Emit `vcf_restored` eventlog event there.

**No auto-launch on restore (Q8):** thread `is_restore=True` through `retry_upload_pipeline` → `reload_vcf_task` → success handler. Success handler skips `auto_run_analyses_for_vcf` when the flag is set. (Alternative: rely on the existing logic if it already no-ops for VCFs with prior analyses — verify during implementation; if so, no flag needed.)

Parallel `restore_gene_coverage_collection(gcc, user)` — locate the existing reload pipeline for gene coverage (likely `genes/tasks/` or `upload/tasks/`) and reuse.

---

## 5. Permissions and entry points (Q6 — revised)

### VCF and GeneCoverageCollection — write-access users via the GUI

Archive/restore for **VCF** and **GeneCoverageCollection** is available to any user with `can_write` on the object (Guardian permission, via `GuardianPermissionsMixin`). Buttons live on the detail page, not in Django admin.

**Detail-page archive/restore buttons** on `view_vcf.html` and the GeneCoverageCollection detail template:
- Show "Archive" button when `obj.can_write(user) and not obj.data_archived` and the precondition check below passes.
- Show "Restore" button when `obj.can_write(user) and obj.data_archived` and the recorded `data_restorable_from` exists on disk.
- Confirmation modal (Bootstrap 4 `data-toggle="modal"`) prompts for `data_archive_reason` on archive.

**Pre-archive precondition check.** Before showing the archive button (and re-checked server-side in the view handler):
- Resolve the uploaded source path: for VCF, `vcf.uploadedvcf.uploaded_file.get_filename()`; for GCC, the equivalent source link.
- If `os.path.exists(path)` is False → button is disabled / replaced with a message:
  > **Uploaded file doesn't exist — cannot archive.** If you want to free up space you can permanently delete this VCF.
  Render a "Delete permanently" button next to the message (links to the existing delete flow that already exists for write-access users).
- If `vcf.uploadedvcf` doesn't exist at all (legacy VCFs imported without an uploaded-file record) → same treatment: archive disabled, permanent-delete offered.

This is the gate that keeps "archive" honest — archive is only a meaningful operation if restore is possible. When the source has gone missing, the right answer is permanent delete, not archive-and-pray.

**Views** (`snpdb/views/views.py`):
- `archive_vcf_view(request, vcf_id)` — POST, `obj.check_can_write(request.user)`, runs the precondition check, calls `archive_vcf(...)`, redirects to detail page.
- `restore_vcf_view(request, vcf_id)` — POST, `obj.check_can_write(request.user)`, calls `restore_vcf(...)`, redirects to detail page.
- Parallel `archive_gene_coverage_collection_view` / `restore_gene_coverage_collection_view`.

URL routes added to `snpdb/urls/urls.py` (and `genes/urls.py`).

### VAV / VTAV — admin only (no action wired in this ticket)

VAV and VTAV are global system-level objects, not user-owned. Archive action lands with #1537 — when it does, it goes through Django admin (superuser-only), since regular write users have no business archiving annotation versions.

For this ticket, the mixin lands and admin list views show the four mixin fields. No action wired.

### Django admin (still useful for sysadmin oversight)

- `snpdb/admin.py` — `VCFAdmin` with `list_display` and `list_filter` covering the mixin fields. No archive action — sysadmins use the GUI like everyone else.
- `annotation/admin.py` — `VariantAnnotationVersionAdmin`, `VariantTranscriptAnnotationVersionAdmin` (mixin fields visible, no actions).
- `genes/admin.py` — `GeneCoverageCollectionAdmin` with `list_display` covering mixin fields.

---

## 6. Eventlog integration (Q7)

Use the project's existing `eventlog` app (per `claude/research/eventlog.md`). Locate the event-creation pattern during implementation (`Event.objects.create(...)` or a `create_event(...)` helper).

Events emitted from inside the archive/restore helpers (post-success, not via model `save()` so emits tie to deliberate operator actions):

- `vcf_archived` — `{vcf_id, vcf_name, reason, restorable_from}`, actor = archive user
- `vcf_restored` — `{vcf_id, restored_from}`, actor = restore user
- `gene_coverage_collection_archived` / `_restored` — same shape

VAV/VTAV: no events in this ticket; #1537 wires them when the archive action lands.

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

### Listing badges (Q11 — keep existing tables)

- VCF list (jqGrid in `snpdb/grids.py`): add an "Archived" badge via cell formatter on an existing column (or new read-only column). No replacement.
- GeneCoverageCollection list: same pattern — locate the existing grid and add a badge.
- Analyses listing (`analysis/grids.py:911` `AnalysesColumns`, already DataTables): add a "Source data archived" badge column annotated via `Exists(...)` subquery walking `analysisnode_set` for archived FK targets. Add a "Hide analyses with archived data" filter wired into `get_q_list` (`:942`). Detail-side analysis editor surfaces archive via `Analysis.get_warnings` for free.

---

## 8. Tests

All `django.test.TestCase`, run with `python3 manage.py test --keepdb`.

### `snpdb/tests/test_data_archive_mixin.py` (new)

- Mixin fields present on VCF, VAV, VTAV, GCC; `data_archived` property reflects `data_archived_date`.
- `archive_vcf` stamps fields, populates `data_restorable_from` from upload path, calls `vcf.delete()` path.
- After archive, `cohortgenotype_set` empty (or partition reflects existing delete behaviour).
- `archive_vcf` is idempotent — second call no-ops.
- Eventlog event emitted with correct payload.

### `snpdb/tests/test_filter_for_user_archive.py`

- `filter_for_user(user)` includes archived by default (detail page works).
- `filter_for_user(user, include_archived=False)` excludes archived.

### `analysis/tests/test_source_node_archive_tolerance.py`

- For each of `SampleNode`, `TrioNode`, `CohortNode`, `PedigreeNode`: archive the source's VCF, assert configuration error mentioning "archived", `_get_node_arg_q_dict` returns `q_none`, queryset returns 0 rows without raising.
- `Analysis.get_errors` surfaces VAV-archived as analysis-level error.
- `Analysis.get_errors` surfaces VTAV-archived as analysis-level error.

### `annotation/tests/test_annotation_queryset_archive.py`

- `get_variant_queryset_for_annotation_version` returns `qs.none()` and emits `report_message` when VAV archived.
- Same for VTAV equivalent.
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
- Calls `retry_upload_pipeline` with `is_restore=True`; success handler skips `auto_run_analyses_for_vcf`.
- Restore button hidden when user lacks `can_write`.

### `genes/tests/test_gene_coverage_archive.py`

- Archive stamps fields, calls delete, eventlog event emitted.
- Restore validates path, calls reload pipeline.
- Queryset chokepoint returns empty when GCC archived.

---

## 9. File-by-file change list

### Mixin / migrations
- `library/django_utils/data_archive_mixin.py` (new)
- `snpdb/models/models_vcf.py` — VCF inherits mixin
- `annotation/models/models.py` — VAV + VTAV inherit mixin
- `genes/models/...` — GeneCoverageCollection inherits mixin (locate file)
- `snpdb/migrations/0XXX_vcf_data_archive.py` (new)
- `annotation/migrations/0XXX_data_archive.py` (new)
- `genes/migrations/0XXX_genecoveragecollection_data_archive.py` (new)

### Settings + path helper
- `variantgrid/settings/components/annotation_settings.py` — `ANNOTATION_ARCHIVE_DIR`
- `annotation/archive_paths.py` (new) — `get_annotation_archive_path()`

### Audit (chokepoint + independent paths)
- `snpdb/models/models_cohort.py` — `Cohort.cohort_genotype_collection` raises `DataArchivedError`
- `analysis/models/nodes/sources/` — shared helper for `SampleNode`/`TrioNode`/`CohortNode`/`PedigreeNode` to catch + surface error
- `annotation/annotation_version_querysets.py:31,52` — VAV archive guards
- VTAV equivalent queryset builder (locate)
- `annotation/models/models.py:1630` — `VariantTranscriptAnnotation.get_overlapping_genes_q` VAV archive guard (covers `VariantGeneOverlap` queries)
- `analysis/models/models_analysis.py:168` — `Analysis.get_warnings` adds VAV/VTAV warnings
- `analysis/models/models_analysis.py` — `Analysis.get_errors` adds VAV/VTAV errors
- `snpdb/variant_sample_information.py:85–120` — filter archived VCFs
- `snpdb/models/models_vcf.py:587–598` — `VCFAlleleSource.get_variant_qs` archived-check
- `snpdb/models/models_vcf.py:120` — `VCF.filter_for_user(include_archived=True)`
- `snpdb/views/views_autocomplete.py` — pass `include_archived=False`
- `analysis/tasks/auto_analysis_tasks.py:9,18` — defensive guards
- `genes/...` — `GeneCoverageCollection` chokepoint guard (locate read entry point)

### Archive/restore helpers + views
- `snpdb/archive.py` (new) — `archive_vcf`, `restore_vcf`, `ArchivePreconditionError`
- `genes/archive.py` (new) — `archive_gene_coverage_collection`, `restore_gene_coverage_collection`
- `snpdb/views/views.py` — `archive_vcf_view`, `restore_vcf_view` (Guardian `check_can_write` gated)
- `genes/views.py` — `archive_gene_coverage_collection_view`, `restore_gene_coverage_collection_view`
- `snpdb/urls/urls.py`, `genes/urls.py` — URL routes
- `upload/tasks/vcf/genotype_vcf_tasks.py` — thread `is_restore` flag through `reload_vcf_task` and success handler

### Admin (oversight only — no archive actions)
- `snpdb/admin.py` — `VCFAdmin` with mixin fields in `list_display` / `list_filter`
- `annotation/admin.py` — VAV / VTAV admin (mixin fields visible, no actions until #1537)
- `genes/admin.py` — `GeneCoverageCollectionAdmin` with mixin fields in list_display

### Templates
- `snpdb/templates/snpdb/data/_data_archived_banner.html` (new — archived banner + restore button)
- `snpdb/templates/snpdb/data/_archive_action.html` (new — archive button + "uploaded file missing" message + permanent-delete fallback)
- `snpdb/templates/snpdb/data/view_vcf.html` — include both partials
- GeneCoverageCollection detail template — include both partials (locate)
- `snpdb/grids.py` — VCF jqGrid archived badge
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

1. **Extract `VCF._delete_underlying_data()` cleanly.** The existing data-cleanup is split across `VCF.delete()` (or its overrides), `vcf_pre_delete_handler` (`snpdb/models/models_vcf.py:223`), and possibly Celery tasks (`snpdb/tasks/soft_delete_tasks.py`). Auditing every site and pulling the partition-drop / zygosity-decrement / node-version-bump steps into one method is the load-bearing refactor of this ticket. Risk: some cleanup may rely on the row already being gone (e.g. cascade FK deletes). Where that's true, archive must implement equivalent SQL-level cleanup that keeps the VCF row intact.
2. **GCC `_delete_underlying_data()` extraction.** Same refactor for `GeneCoverageCollection`. Locate where `GeneCoverage` rows are currently dropped on GCC delete and pull into a standalone method.
3. **Custom (non-VCF) cohort handling.** `Cohort.vcf` is nullable. The chokepoint guard `cohort.vcf and cohort.vcf.data_archived` correctly no-ops for these — confirm during implementation that no other path treats custom cohorts as archivable.
4. **Node-version bump path.** Archive must invalidate the `AnalysisNode` Q-object cache (Redis) and `CompHet` `cache_memoize`. Per Q4, reuse the variantgrid_com#22 sample-delete bump path. The bump must live inside `_delete_underlying_data()` so both delete and archive paths trigger it — if it currently sits in a `pre_delete` signal, move it into the extracted method.
5. **Cross-ticket boundary.** #1536 owns `DataArchiveMixin`, the path-resolution helper (`get_annotation_archive_path`), and the queryset-builder guards. #1537 consumes the helper to write/read partition dumps and wires VAV/VTAV archive actions + eventlog events. Document this in both tickets.
6. **Restore path verification.** During implementation, test that `retry_upload_pipeline` → `reload_vcf_task` actually preserves `vcf_id` (via `ImportCreateVCFModelForGenotypeVCFTask.process_items:42`). Spec assumes this works; verify with a test before relying on it.
7. **`is_restore` flag plumbing.** Threading the flag through `retry_upload_pipeline` → `reload_vcf_task` → success handler may require signature changes on multiple Celery tasks. If the existing `auto_run_analyses_for_vcf` already no-ops when prior analyses exist for a VCF, the flag is unnecessary — verify first.
8. **GCC source-link shape.** "Either uploaded files or backend ones" per spec. Locate the actual link (FK or path field) before deciding `data_restorable_from` semantics. If backend-generated without a stored sha256, archive-time hash computation is needed.
9. **Top-level imports.** Mixin module is safe (no app deps). `snpdb/archive.py` calls `eventlog.create_event(...)` (snpdb→eventlog — confirm this is an existing dependency direction or use lazy import). `analysis` imports `snpdb` already, so source-node helpers using the chokepoint exception are clean.
10. **Bootstrap 4 only** — `data-toggle`, never `data-bs-toggle`.
11. **DataTables vs jqGrid** — keep existing tables, add badge only.

---

## 11. Recommended deferrals

- Per-sample archive (separate axis from VCF; would require repacking CGC arrays).
- Custom (non-VCF) cohort archive.
- VAV/VTAV archive actions + admin wiring (lands with #1537).
- jqGrid → DataTables migration for VCF list pages.
- Auditlog scope expansion (eventlog covers this ticket's needs).

---

## 12. Cross-ticket contracts

**This ticket exposes for #1537:**
- `DataArchiveMixin` on VAV and VTAV (fields available for stamping).
- `get_annotation_archive_path(partition_filename)` helper.
- `ANNOTATION_ARCHIVE_DIR` setting.
- Annotation queryset builders already guard on archive — #1537's archive action only needs to stamp + drop partition + write dump file.
- **Per-VAV partition dump set is three tables**: `VariantAnnotation`, `VariantTranscriptAnnotation`, `VariantGeneOverlap`. All three FK to VAV with `on_delete=CASCADE` and must be included in #1537's pg_dump and restore I/O.

**#1537 owes back:**
- Wiring the `archive_vav` / `archive_vtav` actions (admin actions + eventlog events).
- pg_dump file format + sha256 storage convention.
- Restore-from-dump implementation reading from `data_restorable_from`.
