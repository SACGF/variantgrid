# Issue #1551 — Pre-compute non-missing VariantCollection per sub-cohort

Reference: https://github.com/SACGF/variantgrid/issues/1551 (follow-up to #1546)

## Goal

Replace the wide-cohort EXCLUDE regex seq-scan in
`analysis/models/nodes/cohort_mixin.py` with a hash join against a
pre-computed `VariantCollection` of "variants where at least one sub-cohort
sample is non-missing." Build the VC once at sub-cohort finalisation; the
runtime filter joins to it.

The regex path stays as a correctness fallback for when the VC is not yet
ready, is stale, or build failed — no UI surprises during the build window.

## Trigger surface (clarified by user)

Sub-cohorts are created in two places, both user-driven:

1. **VCF page → "create sub cohort"** (`snpdb/views/views_json.py:66`,
   `create_sub_cohort` view). One-shot: sample list is fixed at creation
   time, no later editing from this surface.
2. **Cohort page → edit samples → "Save cohort"**
   (`snpdb/views/views_json.py:50`, `create_cohort_genotype` view → calls
   `create_cohort_genotype_and_launch_task`). Sample membership is mutated
   via `CohortSample.save()`/`delete()` (which calls `increment_version()`),
   but the user explicitly finalises with the save button.

Both surfaces are explicit, single-shot user actions, so we enqueue the
VC build directly at those two call sites — no debouncing of
`increment_version()` needed, and no signal on `CohortSample`.

## Files to read

- `claude/research/snpdb.md` (research overview of cohorts/VCs).
- `analysis/models/nodes/cohort_mixin.py` — the filter site (§3 below).
- `snpdb/models/models_cohort.py` — `Cohort.create_sub_cohort`,
  `Cohort.increment_version`, `pre_delete_cohort`,
  `CohortGenotypeCollection`, `CohortSample`.
- `snpdb/models/models_variant.py` lines 865-902 — `VariantCollection`,
  `VariantCollectionRecord`, `get_annotation_kwargs`,
  `variant_collection_alias`.
- `snpdb/tasks/cohort_genotype_tasks.py` —
  `create_cohort_genotype_and_launch_task` and `cohort_genotype_task`
  shape/queue patterns.
- `snpdb/tasks/vcf_bed_file_task.py` — concrete VC build pattern (creates
  VC, drives partition fill, sets status).
- `snpdb/views/views_json.py` lines 40-77 — the two trigger views.
- `analysis/management/commands/profile_analysis_nodes.py` lines 641-728 —
  `_profile_cohort_exclude_vc_join` is the SQL-shape reference for the
  build INSERT and the JOIN-time annotation; mirror this in production.

## Changes

### §1 — `CohortVersion` + link model, mirroring `NodeVersion`

Adopt the `NodeVersion` pattern from
`analysis/models/nodes/analysis_node.py:1167` — a per-version row that
version-keyed caches FK to with `on_delete=CASCADE`, so cleanup happens
by deleting stale version rows, not by matching tuples at query time.

**New model `CohortVersion`** in `snpdb/models/models_cohort.py`:

```python
class CohortVersion(TimeStampedModel):
    """ One row per (cohort, version). Version-specific caches FK here
    with on_delete=CASCADE; periodic cleanup of stale rows drops every
    artefact tied to that version. Mirrors analysis.NodeVersion. """
    cohort = models.ForeignKey(Cohort, on_delete=CASCADE)
    version = models.IntegerField(null=False)

    class Meta:
        unique_together = ("cohort", "version")

    @staticmethod
    def get(cohort: Cohort) -> "CohortVersion":
        return CohortVersion.objects.get(cohort=cohort, version=cohort.version)

    def __str__(self):
        return f"CohortVersion({self.cohort.pk} v{self.version})"
```

Wire `Cohort.increment_version()`
(`snpdb/models/models_cohort.py:105-133`) to `get_or_create` the row for
the new version, mirroring `AnalysisNode.save`'s line 1031:

```python
def increment_version(self):
    # ... existing body ...
    self.delete_old_counts()
    CohortVersion.objects.get_or_create(cohort=self, version=self.version)
```

Also call `get_or_create` once in `Cohort.save()` for fresh cohorts (so
version 0 has a row before any bump).

**New model `SubCohortVariantCollection`** in the same file:

```python
class SubCohortVariantCollection(models.Model):
    """ Cached non-missing VariantCollection for a sub-cohort.

    Two staleness axes, each handled by CASCADE:
      - cohort_version: deletes when sub-cohort version row is GC'd
        (Cohort.increment_version creates a new CohortVersion; the
        sweeper task in §5 drops superseded rows).
      - parent_cohort_genotype_collection: deletes when the parent CGC
        is dropped (parent VCF re-import / parent cohort delete).
    """
    cohort_version = models.OneToOneField(
        CohortVersion, on_delete=CASCADE,
        related_name="sub_cohort_variant_collection")
    parent_cohort_genotype_collection = models.ForeignKey(
        CohortGenotypeCollection, on_delete=CASCADE)
    variant_collection = models.OneToOneField(
        "snpdb.VariantCollection", on_delete=CASCADE)
```

`OneToOneField(CohortVersion)` because at most one non-missing VC is
useful per (cohort, version, parent CGC); the parent CGC FK is implied
to be the current one at build time but doesn't need to be in the
unique key — if it changes, that side cascades anyway.

**`post_delete` signal** to drop the orphaned `VariantCollection`
partition when a link row is removed (mirrors `post_delete_node_cache`
at `analysis_node.py:1206-1215`):

```python
@receiver(post_delete, sender=SubCohortVariantCollection)
def post_delete_sub_cohort_variant_collection(sender, instance, **kwargs):
    try:
        if instance.variant_collection_id:
            instance.variant_collection.delete_related_objects()
            instance.variant_collection.delete()
    except VariantCollection.DoesNotExist:
        pass
```

Note: this handler fires on CASCADE deletion of the link as well as
direct deletes — covering both axes.

**Helper on `Cohort`** for the readiness lookup:

```python
def get_non_missing_variant_collection(self) -> Optional["VariantCollection"]:
    """ Returns the SUCCESS-status VC for the current (cohort_version,
    parent CGC); else None. """
    if not self.is_sub_cohort:
        return None
    try:
        cgc = self.cohort_genotype_collection
    except (CohortGenotypeCollection.DoesNotExist, DataArchivedError):
        return None
    link = (SubCohortVariantCollection.objects
            .filter(cohort_version__cohort=self,
                    cohort_version__version=self.version,
                    parent_cohort_genotype_collection=cgc,
                    variant_collection__status=ProcessingStatus.SUCCESS)
            .select_related("variant_collection")
            .first())
    return link.variant_collection if link else None
```

Migration:

1. Create `CohortVersion` table.
2. Data migration: `bulk_create` a `CohortVersion(cohort=c, version=c.version)` for every existing `Cohort` so the post-deploy state has rows for current versions.
3. Create `SubCohortVariantCollection` table.

No fields are added to `Cohort`, `CohortGenotypeCollection`, or
`VariantCollection` — they stay generic.

**Out of scope, worth noting in the PR description**: `CohortGenotypeCollection`
itself currently stores `cohort` + `cohort_version` directly as
columns. A future refactor could migrate CGC to FK `CohortVersion`
instead, unifying staleness handling. Out of scope here — the data
migration would be substantial — but `CohortVersion` is the foundation
that future change would build on.

### §2 — New task: `build_sub_cohort_non_missing_vc_task`

New file `snpdb/tasks/sub_cohort_tasks.py` (or add to
`snpdb/tasks/cohort_genotype_tasks.py` if cohesion fits).

Public helper:

```python
def enqueue_sub_cohort_non_missing_vc(cohort: Cohort) -> Optional[str]:
    """ Sole entry point for triggers in §4 to schedule the build.
    Returns the celery task id, or None if the cohort is not a sub-cohort
    or has no samples. Idempotent: drops any pre-existing VC + queued
    task for this sub-cohort before enqueuing."""
```

Task body (`@celery.shared_task` on `annotation_workers` queue, mirroring
`calculate_vcf_stats`):

1. Resolve `cohort = Cohort.objects.get(pk=cohort_pk)`. Skip if not a
   sub-cohort or sample set empty.
2. Resolve parent CGC via `cohort.cohort_genotype_collection` (this
   already follows `get_base_cohort()`). Resolve / create the
   `CohortVersion` row for `cohort.version`.
3. Drop any pre-existing `SubCohortVariantCollection` for this
   `CohortVersion` (idempotent re-build): `link.delete()` — the
   `post_delete` signal drops the VC partition.
4. `vc = VariantCollection.objects.create(name=f"sub_cohort_{cohort.pk}_v{cohort.version}_non_missing", status=ProcessingStatus.CREATED)`.
5. Build the regex used to detect "any sample non-missing for the
   sub-cohort positions" via the existing
   `cgc.get_sample_zygosity_regex(...)` call — same call shape that the
   filter currently uses, but with `exclude=False` semantics. Confirm the
   regex matches "at least one sub-cohort sample is non-missing" by
   reusing the regex builder (see §3 — we want the *positive* regex, not
   the negated lookahead).
6. INSERT FROM SELECT, mirroring `profile_analysis_nodes._profile_cohort_exclude_vc_join`:

   ```sql
   INSERT INTO <vc_partition>
       (variant_collection_id, variant_id)
   SELECT DISTINCT %s, variant_id
   FROM <uncommon_partition>
   WHERE collection_id = ANY(%s) AND samples_zygosity ~ %s
   UNION
   SELECT DISTINCT %s, variant_id
   FROM <common_partition>
   WHERE collection_id = ANY(%s) AND samples_zygosity ~ %s;
   ```

   `collection_ids = [cgc.pk]` for uncommon, `[cgc.common_collection_id]`
   for common. The `UNION` block for common is only emitted if
   `cgc.common_collection_id` is set.
7. Set `vc.count = cur.rowcount`, `vc.status = ProcessingStatus.SUCCESS`,
   `vc.save()`.
8. Create the link row:
   `SubCohortVariantCollection.objects.create(cohort_version=cohort_version, parent_cohort_genotype_collection=cgc, variant_collection=vc)`.
   No `Cohort.save()` is needed — staleness is encoded entirely in the
   `CohortVersion` FK — so we never accidentally bump `cohort.version`.
9. On exception: drop the partial VC (`delete_related_objects()` + `delete()`),
   `log_traceback()`, re-raise. Filter-time fall-through covers the user.

Queue: `@app.task(queue='annotation_workers')`.

### §3 — Filter-time change in `cohort_mixin.py`

Replace `analysis/models/nodes/cohort_mixin.py:143-148`. Today:

```python
if cohort.is_sub_cohort:
    missing = [Zygosity.UNKNOWN_ZYGOSITY, Zygosity.MISSING]
    sample_zygosities_dict = {s: missing for s in cohort.get_samples()}
    q_sub = cgc.get_zygosity_q(sample_zygosities_dict, exclude=True)
    q_and.append(q_sub)
```

After:

```python
if cohort.is_sub_cohort:
    vc = cohort.get_non_missing_variant_collection()  # helper from §1
    if vc is not None:
        q_and.append(Q(**{f"{vc.variant_collection_alias}__isnull": False}))
        # annotation_kwargs registration handled in
        # _get_annotation_kwargs_for_node — see below
    else:
        missing = [Zygosity.UNKNOWN_ZYGOSITY, Zygosity.MISSING]
        sample_zygosities_dict = {s: missing for s in cohort.get_samples()}
        q_sub = cgc.get_zygosity_q(sample_zygosities_dict, exclude=True)
        q_and.append(q_sub)
```

`Q` over an alias requires the alias to be registered in
`_get_annotation_kwargs_for_node`. Update the existing override at
`cohort_mixin.py:52-56`:

```python
def _get_annotation_kwargs_for_node(self, **kwargs) -> dict:
    annotation_kwargs = super()._get_annotation_kwargs_for_node(**kwargs)
    if cgc := self.cohort_genotype_collection:
        annotation_kwargs.update(cgc.get_annotation_kwargs(**kwargs))
    cohort = self._get_cohort()
    if cohort and cohort.is_sub_cohort:
        if vc := cohort.get_non_missing_variant_collection():
            annotation_kwargs.update(vc.get_annotation_kwargs(**kwargs))
    return annotation_kwargs
```

Both call sites use the single helper `Cohort.get_non_missing_variant_collection()`
from §1 — readiness logic lives in one place.

`_get_cache_key` (`cohort_mixin.py:30-36`) already includes the CGC pk.
That covers parent VCF re-imports (new CGC pk → new cache key) but not
the moment a VC build *completes* under the same CGC. To invalidate the
analysis cache when the VC transitions from missing → ready, append the
VC pk (or 0 if missing) to the cache key:

```python
def _get_cache_key(self) -> str:
    cache_key = super()._get_cache_key()
    cgc_id = 0
    vc_id = 0
    if cgc := self.cohort_genotype_collection:
        cgc_id = cgc.pk
    cohort = self._get_cohort()
    if cohort and cohort.is_sub_cohort:
        if vc := cohort.get_non_missing_variant_collection():
            vc_id = vc.pk
    return "_".join((cache_key, str(cgc_id), str(vc_id)))
```

### §4 — Trigger sites

Both triggers call `enqueue_sub_cohort_non_missing_vc(cohort)` from §2.

**Trigger 1 — VCF page**: `Cohort.create_sub_cohort` in
`snpdb/models/models_cohort.py:230-243`. After all `CohortSample` rows
have been created (and `increment_version()` has fired per save), call
the helper before returning. Rationale: per the user's clarification,
this surface is one-shot — sample list is fixed and the user expects the
sub-cohort to be ready to query immediately.

**Trigger 2 — Cohort page**: `create_cohort_genotype_and_launch_task` in
`snpdb/tasks/cohort_genotype_tasks.py:20`. This function already handles
the case where the to-be-built cohort is a subset of an existing VCF
cohort (`containing_cohort` path) — it sets `parent_cohort` and
`import_status`. After that branch (and on the "already a sub-cohort"
case where containing_cohort matches existing parent), call the helper.

Concretely — in the `containing_cohort` branch (`models_cohort.py`-side
flow today sets parent_cohort + status, returns `("SUCCESS", None)`),
add `enqueue_sub_cohort_non_missing_vc(cohort)` before returning.

### §5 — Invalidation

All staleness propagates via FK CASCADE; deletion of the link row fires
the `post_delete` handler in §1 which drops the VC partition.

| Trigger | Action |
|---|---|
| `Cohort.create_sub_cohort()` | §4 trigger 1 enqueues build. Build task is idempotent (step 3 in §2). |
| `create_cohort_genotype_and_launch_task` for a sub-cohort | §4 trigger 2 enqueues build. |
| `CohortSample` save/delete on a sub-cohort | `Cohort.increment_version()` creates a new `CohortVersion`. The OLD `CohortVersion` row lingers (it has the old link FK pointing at it) until the periodic sweeper drops it — exactly like `NodeVersion`. Filter-time lookup queries by `cohort_version__version=self.version` so the old link is invisible. No code change needed. |
| Parent VCF re-import (new parent CGC) | The old parent CGC may be deleted (via `delete_old_counts` → `marked_for_deletion=True` → `delete_old_cohort_genotypes_task`). When that runs, `on_delete=CASCADE` on `SubCohortVariantCollection.parent_cohort_genotype_collection` drops the link → `post_delete` handler drops the VC. Lookup at filter time finds nothing → regex fallback until the user re-Saves. |
| `pre_delete_cohort` (parent VCF cohort delete) | The existing handler at `snpdb/models/models_cohort.py:301-308` resets `sub_cohort.import_status = CREATED`. Cascade chain: parent CGC delete → `SubCohortVariantCollection` link delete → `post_delete` → VC partition cleanup. No additional code in the handler. |

**Sweeper task** for stale `CohortVersion` rows, mirroring
`delete_analysis_old_node_versions` at
`analysis/tasks/node_update_tasks.py:143-150`:

```python
@celery.shared_task(ignore_result=False)
def delete_old_cohort_versions():
    """ Drops every CohortVersion that isn't the current version of its
    Cohort. CASCADE removes any version-keyed caches FK'd to it
    (currently SubCohortVariantCollection). """
    qs = CohortVersion.objects.all()
    latest = qs.filter(cohort_id=OuterRef('cohort_id')).order_by('-version')
    sub_query = Subquery(latest.values('pk')[:1])
    qs.annotate(latest_pk=sub_query).exclude(pk=F('latest_pk')).delete()
```

Run on the same trigger profile as the analysis sweeper — for our scope,
calling it from the build task right after `link.create()` is sufficient
(opportunistic cleanup; cheap because the qs is small per cohort). A
periodic Celery beat schedule is optional follow-up.

No signal on `CohortSample` save/delete: per the user, sub-cohort
membership only changes via the explicit user action that calls Trigger 2.

### §6 — Backfill (optional, recommended)

Add a one-shot management command
`snpdb/management/commands/build_sub_cohort_non_missing_vcs.py`:

- Iterate `Cohort.objects.filter(parent_cohort__isnull=False, import_status=ImportStatus.SUCCESS)`.
- For each, call `enqueue_sub_cohort_non_missing_vc` synchronously
  (`.apply()`) so the operator can monitor progress.
- Skip cohorts where `cohort.get_non_missing_variant_collection()`
  already returns non-None — already done.

This avoids forcing the first user of each existing sub-cohort to pay
the regex cost while waiting for the lazy-build chain that doesn't exist
(we don't auto-enqueue at filter time — keeping request paths
read-only).

### §7 — Tests

Add to `analysis/tests/test_zygosity_nodes.py` (already has a sub-cohort
fixture at line 22):

1. **VC built on create_sub_cohort**: after creation,
   `cohort.get_non_missing_variant_collection()` returns a VC with
   `status=SUCCESS`; the `CohortVersion` for `(cohort, cohort.version)`
   exists, and a `SubCohortVariantCollection` link FKs to it.
2. **VC contents correct**: VC contains exactly the variants where at
   least one sub-cohort sample has zygosity in {R, E, O, U}; reuse the
   existing fixture and assert against a known set.
3. **Filter uses VC join when ready**: the SQL produced by the cohort
   node references `variantcollection_<pk>` and not `samples_zygosity ~`.
4. **Fall-through when no link row**: delete the
   `SubCohortVariantCollection` row and assert the filter falls back to
   regex (SQL contains `samples_zygosity`).
5. **Stale on sub-cohort version bump**: add a `CohortSample`
   (bumping `cohort.version`), assert `get_non_missing_variant_collection()`
   returns None and filter falls through to regex.
6. **Cascade on parent CGC delete**: delete the parent CGC, assert the
   `SubCohortVariantCollection` row and its VC partition are gone (via
   CASCADE + `post_delete` handler).
7. **Sweeper drops stale CohortVersion + cascades**: bump version, run
   `delete_old_cohort_versions`, assert the old `CohortVersion` row,
   the link row, and the VC partition are all gone.
8. **pre_delete_cohort cascades correctly**: delete the parent cohort,
   assert the cascade chain leaves no orphaned VC partition.

For task-level coverage, the existing `cohort_genotype_task` tests are a
template — small synthetic cohort + assert row count post-build.

### §8 — Out of scope (deliberately)

- Sub-cohort grid panel het/hom/ref/unk counts (slow on wide cohorts via
  the substring fallback in `cohort_node.py:135-164`). Issue body calls
  this out as a clean follow-up; it can re-use the now-pre-computed VC
  to filter the queryset before running `calculate_cohort_stats`. New
  issue, not this PR.
- Multi-VCF custom cohorts (`vcf__isnull=True`,
  `parent_cohort__isnull=True`). They already build their own CGC and
  don't suffer the regex EXCLUDE problem.
- Caching VCs across sub-cohort versions. We drop and rebuild on every
  finalise; if that turns out to be expensive in practice, a content-hash
  reuse layer is a future optimisation.

## Expected impact (from issue body, prod cohort 659)

| metric | today | with VC cache |
|---|---:|---:|
| Grid count for EXCLUDE filter | 29 s | ~9 s |
| One-time build at finalise | n/a | ~30-60 s on `annotation_workers` |

Build cost amortises across all later queries on the sub-cohort.

## Sequencing

1. `CohortVersion` model + data migration (§1) — landable on its own;
   nothing reads it yet, so no behaviour change.
2. `Cohort.increment_version` / `Cohort.save` hook to `get_or_create`
   the row (§1) + sweeper task (§5).
3. `SubCohortVariantCollection` model + `post_delete` signal +
   `Cohort.get_non_missing_variant_collection` helper (§1).
4. Build task + enqueue helper (§2).
5. Filter-time change (§3) and cache-key tweak.
6. Trigger sites (§4).
7. Tests (§7) — before backfill so harness is proven.
8. Management command for backfill (§6).

§3 + §4 land together; until §4 enqueues, §3 is a no-op (helper always
returns None, regex path stays in use).
