# #1432 Grid write-permission batching: benchmark at scale

Written by Claude Opus 5.5 (claude-opus-5-5), 2026-09-24
Status: in progress - `pks` implemented (Outcome); awaiting the vg-test2 re-run under **Verification**

## Goal

Measure whether resolving a grid page's write permissions in one query beats checking each row, on a database
with realistic permission volumes. Run it on vg-test2 (read-only). The result decides whether #1432 keeps the batched
check on the VariantTag detail grid and whether the batched form needs scoping to the page.

## Background

Grids resolve "can this user write this row?" for a whole page through
`snpdb/views/datatable_view.py:DatatableConfig._writable_pks_for_page`, which runs
`Model.filter_writable_for_user(user).filter(pk__in=page_pks)`. That pattern has been on master since 59338a250
(2026-08-30), behind `render_delete` on these grids:

| Grid | Model |
|---|---|
| `variantopedia/grids.py:VariantTagsColumns` | VariantTag |
| `snpdb/grids.py:SamplesListColumns` | Sample |
| `snpdb/grids.py:VCFListColumns` | VCF |
| `snpdb/grids.py:CohortListColumns` | Cohort |
| `patients/grids.py:PatientListColumns` | Patient |
| `analysis/grids.py:AnalysesListColumns` | Analysis |

#1432 also moves `variantopedia/grids.py:VariantTagDetailColumns` (one variant + one tag, usually a handful of
rows) from a per-row `VariantTag.objects.get(pk).can_write(user)` onto the same batched check.

The concern: `analysis/models/models_variant_tag.py:VariantTag.filter_writable_for_user` (like the Sample, Cohort
and Patient overrides) composes `Q(pk__in=own)`, where `own` comes from
`library/django_utils/guardian_permissions_mixin.py:GuardianPermissionsMixin.filter_writable_for_user`, a
Guardian lookup over the whole table. On a dev database with 39 tags, `EXPLAIN ANALYZE` shows the page's
`pk__in` staying on the outer query while that subquery seq-scans all of `analysis_varianttag` and hashes every
VariantTag permission row for the user and their groups. SA Path has ~412k VariantTags, and a lab group can hold a
permission row on most of them, so the batched query's cost may grow with the table rather than the page.

## Method

Everything below is SELECTs through the ORM and uses code already on master, so it runs on vg-test2 as deployed.

1. `python3 manage.py vg status` and `git log -1 --format=%h`. Confirm
   `analysis/models/models_variant_tag.py:VariantTag.filter_writable_for_user` exists in the deployed checkout.
2. Save the script below to the session scratchpad and run it with `python3 manage.py shell < bench.py`.
3. Report back in the format under **Report**.

For each model it picks the *heavy user*: the active non-superuser in the group holding the most object permission
rows for that model (superusers short-circuit to `objects.all()` and measure nothing). It builds two page shapes:

- `page100`: the 100 newest rows the heavy user can view (`filter_for_user`), which is the list grid's page.
- `detail` (VariantTag only): every tag on the (variant, tag) pair with the most taggings, which is the detail grid's page.

It then times three ways of answering "which of these can the user write?":

- `per_row`: `Model.objects.get(pk=pk).can_write(user)` for each pk (the pre-batching behaviour).
- `batched`: `Model.filter_writable_for_user(user).filter(pk__in=pks)` (master's `_writable_pks_for_page`).
- `scoped` (VariantTag only): the same rule with Guardian starting from the page's rows, as a candidate fix.

Each is warmed once, then timed as the best of 5, and all three must return the same set.

```python
import time

from django.contrib.auth.models import User
from django.contrib.contenttypes.models import ContentType
from django.db import connection
from django.db.models import Count, Q
from django.test.utils import CaptureQueriesContext
from guardian.models import GroupObjectPermission, UserObjectPermission
from guardian.shortcuts import get_objects_for_user

from analysis.models import Analysis, VariantTag
from patients.models import Patient
from snpdb.models import VCF, Cohort, Sample

MODELS = [VariantTag, Sample, VCF, Cohort, Patient, Analysis]
PAGE = 100
REPS = 5


def timed(fn):
    fn()
    best = None
    for _ in range(REPS):
        with CaptureQueriesContext(connection) as ctx:
            start = time.perf_counter()
            result = fn()
            elapsed = (time.perf_counter() - start) * 1000
        best = elapsed if best is None else min(best, elapsed)
    return result, best, len(ctx.captured_queries)


def heavy_user(model):
    ct = ContentType.objects.get_for_model(model)
    top = (GroupObjectPermission.objects.filter(content_type=ct).values("group")
           .annotate(n=Count("id")).order_by("-n").first())
    if top is None:
        return None, 0
    user = User.objects.filter(groups=top["group"], is_superuser=False, is_active=True).first()
    return user, top["n"]


def per_row(model, user, pks):
    return {pk for pk in pks if model.objects.get(pk=pk).can_write(user)}


def batched(model, user, pks):
    return set(model.filter_writable_for_user(user).filter(pk__in=pks).values_list("pk", flat=True))


def scoped_variant_tag(user, pks):
    page_qs = VariantTag.objects.filter(pk__in=pks)
    own = get_objects_for_user(user, VariantTag.get_write_perm(), klass=page_qs, accept_global_perms=False)
    qs = page_qs.filter(Q(analysis__in=Analysis.filter_writable_for_user(user)) |
                        Q(analysis__isnull=True, pk__in=own))
    return set(qs.values_list("pk", flat=True))


for model in MODELS:
    ct = ContentType.objects.get_for_model(model)
    print(f"## {model.__name__}: rows={model.objects.count()} "
          f"group_perms={GroupObjectPermission.objects.filter(content_type=ct).count()} "
          f"user_perms={UserObjectPermission.objects.filter(content_type=ct).count()}")
    user, group_perm_rows = heavy_user(model)
    if user is None:
        print("  no group permissions - skipped")
        continue
    print(f"  heavy user id={user.pk} (group holds {group_perm_rows} perm rows)")

    shapes = {"page100": list(model.filter_for_user(user).order_by("-pk").values_list("pk", flat=True)[:PAGE])}
    if model is VariantTag:
        pair = (VariantTag.objects.values("variant", "tag").annotate(n=Count("id")).order_by("-n").first())
        shapes["detail"] = list(VariantTag.objects.filter(variant=pair["variant"], tag=pair["tag"])
                                .values_list("pk", flat=True))

    for shape, pks in shapes.items():
        runs = {"per_row": lambda: per_row(model, user, pks),
                "batched": lambda: batched(model, user, pks)}
        if model is VariantTag:
            runs["scoped"] = lambda: scoped_variant_tag(user, pks)
        results = {}
        for name, fn in runs.items():
            result, ms, queries = timed(fn)
            results[name] = result
            print(f"  {shape:8s} rows={len(pks):4d} {name:8s} {ms:9.1f} ms {queries:5d} queries "
                  f"writable={len(result)}")
        print(f"  {shape:8s} same answer: {len({frozenset(r) for r in results.values()}) == 1}")

    if model is VariantTag:
        pks = shapes["page100"]
        qs = VariantTag.filter_writable_for_user(user).filter(pk__in=pks).values_list("pk", flat=True)
        print(qs.explain(analyze=True, buffers=True))
```

## Report

Paste back:

1. Host, deployed sha, and the `## Model: rows=… group_perms=… user_perms=…` line for each model.
2. The timing lines for each model and shape, as printed.
3. Every `same answer:` line. A `False` is a correctness bug, so give the model, shape and the three `writable=` counts.
4. The VariantTag `EXPLAIN (ANALYZE, BUFFERS)` output in full.
5. For each model and shape: which of `per_row` / `batched` / `scoped` was fastest, and by how much.

## Decision

- `batched` at or below `per_row` for every model and shape: the pattern stands, and #1432's VariantTag detail
  grid change ships as is.
- `batched` above `per_row` for VariantTag, with `scoped` fast: `filter_writable_for_user` gains an optional
  `queryset` argument, mirroring `filter_for_user`, so the Guardian lookup starts from the page's rows. The mixin
  and every override (VariantTag, Sample, Cohort, Patient, Analysis) take it, and
  `snpdb/views/datatable_view.py:DatatableConfig._writable_pks_for_page` passes the page's queryset.
- Another model showing the same growth as VariantTag goes into the same change.
- If vg-test2's VariantTag count is well under SA Path's ~412k, the VariantTag verdict needs the same script run on
  an SA Path host, which the user arranges.

## Results (vg-test2, aee138551, 2026-09-24)

Data: VariantTag 200,019 rows (200,015 synthetic, spread 2021-2026), one group holding a `view` and `change` row
on each (400,000 group permission rows); heavy user 2, non-superuser. Sample has no group permissions and was
skipped. VCF, Cohort, Patient and Analysis tables are tiny (17-51 rows). Every `same answer:` line was `True`.

`analysis_varianttag` and `guardian_groupobjectpermission` had **never been analyzed** (pg_stat said 15 and 6,116
rows). The first run was on those stats; the tables were then `ANALYZE`d and the script re-run. Timings in ms,
best of 5; query counts past the first VariantTag shape are unreliable (Django's 9,000-query log cap).

| Model / shape | rows | per_row | batched | scoped (stale stats) | scoped (analyzed) | literal |
|---|---|---|---|---|---|---|
| VariantTag page100 | 100 | 335 | 398 | 1,463 | 11.5 | 10.4 |
| VariantTag detail | 555 | 1,857 | 408 | 7,919 | 7,851 | 18.7 |
| VCF page100 | 14 | 51.8 | 3.8 | - | - | - |
| Cohort page100 | 16 | 85.2 | 8.0 | - | - | - |
| Patient page100 | 17 | 55.9 | 4.1 | - | - | - |
| Analysis page100 | 18 | 74.2 | 5.0 | - | - | - |

(Analysis's heavy user had writable=0, so it only exercised the "no" path.)

Findings:

- `batched` is flat at ~400 ms for VariantTag whatever the page size - the cost follows the table, as feared.
  `EXPLAIN (ANALYZE, BUFFERS)`: 544-566 ms, almost all in the `analysis__isnull=True, pk__in=own` branch, which
  seq-scans all of `analysis_varianttag`, hash-aggregates 200k `id::varchar` values and hashes all 200k of the
  group's `change_varianttag` rows. ANALYZE did not change it. At SA Path's ~412k tags expect roughly double.
- `scoped` (Guardian `get_objects_for_user(klass=page_qs)`) is **not** a fix. Guardian joins
  `object_pk::text = id::varchar::text`, which the planner can neither estimate (rows=1) nor drive the
  `(group_id, permission_id, object_pk)` index from, so it still reads every permission row the group holds. On
  page100 it happened to get a good plan after ANALYZE; on the detail shape (555 pks) it chose a nested-loop semi
  join over 200k x 555 (110M join-filter rows, 16 s under EXPLAIN).
- `literal` - Guardian's rows filtered by the page's pks as a **text list**, which the unique index answers
  directly - is fast and stable for both shapes and returns the same set:

```python
def literal_variant_tag(user, pks):
    perm = Permission.objects.get(content_type=ContentType.objects.get_for_model(VariantTag),
                                  codename="change_varianttag")
    str_pks = [str(pk) for pk in pks]
    own = set(UserObjectPermission.objects.filter(user=user, permission=perm, object_pk__in=str_pks)
              .values_list("object_pk", flat=True))
    own |= set(GroupObjectPermission.objects.filter(group__user=user, permission=perm, object_pk__in=str_pks)
               .values_list("object_pk", flat=True))
    qs = VariantTag.objects.filter(pk__in=pks).filter(Q(analysis__in=Analysis.filter_writable_for_user(user)) |
                                                      Q(analysis__isnull=True, pk__in=[int(p) for p in own]))
    return set(qs.values_list("pk", flat=True))
```

  (4 queries as written; the permission lookup and the two Guardian reads can fold into subqueries of one.)

## Outcome

Neither Decision branch applies as written: `batched` loses to `per_row` on VariantTag page100 (398 vs 335 ms), and
`scoped` is not reliably fast. The detail grid change is still a clear win (408 vs 1,857 ms), and the small models
are 11-15x faster batched.

Recommended change for #1432, in place of the `queryset` argument:

- `filter_writable_for_user(user, pks=None)`: with `pks`, the Guardian part filters
  `UserObjectPermission` / `GroupObjectPermission` by `object_pk__in=[str(pk) for pk in pks]` (plus
  `content_type`/`permission`, and user or user's groups) instead of calling `get_objects_for_user` over the whole
  table, and the outer queryset is limited to `pk__in=pks`. `pks=None` keeps today's behaviour.
- The mixin and every override (VariantTag, Sample, Cohort, Patient, Analysis) take `pks`; overrides that
  delegate to a permission object (VariantTag -> Analysis) keep delegating unscoped, since Analysis is small.
- `snpdb/views/datatable_view.py:DatatableConfig._writable_pks_for_page` passes its `pks`.
- Implemented as described: `library/django_utils/guardian_permissions_mixin.py:GuardianPermissionsMixin._object_permission_for_pks_qs`
  holds the scoped Guardian lookup, and Analysis and AnalysisTemplate take `pks` too.
- An SA Path run is still worthwhile for scale, but no longer decides the design: the literal form's cost follows
  the page, not the table.

## Verification

On vg-test2, with the same synthetic data still in place:

1. Upgrade to master (`scripts/upgrade.sh --quick`, then the usual restart and `vg status`) and confirm the deployed
   checkout has `filter_writable_for_user(cls, user, pks=None)` in `library/django_utils/guardian_permissions_mixin.py`.
2. Re-run the Method script with one change: in `runs`, add
   `"pks": lambda: set(model.filter_writable_for_user(user, pks=pks).values_list("pk", flat=True))`, and point the
   closing `EXPLAIN` at `VariantTag.filter_writable_for_user(user, pks=pks)`.
3. Report as before, plus the `pks` line for every model and shape.

Pass: `pks` lands at 10-20 ms for VariantTag page100 and detail, at or below `batched` for every other model, and
every `same answer:` is `True`; the EXPLAIN reads `guardian_groupobjectpermission` through an index on
`object_pk` rather than a scan of the group's rows. On a pass the plan's Status becomes landed and the plan is
deleted.
