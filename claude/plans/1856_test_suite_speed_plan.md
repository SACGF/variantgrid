# #1856 — Speed up the unit test suite

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-09
Status: draft

[#1856](https://github.com/SACGF/variantgrid/issues/1856) profiled the suite (3,245 tests, 8-core laptop,
Django 6.1, PostgreSQL 18): 128 s wall with `--keepdb --parallel 4`, 470 s summed CPU, 11 min serial.
The cost is flat across modules, so the wins are in the shared fixtures every class pays for, in this order:

| Cost centre | Serial CPU | Phase |
|---|---|---|
| PBKDF2 password hashing (75 hashes, ~1 s each) | 68 s | 1 |
| `get_fake_annotation_version` (240 classes × 0.32 s, ontology re-import inside it) | 77 s | 2 |
| `assign_permission_to_user_and_groups` (1,562 calls × 27 ms) | 41 s | 3 |
| `create_fake_cohort` (133 calls × 0.26 s, mostly the three permission calls above) | 35 s | 3 |

Target after phases 1-3: 4-worker wall from ~128 s to ~85-90 s. Each phase lands on its own, and the
acceptance for each is the same measurement, recorded in the table at the end:

```bash
/usr/bin/time -f "wall %e s, cpu %U+%S s" python3 manage.py test --keepdb --parallel 4
```

## Data

No model changes. Phase 2 changes *when* fixture rows are created (once per run instead of once per
class), not what they are. The one fixture row that changes shape is `GeneAnnotationVersion`:
`annotation/fake_annotation.py:get_fake_annotation_version` currently passes `gnomad_import_date=timezone.now()`
as a lookup key to `get_or_create`, so every call creates a new row; the date moves into `defaults`.

## Phase 1 — MD5 password hasher under `UNIT_TEST`

In the `UNIT_TEST` block of `variantgrid/settings/components/default_settings.py`, alongside the
`LocMemCache` override:

```python
if UNIT_TEST:
    # Django's PBKDF2 is ~1 s per hash (1.5M iterations); the suite creates ~75 users
    PASSWORD_HASHERS = ["django.contrib.auth.hashers.MD5PasswordHasher"]
```

Measured: `User.objects.create_user(password=...)` 1,268 ms → 2 ms. Nothing in the suite asserts on
the hash format and the `client.login()` calls keep working because the same hasher verifies. Add a line to
`claude/guides/testing.md` under the `UNIT_TEST` settings paragraph.

Expected: ~70 s CPU off the run (15%).

## Phase 2 — seed the fake annotation version once per run

### Why per-class is the wrong granularity

Every `TestCase` runs its `setUpTestData` inside a transaction that is rolled back at the end of the class,
so the 240 classes that call `get_fake_annotation_version` each re-import the ontology OWL (`load_hpo`,
0.18 s, pronto parse + chardet) and recreate the annotation version chain. The rows are identical every
time. There is no `TransactionTestCase` in the suite (nothing flushes tables between classes), so rows
created before the first test survive the whole run.

### Where the seed goes

`variantgrid/test_runner.py:VariantGridTestRunner.setup_databases` already wraps Django's database setup.
Django creates the main test database and then, when `parallel > 1`, clones it for each worker inside the
same call, so the seed has to run between the two steps:

1. Call `super().setup_databases()` with `self.parallel` temporarily set to 1 (main database created and
   migrated, `_check_kept_test_db_matches_disk` runs as now).
2. Seed on the `default` connection: `get_fake_annotation_version(build)` for GRCh37 and GRCh38 (the two
   builds the suite uses; `builds_with_annotation()` is not used by any fixture).
3. Clone `self.parallel` times with `connection.creation.clone_test_db(suffix=str(index + 1), verbosity=..., keepdb=self.keepdb)`,
   the loop Django's `setup_databases` would have run. `_drop_test_db_clones` keeps running first, as now.

Under `--keepdb` the seeded rows persist in `test_snpdb` between runs; the seed is idempotent
(below) so a second run is a handful of `get_or_create` lookups. Without `--keepdb` it runs against the
fresh database. The seed stays outside any transaction, so a failure raises and aborts the run with the
real traceback rather than surfacing as 240 class errors.

### Making the fixture idempotent and cheap on re-entry

The per-class calls stay in the tests - a class must still be able to run alone under
`manage.py test <label>`, and against a database seeded by the runner they cost a few queries. For that:

- `annotation/fake_annotation.py:get_fake_annotation_version`: move `gnomad_import_date` from the
  `get_or_create` lookup into `defaults` (see Data). The other `get_or_create` calls already key on stable
  values.
- `ontology/tests/test_data_ontology.py:create_ontology_test_data`: call `load_biomart` / `load_hpo` with
  `force=False` and catch `ontology/ontology_builder.py:OntologyBuilderDataUpToDateException`, which
  `ensure_hash_changed` raises when the file hash matches the previous `OntologyImport`. That turns the
  0.18 s parse into one md5 of a small file and one query.

`create_test_ontology_version` is already a lookup-first function.

### What can break, and how to find it

- A test that asserts on a count or `first()` of `AnnotationVersion`, `VariantAnnotationVersion`,
  `GeneAnnotationVersion`, `OntologyVersion`, `OntologyImport` or `OntologyTerm` now sees the seeded rows.
  `grep` finds no such count assertion today; the run tells.
- A test that creates its own `VariantAnnotationVersion` with `status=ACTIVE` for a build now has two, and
  `AnnotationVersion.latest(build)` may return the seeded one. Fix the test to use the fixture (or a
  non-active status), not the seed.
- A test that mutates the shared rows in its transaction is rolled back at class end - fine.

Run the whole suite serially once (`--keepdb`, no `--parallel`) as well as parallel: serial catches an
order dependency a worker split can hide.

Expected: ~80 s CPU off the run. If the seed turns out to break more tests than is worth fixing, the fallback
is smaller and self-contained: cache the parsed `pronto.Ontology` in-process keyed by file hash inside
`load_hpo` (~27 s), and keep per-class creation.

## Phase 3 — batch `assign_permission_to_user_and_groups`

`library/guardian_utils.py:assign_permission_to_user_and_groups` is production code (16 call sites: VCF
import per sample, cohort, analysis, classification…), so the saving lands in real imports too. Today one
call is ~27 ms:

- `assign_perm(...)` × (2 + read groups + write groups); each guardian call is a `ContentType` lookup, a
  `Permission.objects.get`, and a `get_or_create`.
- `UserSettings.get_for_user(user)` → `get_settings_overrides`: `get_or_create` on user, lab and
  organisation overrides plus `GlobalSettings`, then `initial_perm_read_and_write_groups` runs one
  `SettingsInitialGroupPermission` query per override. About 16 ms and 5+ queries.

New shape, same behaviour and same rows written:

1. Resolve the content type once and the two `Permission` rows in one query
   (`Permission.objects.filter(content_type=ctype, codename__in=[read, write])`). guardian's
   `_ensure_permission` accepts a `Permission` instance and skips its lookup.
2. Write the user rows with `UserObjectPermission.objects.bulk_create([...], ignore_conflicts=True)` and the
   group rows with `GroupObjectPermission.objects.assign_perm_to_many(perm, groups, obj, ignore_conflicts=True)`,
   one `bulk_create` per permission. `ignore_conflicts` replaces the per-row `get_or_create`; the unique
   constraint on (user/group, permission, content_type, object_pk) is what guardian relies on anyway.
3. `UserSettings.get_initial_perm_read_and_write_groups`: one query over
   `SettingsInitialGroupPermission.objects.filter(settings__in=overrides)`, applied in override order in
   Python, instead of a query per override.

Target: under 10 ms and ~6 queries per call. The existing inline `from snpdb.models import UserSettings`
in that function is a pre-existing cycle (`snpdb.models` imports `library.guardian_utils`) and stays as it
is in this change.

Tests: `library/tests/` gets one test that the batched function writes exactly the user and group rows the
old loop did for a user with one read-only and one read-write initial group (that is the rule that is easy
to get wrong; guardian's own row writing is theirs to test). Then `scripts/vg tests --explain` for the
callers.

Expected: ~25 s CPU off the run; `create_fake_cohort` drops with it.

## Phase 4 — docs

`claude/guides/testing.md`:

- Serial run is ~11 min, not ~25; the suite is 3,245 tests (the file says 2,741).
- `--parallel 4` on an 8-core box; 8 workers doubled the CPU (990 s) and lost wall time (142 s) from Postgres
  contention.
- The runner seeds the fake annotation versions for GRCh37/GRCh38 before cloning, and what that means for a
  test that counts annotation or ontology rows (Phase 2).
- The `UNIT_TEST` settings paragraph gains the password hasher (Phase 1).

`scripts/vg docs check` after the edit.

## Out of scope, noted for later

- **`TestData` deep-copies** (23 s): Django copies every `setUpTestData` attribute on access, and a
  `GenomeBuild` / `VariantAnnotationVersion` carrying populated `cached_property` values costs 13-22 ms per
  access. Worth its own measurement once phases 1-3 have landed and the profile is re-run.
- **Audit tests** (~20 s, `library/tests/test_decorator_audit.py` AST-parses every source file,
  `library/tests/test_signal_receiver_registration.py` imports every app): each blocks one worker for its
  duration; caching the scan or moving it to lint is a small separate change.
- **Per-test timing runner**: the JSONL timing runner used for the profile could become a
  `VG_TEST_TIMING=<file>` option on `VariantGridTestRunner` if the suite needs re-profiling regularly.
- `snpdb.tests.test_query_counts.ViewSampleScalingTest.test_view_sample_query_count_flat_with_more_trios`
  failed once serially only (60 vs 59 queries); order-sensitive, unrelated to this plan but the serial run
  in Phase 2 will show whether it recurs.

## Measurements

| After | Wall (`--parallel 4`) | Summed CPU | Notes |
|---|---|---|---|
| baseline (issue) | 128 s | 470 s | |
| Phase 1 | | | |
| Phase 2 | | | |
| Phase 3 | | | |
