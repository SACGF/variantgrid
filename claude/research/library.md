# library — research notes

Verified against 7c4408c62 on 2026-09-06

`library/` is the model-free plumbing every Django app in VariantGrid stands on: object permissions over
django-guardian, the notification and error-reporting path to Slack, Rollbar and the event log, hover previews and
search summaries, health and uptime checks, per-instance Postgres partitions, the two caching managers, the "major
operation" throttle, the grid filter vocabulary, the test base classes, `library/utils/` (the `library.utils` facade),
and `library/vg/`, the introspection tooling behind `manage.py vg`. It is not an installed app - it has no models,
migrations or URLs of its own, only two abstract mixins that other apps' models inherit - which is why it sits below
`snpdb` in the dependency order and why its upward imports - each one listed in the `ignore_imports` of
`.importlinter` - are worth a paragraph below. `library/AGENTS.md` holds the rules; this is the story behind them. Setting values are in the generated
[settings map](../maps/settings.md), commands in [commands](../maps/commands.md), and `claude/domain.md` has the
vocabulary.

## Flows

### A permission check, from view to guardian

A user-owned model inherits `library/django_utils/guardian_permissions_mixin.py:GuardianPermissionsAutoInitialSaveMixin`
before `models.Model`, so its `save` runs first: on the initial save (no pk yet, or `assign_permissions=True` passed
explicitly) it calls `library/guardian_utils.py:assign_permission_to_user_and_groups`, which grants the owner read and
write on the object and then reads `UserSettings.get_for_user(user).initial_perm_read_and_write_groups` to grant the
same to the lab or organisation groups the user has asked new objects to be shared with. That `UserSettings` import
sits inside the function - `library` cannot import `snpdb` at module level without a cycle. The other function-level
imports of an app are `UserSettingsManager` in `library/utils/export_utils.py:ExportSettings.get_for_request` and
`PartitionArchive` in `RelatedModelsPartitionModel._warn_if_no_archive`; `preview_request` (`genes.enums`,
`variantgrid.perm_path`) imports its at module level. Permission strings are never spelled out by hand:
`library/guardian_utils.py:DjangoPermission.perm` builds `view_<model>` / `change_<model>` from a model, a queryset or
a list, so a rename of the model cannot leave a stale string behind.

A view then asks the object, not guardian. `library/django_utils/guardian_permissions_mixin.py:GuardianPermissionsMixin.can_view`
and `can_write` accept a User or a Group (groups go through `get_group_perms`, users through `user.has_perm` with the
object, so guardian's object backend is consulted), and `get_for_user` is `get_object_or_404` plus `check_permission`,
which raises `PermissionDenied` - the 403 a page test expects. A model that has no permissions of its own delegates:
`get_permission_object` (an instance, e.g. a Sample answering with its VCF) and `get_permission_class` (a class, for
queryset filtering) both default to `self`/`cls`, and every check first follows that pointer. The standard groups -
`library/guardian_utils.py:all_users_group`, `public_group`, `bot_group` - are `get_or_create`d on demand from
`settings.LOGGED_IN_USERS_GROUP_NAME` / `PUBLIC_GROUP_NAME`, and `library/guardian_utils.py:admin_bot` is the system
user that owns automated writes.

The batch forms are where the design shows. `GuardianPermissionsMixin.filter_for_user` resolves the permitted pks
off the bare model in `_permitted_for_user_qs` (falling back to the public group when there is no authenticated user)
and applies them as `pk__in`, skipping even that when the permitted queryset has no filters (a superuser, or a global
model permission). It does not hand the caller's queryset to `get_objects_for_user`: guardian embeds whatever it is
given as a subquery in both its user-perm and its group-perm lookups, so a grid's annotated queryset would have its
joins planned three times over. `filter_writable_for_user` (#1794) is the same idea for the delete column
of a grid page: one query instead of two guardian lookups per row, superusers short-circuited, and
`accept_global_perms=False` because that is how `can_write` asks. A class that overrides one must override the other.
Given the page's `pks`, it filters `UserObjectPermission` / `GroupObjectPermission` by `object_pk__in` the pks as text
(`_object_permission_for_pks_qs`) rather than calling `get_objects_for_user`, whose `object_pk = id::varchar` join
reads every permission row the group holds - flat ~400ms a page at 200k VariantTags on vg-test2, 10-19ms scoped (#1432).
`allow_group_permission_delete` is the third knob: the generic group-permissions delete view will hard-delete anything
whose class returns True, so the base returns False and only the auto-initial-save mixin (user-created objects) opts in,
keeping audit rows like `ClassificationModification` out of its reach.

### Reporting a problem

`library/log_utils.py:report_exc_info` sends `sys.exc_info()` to Rollbar
with the exception message as extra data and prints the traceback; `library/log_utils.py:report_message` is the
non-fatal variant. Both find the current request through `library/request_context.py` (`RequestContextMiddleware` in
`MIDDLEWARE`), so a caller deep in a model never has to thread a request through. Uncaught view exceptions take a
different road: `library/django_utils/rollbar_middleware.py:CustomRollbarNotifierMiddleware` attaches the message and
skips anything derived from `RollbarIgnoreException` (analysis node errors that are the user's problem, not ours).
The nightly health check is an `eventlog` notification: `variantopedia/tasks/server_status_tasks.py:notify_server_status_now`
builds a `NotificationBuilder` and `variantopedia/server_status.py:populate_health_check` fills it from every
`health_check_signal` / `health_check_overall_stats_signal` receiver, sorting the returned
`library/health_check.py:HealthCheckStat` objects into a "last 24 hours" and an "Overall" section. The receivers are
called with `send_robust`, so one app's broken check becomes a line in the message rather than a lost message.
`library/uptime_check.py:retrieve_uptime_response` is the public-safe cousin (a database probe plus
`uptime_check_signal`), rendered on the status page for an external monitor.

### A preview, from a hover to JSON

A model becomes previewable by mixing in `library/preview_request.py:PreviewModelMixin`: `preview_category` defaults to
the verbose name, `preview_icon` to a FontAwesome class, and `preview_enabled` consults
`variantgrid/perm_path.py:get_visible_url_names` through `preview_if_url_visible` so a model whose pages this
deployment does not mount (patients on Shariant) is invisible to search as well. `PreviewModelMixin.preview_with` is the
builder: it folds a `summary` string into `summary_extra` as a dedicated-row `PreviewKeyValue` and hands everything to
`library/preview_request.py:PreviewData.for_object`, which fills the rest from the object - a string pk becomes the
identifier, `get_absolute_url` the internal link (swallowing `NoReverseMatch` for unmounted pages), and
`genome_build(s)` / `annotation_consortium(-ia)` are lifted so search can filter results by build. Other apps add rows
without touching the model: `PreviewData.summary_all` sends `preview_extra_signal` with the object as sender and
appends every `PreviewKeyValue` it gets back, which is how classification counts appear on an Allele preview and
zygosity counts on a Variant's. The hover endpoint is `library/preview_request.py:preview_view` (mounted by snpdb as
`preview_data`), which wraps `db:idx` in a `library/preview_request.py:PreviewRequest` and sends
`preview_request_signal`; more than one answer is reported, the first is used. Search (`snpdb/search.py`) builds its
result rows from the same `PreviewData`, so the icon, category and summary a search result shows are exactly what the
hover shows. `PreviewData.__hash__` exists because search deduplicates results across genome builds (#1115) - a new
field must be added to the hash or two results that differ only there collapse into one.

### A partitioned table, created and dropped

`library/django_utils/django_partition.py:RelatedModelsPartitionModel` is inherited by the model that *holds* a
collection - `CohortGenotypeCollection`, `VariantCollection`, `VariantZygosityCountCollection`, `GeneCoverageCollection`
- and names, in `RECORDS_BASE_TABLE_NAMES` and `RECORDS_FK_FIELD_TO_THIS_MODEL`, the record tables whose rows point at
it. `RelatedModelsPartitionModel.save` detects the first save and calls `create_partition`, whose SQL is
`CREATE TABLE <base>_<label>_<pk> (LIKE <base> INCLUDING INDEXES, CHECK (fk = pk)) INHERITS (<base>)` followed by
`ALTER COLUMN id SET DEFAULT nextval(...)`, because Postgres does not inherit an identity column's default; the
sequence name is asked from `pg_get_serial_sequence` rather than assumed to be `<table>_id_seq`. Reads and writes reach the child through `sql_partition_transformer`, a string replace of the
double-quoted base table name that `library/django_utils/django_queryset_sql_transformer.py:get_queryset_with_transformer_hook`
applies to compiled SQL (the annotation app's `SubVersionPartition` reuses the same hook), or through
`library/django_utils/django_partition.py:temporary_db_table`, which swaps `_meta.db_table` for one query and clears
every field's `cached_col` on the way in and out, or later queries in the process render a mix of both table
names. `delete_related_objects` and `truncate_related_objects` run `DROP` / `TRUNCATE`
inside their own savepoint so a missing partition cannot poison an enclosing `atomic()`, and before a drop
`_warn_if_no_archive` logs when no COMPLETE `snpdb/models/models_partition_archive.py:PartitionArchive` covers the
table. That archive is the #1537 pipeline - `snpdb/partition_archive.py:archive_partitioned_model` schedules a
pg_dump of the child tables, then the drop, then stamps the holder through
`library/django_utils/data_archive_mixin.py:DataArchiveMixin` (#1536: `data_archived_date`, `_by`, `_reason`,
`data_restorable_from`) so the row survives with a pointer to its dump; the warning is telemetry, since VCF re-import
and coverage restore legitimately drop without archiving because the source file is canonical for them.

### Caching managers, in production and in tests

`library/django_utils/django_object_managers.py:ObjectManagerCachingImmutable` and `ObjectManagerCachingRequest` swap
the manager's queryset class for one whose `get()` consults a cache first: the immutable one keeps results in a
module-level dict for the life of the process (GenomeBuild, FlagType), the request one stores them as a
`library/request_context.py` request variable that dies with the request (Allele, Lab, Organization, GeneSymbol, ResolvedVariantInfo). Keys are
`(model, args, frozendict(kwargs))` and only a `get()` on the bare manager queryset uses them
(`QuerySetCachingBase._is_bare_queryset`: no filters, annotations, `values()`, locking or `using()`), so an unhashable
argument or a filtered queryset silently bypasses the cache. Both constructors check
`settings.UNIT_TEST` and leave the plain `QuerySet` in place when it is set, because a cached instance
outlives the test transaction that created it - the same reason `admin_bot` skips its `lru_cache` under `UNIT_TEST`. The consequence for query counting is that a test sees reads production never issues, which is what
`library/django_utils/unittest_utils.py:production_query_count` subtracts: savepoints and any `SELECT ... FROM` one of
`PRODUCTION_CACHED_TABLES`. `library/cache.py:timed_cache` is the function-level equivalent - per-process memo keyed on
positional and keyword arguments with an optional TTL and size cap - used for `GenomeBuild.get_name_or_alias`,
`get_visible_url_names`, evidence keys and feature tips; Redis (`django.core.cache`) is the only cache that crosses
processes, and `CACHE_VERSION` in `default_settings.py` is how its keys are invalidated on deploy.

### Throttling a heavy request

`library/django_utils/major_operation.py:major_operation` is a context manager around any request that can hammer
Postgres (an analysis node grid over millions of variants). It claims a slot with `cache.add` + `cache.incr` on a
per-user Redis counter that carries a safety TTL (`MAJOR_OPERATION_SLOT_EXPIRE_SECONDS`) so a crashed request frees its
slot, raises `TooManyMajorOperationsError` above `MAJOR_OPERATION_MAX_CONCURRENT_PER_USER`, and lowers the connection's
`statement_timeout` to `MAJOR_OPERATION_STATEMENT_TIMEOUT_SECONDS` for the duration, restoring the global value
`variantgrid/wsgi.py` set (connections are reused under `CONN_MAX_AGE`). `_release_slot` tolerates the key having
expired mid-operation - the slow operation is exactly the kind this exists for, and an unguarded `decr` in the
`finally` would turn its result into a 500. `MajorOperationViewMixin` wires it into class-based views by
`major_operation_name`, answering 503 because DataTables ajax POSTs carry their parameters in the body and cannot be
redirect-and-retried the way the analysis node grid is.

### The `vg` tooling

`library/vg/__init__.py` owns the implementation behind `snpdb/management/commands/vg.py`, which only parses arguments,
and `scripts/vg` is the thin dispatcher that runs the Django-free subcommands (`tests`, `outline`, `docs`, the static
maps) without booting the project and delegates the rest to `manage.py vg` under the repo's venv. The split is the
design: `library/vg/repo.py`, `import_graph.py`, `test_selection.py`, `outline.py`, `settings_chain.py`, `docs.py` and
the signals / settings / tasks map generators are AST and filesystem only, so `scripts/vg tests --explain` answers in
well under a second and CI runs `vg map` and `vg docs check` without a database. `library/vg/maps/__init__.py`
renders each map to `claude/maps/<name>.md`, which is gitignored: the SessionStart hook and the `agent-maps.yml` job
rebuild them rather than a committed copy being kept in step. They render under the canonical CI settings module,
because `models` and `urls` depend on `INSTALLED_APPS` and settings-gated URL includes. `library/vg/test_selection.py:select_tests` maps changed
files to test labels through the first-party import graph (a changed module selects every test module that
transitively imports it; templates, JS and `urls.py` select the app's `tests.test_urls`; a migration the app's whole
`tests` package). `library/vg/page.py:render_page` logs the `claude_agent` user in through the test client inside an
`atomic()` that is always rolled back, tracks templates the way Django's test runner does, and reuses
`production_query_count` for `--queries`; `library/vg/inspect/__init__.py:inspect` does the same rollback around one
object's graph by domain kind. `library/vg/docs.py:check_docs` is what keeps this file honest: every backticked path
and `path:Symbol` is resolved by AST against the working tree, live plans included, generated maps excluded. `library/vg/deprecations.py:find_deprecations` runs basedpyright's `reportDeprecated` as
`LOOKAHEAD_PYTHON`, so a call typeshed marks deprecated in a Python we haven't deployed shows before the upgrade, and
splits the hits from our own `@deprecated` markers by name; `manage.py check` under `-W always::DeprecationWarning`
adds the `warnings.warn` deprecations Django and packages raise at import, which no stub records.
`library/vg/health.py:render_health` lists the pylint errors from `lint.txt`; `agent-maps.yml` runs
`scripts/linting/run_pylint.sh --errors-only` on every push and turns each into a warning annotation, not a failure.
`library/vg/health.py:ruff_health` runs `ruff check` live for the count per rule, and the same workflow writes
ruff's per-rule statistics to the job summary, annotating only the pyflakes and bugbear families.
`scripts/linting/pylint_model_utils.py` is a pylint plugin that drops django-model-utils' `TYPE_CHECKING` manager
stubs, without which every query through an `InheritanceManager` reads as assignment-from-no-return.

## Why it is shaped this way

The permission API lives on the model rather than in views because the same question - can this user see this
Sample - is asked by pages, grids, autocompletes, celery tasks and the search index, and each would otherwise grow
its own guardian call with its own bugs. Delegation (`get_permission_object`) exists so the permission rows live on
one object per family (the VCF, not each Sample; the Classification, not each modification) and sharing a VCF is one
`assign_perm`, not hundreds. `filter_for_user` taking a `queryset=` rather than always starting from `cls.objects`
is what lets grids pass their annotated queryset through, and the pk-list shape is the compromise that keeps guardian
out of the planner's way.

Notifications go through a builder rather than a Slack call because the same content has three destinations with
three renderings, and because every send is also an `Event` row: the event log is the audit trail and the overflow
buffer for Slack's size limit. Reporting through the request context rather than passing requests is the price of having
model methods that can report problems from a view or a task alike.

Partitions are inheritance children rather than declarative partitions because the design needs no range or hash
routing: each collection is its own table so that deleting a cohort or an
annotation version is a `DROP TABLE`, not millions of row deletes, and the base table stays empty. That is why a
plain ORM query on the base model returns nothing and why every read must go through the transformer or
`temporary_db_table` (`claude/guides/operations.md#scale`).

The caching managers exist because GenomeBuild, Lab and Allele are fetched thousands of times per request from code
that has no natural place to pass them down; process-wide caching is reserved for tables never updated in place and
request caching for the rest, and neither invalidates, because invalidation would only work in the process that
did the write. They are off under `UNIT_TEST` rather than cleared between tests because a cached instance (a User,
say) survives the test's rollback, which no clear-on-teardown reliably catches.

`vg` lives in `library/` rather than `scripts/` so its Django-dependent halves can import models, and its
Django-free halves are kept import-light on purpose: the agent loop calls `scripts/vg tests` and `docs check` on every
edit, and a second of Django boot per call would make that unusable.

## Traps

`assign_permission_to_user_and_groups` requires a saved `UserSettings` chain to answer
`initial_perm_read_and_write_groups`; a freshly created test user without the usual groups gets only personal
permissions, which is fine until a test expects the lab to see the object. `GuardianPermissionsAutoInitialSaveMixin`
must come before `models.Model` in the bases and the model must have a `user` attribute, or `save` raises
`ValueError` after the row is already written. `filter_for_user(user, queryset=...)` with a queryset from a
*different* model than `cls` silently filters the wrong pks - the delegation path expects the queryset to be of `cls`.

`report_exc_info` outside an `except` reports `None`.

`RelatedModelsPartitionModel.get_partition_table()` with no argument raises once a holder has more than one base table
(VariantZygosityCountCollection), and `save()` on an existing holder never recreates a dropped partition - re-import
code calls `create_partition` explicitly. A `ProgrammingError` from a drop is logged and swallowed unless
`LOG_PARTITION_WARNINGS` is False (URLTestCase sets it so), so a typo in `RECORDS_BASE_TABLE_NAMES` shows up as a
warning, not a failure. `sql_partition_transformer` is a string replace on the quoted table name: a query that joins
the base table twice gets both occurrences rewritten, which is usually right and occasionally not.

`ObjectManagerCachingImmutable` on a model whose rows are ever updated serves stale instances until every worker
restarts; `Meta.base_manager_name = 'objects'` is needed for the cache to cover related-object access too. The request
cache needs `RequestContextMiddleware`; in a celery task there is no request and every `get` hits the database, which is
correct but surprising in a query count. `timed_cache`'s `cache_clear` only clears the calling process, and keys
on the positional arguments plus the sorted keyword arguments, so `f(1)` and `f(x=1)` are cached separately.

`major_operation` counts against a single per-user key shared by every operation type; a user with twenty grids open
in twenty tabs is at the limit for the analysis grid too. The shorter `statement_timeout` goes through
`library/django_utils/database_utils.py:pg_settings`, which reads each setting first and puts it back on exit, so
nested major operations restore the outer one's value. Under `MAJOR_OPERATION_LIMITS_ENABLED=False` neither the cap nor
the shorter timeout applies.

`library/vg/inspect/__init__.py:inspect` and `library/vg/page.py:render_page` call `set_rollback(True)` only after the
work is done; marking rollback first makes Django refuse every query in the block. `vg map` under any settings module
but `variantgrid.settings.env.github_actions` produces maps that differ from CI's.
`library/vg/docs.py` skips fenced code blocks and treats a bare filename such as `default_settings.py` as a citation, so a filename mentioned in
prose must exist somewhere in the tree.
