# library — agent notes
Owns: model-free plumbing — Guardian permissions, notifications/error reporting, previews, health checks, caching
managers, partitioning, grid filter rules/CSV export, test base classes, `library/utils/` (re-exported as `library.utils`).
Start with: `guardian_utils.py` (DjangoPermission, group getters, assign_permission_to_user_and_groups) ·
`django_utils/guardian_permissions_mixin.py` (GuardianPermissionsMixin) · `log_utils.py` (NotificationBuilder,
report_exc_info) · `preview_request.py` (PreviewModelMixin, PreviewData, the preview signals) ·
`django_utils/unittest_utils.py` (URLTestCase, query profiling) · `django_utils/django_object_managers.py`
Patterns here:
- Give a new user-owned model permissions by inheriting
  `library/django_utils/guardian_permissions_mixin.py:GuardianPermissionsAutoInitialSaveMixin` before `models.Model`;
  it needs a `user` field and calls `library/guardian_utils.py:assign_permission_to_user_and_groups` on first save,
  which reads the user's initial read/write groups from `UserSettings`.
- Check access through `library/django_utils/guardian_permissions_mixin.py:GuardianPermissionsMixin`
  (`can_view`, `check_can_write`, `filter_for_user`, `filter_writable_for_user`, `get_for_user`) rather than guardian
  shortcuts; a model may delegate to another object's permissions via `get_permission_class` /
  `get_permission_object`, and `allow_group_permission_delete` is opt-in (default False).
- Build permission strings with `library/guardian_utils.py:DjangoPermission.perm(obj, DjangoPermission.READ|WRITE)`;
  standard groups are `library/guardian_utils.py:all_users_group` / `public_group` / `bot_group`, the system user `admin_bot`.
- Notify humans with `library/log_utils.py:AdminNotificationBuilder` (Slack, plus superuser email when
  `is_communication=True`) or snpdb's `LabNotificationBuilder`; add blocks with `add_header`/`add_field`/`add_markdown`
  (Markdown, converted to HTML and Slack) and call `send()`. `library/log_utils.py:send_notification` is the raw
  Slack hook, falls back to Rollbar when Slack is unconfigured, and truncates blocks at `SLACK_CHARACTER_LIMIT`.
- Report caught exceptions with `library/log_utils.py:report_exc_info` (Rollbar + traceback) and non-fatal problems
  with `library/log_utils.py:report_message`; both find the current request via django-threadlocals.
- Make a model hover-previewable and searchable by implementing `library/preview_request.py:PreviewModelMixin`
  (`preview_category`, `preview_icon`, `preview`) and returning `library/preview_request.py:PreviewData` via
  `preview_with`; contribute extra key/values from another app with `@receiver(preview_extra_signal, sender=Model)`
  returning `library/preview_request.py:PreviewKeyValue` (use `PreviewKeyValue.count` for counts).
- Register an app health check with `@receiver(signal=health_check_signal)` returning `library/health_check.py:HealthCheckRecentActivity`
  / `HealthCheckTotalAmount` / `HealthCheckAge` stats; uptime probes go through `library/uptime_check.py:uptime_check_signal`.
- Wrap DB-heavy request handlers in `library/django_utils/major_operation.py:major_operation` (per-user concurrency cap
  plus a lower `statement_timeout`), handling `TooManyMajorOperationsError`; grids use `MajorOperationViewMixin`.
- Grid column filters and FilterNode share one vocabulary: `library/django_utils/filter_rules.py:FILTER_OPERATIONS`
  and `rules_to_q`; streaming CSV downloads go through `library/django_utils/grid_export.py:grid_export_csv`.
- Before writing a helper, check `library/utils/`: `collection_utils.py` (`batch_iterator`, `group_by_key`, `first`,
  `get_single_element`, `sorted_nicely`, `invert_dict`, `LazyAttribute`) · `database_utils.py` (`queryset_to_sql`,
  `dictfetchall`, `sql_delete_qs`) · `text_utils.py` (`pretty_label`, `limit_str`, `format_percent`) ·
  `html_utils.py` (`html_id_safe`, `sanitize_html`, `html_to_text`) · `json_utils.py` (`force_json`, `strip_json`,
  `JsonDiffs`) · `hash_utils.py` (`md5sum_str`, `sha256sum_str`, `stable_dict_hash`) · `file_utils.py`
  (`open_handle_gzip`, `mk_path_for_file`, `name_from_filename`) · `date_utils.py` (`calculate_age`, `parse_yymm`) ·
  `diff_utils.py` (`diff_text`, `MultiDiff`) · `export_utils.py` (`ExportRow`, `export_column`) · `model_utils.py`
  (`ArrayLength`, `model_has_field`) · `django_utils.py` (`is_ajax`, `render_ajax_view`, `refresh_for_update`) ·
  `class_utils.py` (`import_class`, `get_all_subclasses`) · `os_utils.py` (`execute_cmd`) · `timer_utils.py`
  (`get_timer`) · `misc_utils.py` (`empty_to_none`, `ChoicesEnum`, `iter_http_lines`) · `xml_utils.py` (`XmlParser`);
  plus `library/django_utils/__init__.py:require_superuser` and `library/django_utils/django_postgres.py:copy_from_file`.
Gotchas:
- `library/django_utils/django_object_managers.py:ObjectManagerCachingImmutable` caches `.get()` results
  process-wide forever and `ObjectManagerCachingRequest` per request; both are plain managers under
  `settings.UNIT_TEST`, so tests see queries production never runs (GenomeBuild, GeneSymbol, FlagType, Allele, Lab,
  Organization, ResolvedVariantInfo). Put a caching manager only on tables never updated in place.
- `library/guardian_utils.py:assign_permission_to_user_and_groups` imports `snpdb.models.UserSettings` inside the
  function: `library` depends on `snpdb` at runtime, and that inline import is the one sanctioned cycle-break here.
- `library/guardian_utils.py:admin_bot` skips its lru_cache under `UNIT_TEST` because a cached User outlives the
  test transaction rollback; do the same for any module-level cache of a model instance.
- `library/django_utils/guardian_permissions_mixin.py:GuardianPermissionsMixin.filter_for_user` resolves permitted
  pks off the bare model and applies them as `pk__in`, so pass the caller's annotated queryset as `queryset=` rather
  than the class; Guardian embeds whatever it is given in both its lookups. Override `filter_writable_for_user`
  alongside `can_write` or the two drift.
- `library/django_utils/django_partition.py:temporary_db_table` must bracket any query against a partition table;
  swapping `_meta.db_table` by hand leaves `Field.cached_col` pointing at the partition for the life of the process.
- `library/cache.py:timed_cache` keys on `(*args, kwargs.items())` and holds results in module memory, so each celery
  and web worker process has its own copy. `clear_cached_property` there and
  `library/utils/misc_utils.py:invalidate_cached_property` do the same job — use either, add no third.
- Signal receivers under `<app>/signals/` only register when imported from the app's `ready()`;
  `library/tests/test_signal_receiver_registration.py` fails the suite if such an import is dropped, and
  `library/tests/test_decorator_audit.py` fails it on `@classmethod` stacked with `@property` (removed in 3.13).
- A plain script needs `DJANGO_SETTINGS_MODULE` set explicitly (e.g. `variantgrid.settings.env.vgtest2`);
  `variantgrid/settings/__init__.py` otherwise derives the module from the hostname (`vg-test2` → `vgtest2`,
  `env_developers/` checked before `env/`) and only logs an error when no file matches.
Tests: `TEST_RUNNER` is `variantgrid/test_runner.py:VariantGridTestRunner`, which swaps `ClinGenAlleleRegistryAPI`
and `TranscriptSequenceFetcher` for recorded mocks that raise, naming the fixture to add, when asked for anything
unrecorded; `FastaRecordingRunner` regenerates the sparse test fastas. Page tests extend
`library/django_utils/unittest_utils.py:URLTestCase` (Celery eager, plain static storage, annotation web resources
off) and use `_test_urls` / `_test_datatable_urls` / `_test_autocomplete_urls` / `_test_datatables_grid_urls_contains_objs`.
Count queries with `library/django_utils/unittest_utils.py:production_query_count`, which drops savepoints and the
production-cached tables above. `VG_QUERY_PROFILE=<file>` makes every `URLTestCase` GET append a JSON line of query
count, duplicates and timings via `library/django_utils/unittest_utils.py:QueryProfilingClient`; `VG_QUERY_TRACE=<regex>`
also writes Python stacks for matching SQL to `<file>.trace`. Under `UNIT_TEST` the cache is in-process LocMem,
Rollbar and axes are off, and Postgres JIT is disabled. Library's own tests: `library/tests/`, `library/django_utils/tests/`.
Deep reference: __library_readme.md · claude/research/library.md
