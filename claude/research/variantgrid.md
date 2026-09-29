# variantgrid — research notes

Verified against a96540a68 on 2026-09-25

`variantgrid/` is the Django project package rather than an app: it has no models, and owns the settings, the root
URLconf, the Celery app and beat schedule, the WSGI entry point, the test runner, the static-files and template trees
for every site, and a handful of project-wide views (login redirect, error pages, version, the capabilities endpoint).
Its tasks and URLs are in the generated maps ([tasks](../maps/tasks.md#variantgrid), [urls](../maps/urls.md#variantgrid),
[settings](../maps/settings.md), [commands](../maps/commands.md)); which host loads which settings file, the queues and
the authentication surface are in `claude/guides/operations.md`, which this doc does not repeat.

## Flows

### Settings: hostname to one flat module

`variantgrid/settings/__init__.py` runs whenever anything imports `variantgrid.settings`. With the default
`DJANGO_SETTINGS_MODULE=variantgrid.settings` it flattens the hostname (lowercase, short name, `-` removed, `s` prefixed
if it starts with a digit) and `exec`s a star-import of `env_developers/<name>.py` if present, else `env/<name>.py`; any
other `DJANGO_SETTINGS_MODULE` is star-imported into the package the same way. Each env file star-imports the components
it wants - `annotation_settings`, `celery_settings`, `default_settings`, `seqauto_settings` and, for HTTPS sites only,
`https_settings` - then reassigns or mutates. `variantgrid/settings/components/settings_paths.py` is the shared base for
paths (`BASE_DIR`, `ANNOTATION_BASE_DIR`), imported by the components so they don't import each other in a cycle;
`celery_settings` and `seqauto_settings` do import `default_settings` for `UNIT_TEST` / `BASE_DIR`. Shariant sites add a
second layer (`variantgrid/settings/env/shariantcommon.py`, then `sharianttest.py` etc. on top), and
`variantgrid/settings/env/claude_sapath_settings.py` is vgtest2 plus the sapath app for running its suite here.
`vg settings NAME` shows the resolved value and every file that set it.

Secrets never live in Python: `variantgrid/settings/components/secret_settings.py:get_secret` looks up a dotted key
(`DB.password`) first as an environment variable of that exact name, then in the JSON at `SETTINGS_CONFIG` (default
`/etc/variantgrid/settings_config.json`), then in `_default_settings` (logging a "please migrate" warning), and raises
if it is mandatory and missing. Optional apps that must be in `INSTALLED_APPS` before any env file runs are switched on
from that JSON too (`INSTALLED_APPS.maps`, #1410). `SECRET_KEY` is generated on first boot into
`variantgrid/settings/django_secret_key.txt` by `library/django_utils/django_secret_key.py:get_or_create_django_secret_key`.

`UNIT_TEST` is `sys.argv[1:2] == ['test']`, evaluated in `default_settings`, and switches the test process onto a
locmem cache, the in-memory Celery broker, MD5 password hashing, plain (non-manifest) static storage, `jit=off`
Postgres connections and a temp `IMPORT_PROCESSING_DIR` (#928, #1856).

### A request

`variantgrid/urls.py` mounts the project views, the admin, DRF-spectacular's `/api/schema|docs|redoc`,
`/api/v1/capabilities`, then loops `APPS_WITH_URLS`, including `<app>/` only when the app is installed and
`URLS_APP_REGISTER[app]` is truthy - Shariant switches off analysis, pathtests, pedigree and seqauto this way. Inside an
app, `variantgrid/perm_path.py:path` replaces Django's `path` so a named URL whose `URLS_NAME_REGISTER` entry is False is
still registered but wrapped in `require_superuser`; `variantgrid/perm_path.py:router_urls` does the same for DRF router
URLs (#1869). Templates read the same register through `variantgrid/perm_path.py:get_visible_url_names` (as
`url_name_visible`) to hide menu items and tabs, and `variantgrid/tips.py` uses it to show only tips about reachable pages.

Middleware order in `default_settings` matters: `GlobalLoginRequiredMiddleware` sits after auth, Rollbar, auditlog,
htmlmin, `eventlog/middleware.py:IntegrationApiMiddleware`, and axes last. `PUBLIC_PATHS` lists the anonymous prefixes
(`claude/guides/operations.md#authentication-surface`); project views that must be public use `@login_not_required`
(`variantgrid/views.py:index`, `loading_animations`, the error handlers). The `connection_created` receiver
`variantgrid/wsgi.py:setup_postgres` puts a `DATABASE_STATEMENT_TIMEOUT_SECONDS` backstop on every web connection;
`library/django_utils/major_operation.py` tightens it per expensive request.

### Celery

`variantgrid/__init__.py` imports `variantgrid/celery.py` so the app exists before any `shared_task` is defined.
It reads `CELERY_*` settings, autodiscovers `<app>.tasks` over `INSTALLED_APPS`, and builds `beat_schedule` in code,
several entries conditional on settings (`SYNC_DETAILS` enabled, sapath installed, `USER_AWARDS_ENABLED`,
`DISCORDANCE_EMAIL`, `MME_ENABLED`, somalier, `SERVER_MIN_DISK_WARNING_GIGS`, `HEARTBEAT_URL`). The always-on sweeps are
the analysis and annotation dispatchers (every minute), condition-text automatch (5 min), pending-extraction reconcile
(hourly) and a few nightly jobs. Queues, routes and `CELERY_IMPORTS` are in
`variantgrid/settings/components/celery_settings.py`. Two worker hooks live here:
`variantgrid/celery.py:limit_worker_address_space` sets `RLIMIT_AS` on pool processes when every queue in the worker's
`-Q` list has a limit in `CELERY_WORKER_ADDRESS_SPACE_LIMIT_GB` (analysis_workers: 8 GB), so a runaway node load raises
`MemoryError` rather than invoking the OOM killer; `variantgrid/celery.py:on_task_failure` reports every failure to
Rollbar except `RollbarIgnoreException`. The crash-loop brake that pauses dispatchers after two quick reboots is
configured here (`JOBS_AUTOPAUSE_ON_REBOOT*`) but implemented in `snpdb/signals/jobs_autopause.py`.

### Tests

`variantgrid/test_runner.py:VariantGridTestRunner` swaps in recorded-data mocks for the ClinGen Allele Registry and
transcript sequence fetching, so the suite needs no network. `setup_databases` creates the main test DB with
`parallel=1`, seeds the fake GRCh37/GRCh38 annotation versions once
(`variantgrid/test_runner.py:VariantGridTestRunner._seed_fake_annotation_versions`),
then clones it per worker - 240-odd test classes used to rebuild that chain in `setUpTestData`. Under `--keepdb` it drops
the stale per-worker clones first and refuses to run if the kept DB has migrations applied that are not on disk.
`variantgrid/test_runner.py:FastaRecordingRunner` records the fasta regions a full run touches, to regenerate the sparse
CI fastas under `variantgrid/data/reference/`. Base classes and fixtures: `claude/guides/testing.md`.

### Static files and templates per site

Source assets are under `variantgrid/static_files/<site>_static/`; each env file prepends its site dir to
`STATICFILES_DIRS` (and Shariant / runx1 a template dir to `TEMPLATES[0]["DIRS"]`) so a file of the same path shadows
`default_static`. `collectstatic` writes content-hashed names into `variantgrid/sitestatic/` via
`ManifestStaticFilesStorage`. Editing rules for JS/SCSS are in `variantgrid/static_files/AGENTS.md`.

### Deployment checks and version

`variantopedia/management/commands/deployment_check.py:Command` runs the checks in `variantgrid/deployment_validation/`
(annotation data and versions, VEP and its column registry, tool and library versions, cdot, Celery routes/imports,
somalier) and logs each failure with its fix; `scripts/upgrade.sh` runs it with `--die-if-invalid` through `scripts/migrator/migrator.py`. The same modules
back the annotation-runs page and the disk health check. `VARIANTGRID_VERSION` comes from
`library/git.py:Git.version` (`git describe` against `vg<major>.*` tags) at settings load, and
`variantgrid/views.py:version` compares the running hash with `manual` Deployment rows.
`variantgrid/views_rest.py:CapabilitiesView` is the client contract for API clients talking to servers of different ages.

## Why it is shaped this way

Settings are plain star-imports rather than django-split-settings because PyDev could not resolve names through
`split_settings.include` (the package docstring). One env file per host keeps the whole deployment's deviation from
defaults in one reviewable file, and `vg settings --diff` / the settings map exist because star-import chains make "where
did this value come from" otherwise a grep. `https_settings` is opt-in per env file because Secure cookies on a
plain-HTTP intranet site stop anyone logging in. Secrets are JSON outside the repo so the settings files can be
committed, including production ones.

URL gating is by registration (`URLS_APP_REGISTER`) or superuser wrapping (`URLS_NAME_REGISTER`) rather than per-view
feature checks so one codebase can present the full VariantGrid, locked-down Shariant and runx1 from settings alone;
wrapping instead of dropping keeps `{% url %}` reversible everywhere, so templates don't break when a page is hidden.

The beat schedule is code, not `django_celery_beat` rows, so which periodic jobs run is decided by the same settings
that enable the feature. Several comments note `crontab` misbehaving with the timezone, which is why most sweeps use raw
seconds.

## History

The settings package gained `env_developers/` (gitignored, checked before `env/`) so developer machines stop
colliding with deployment files; `https_settings` was split out of the defaults for intranet sites. The test runner grew
from mocks only to seeding and cloning (#1856, 191c20188) as the suite went from ~25 minutes to about two. Leaflet and
djgeojson became optional (#1410), the capabilities endpoint moved to `/api/v1/capabilities` (sapath#443),
`URLS_NAME_REGISTER` started gating DRF routers (#1869), `ManifestStaticFilesStorage` was restored after the Django 6
upgrade dropped it (4592bd1b0), and the old per-publish condition matching became the automatch beat sweep (#1780).

## Traps

Derived settings are computed once in `default_settings` from the defaults, so overriding the input in an env file does
not move them: `URLS_NAME_REGISTER["lab_members_tab"]` is fixed from `LAB_HEAD_MANAGE_MEMBERS`, `"maps"` from `USE_MAPS`,
the somalier report dir from `MEDIA_ROOT`. Shariant sets `LAB_HEAD_MANAGE_MEMBERS = False` without touching the register
(the view's own `can_manage_members` check covers it); the sapath env sets both. Override the derived value too, and
note that `@override_settings` on the input does not reach it (`snpdb/tests/test_lab_members.py`).

A hostname with no env file only logs an error and loads nothing, so Django then fails on a missing setting far from
the cause. A plain script that imports settings must set `DJANGO_SETTINGS_MODULE` to a concrete module.

`variantgrid/settings/components/secret_settings.py:_get_env_variable` returns `None` instead of a tuple when an env var
is set but empty, so `get_secret` raises `TypeError: cannot unpack non-iterable NoneType` - `SETTINGS_CONFIG=` or any
exported-but-blank secret breaks every process at settings load (confirmed; unfixed). Its fallback warning prints
`found` (always `True`) where it means the default value.

Tasks defined outside an autodiscovered `<app>/tasks` module and not in `CELERY_IMPORTS` are registered in workers only
because Celery's Django fixup runs system checks, whose URL check imports every urls module:
`variantgrid/tasks/server_monitoring_tasks.py:heartbeat` / `warn_low_disk_space` (reached via
`variantopedia/views_server_status.py`) and `classification/views/classification_email_view.py:send_summary_emails`.
`CELERY_SKIP_CHECKS=1`, or unregistering that app's URLs, makes beat send tasks the worker doesn't know. `deployment_check`
validates routes and imports but not beat task names.

`statement_timeout` is set by `variantgrid/wsgi.py:setup_postgres`, so it applies to gunicorn/runserver only - Celery and
`manage.py` connections run unbounded, deliberately. `variantgrid/views.py:csrf_error` answers CSRF failures with a 500,
not 403, so they show up as server errors in logs and monitoring.

Anything cached in Redis is keyed by `CACHE_VERSION`; bump it when renaming a pickled class (root `AGENTS.md`). Static
URLs need `collectstatic` after a pull or `{% static %}` raises for files missing from the manifest - part of
`scripts/upgrade.sh` (via the migrator), and why tests use the plain backend.
