# Operations guide

The L4 facts from [agent_system.md](../plans/agent_system.md) §4.4: what is running where, how a deploy happens, and the
scale that shapes how a change must be written. Live values come from `python3 manage.py vg status`; this file holds
what the command cannot know. Verified against vg-test2 on 2026-09-06.

## Deployments

| Deployment | Hostname → settings file | Notes |
|---|---|---|
| vg-test2 (test.variantgrid.com) | `variantgrid/settings/env/vgtest2.py` | The lab box this repo is usually driven from; see "This box" in `CLAUDE.md` |
| variantgrid.com | `variantgrid/settings/env/vgaws.py` | Public instance on AWS |
| Shariant (test / demo / prod) | `variantgrid/settings/env/shariantcommon.py` + `sharianttest.py` / `shariantdemo.py` / `shariant.py` / `shariantsecurity.py` alongside it | Australian classification sharing; patients and analysis URLs unregistered |
| SA Pathology | private repo `variantgrid_sapath` (settings and site-specific apps live there) | Largest production data; clinical use |
| runx1db | `variantgrid/settings/env/runx1db2.py` | Gene-specific public database |
| CI | `variantgrid/settings/env/github_actions.py` with `config/ci/settings_config.json` | Also the canonical module for `vg map` |

Settings resolve by hostname: `variantgrid/settings/__init__.py` lowercases the short hostname, strips `-`, and imports
`env_developers/<name>.py` if it exists, else `env/<name>.py`. Each env file star-imports the components (`default_settings`,
`celery_settings`, `annotation_settings`, `seqauto_settings`, `https_settings`) then overrides. `python3 manage.py vg settings NAME`
shows the resolved value and every file that assigned it; `vg settings --diff` lists what this box's env file changes; a plain
script must set `DJANGO_SETTINGS_MODULE` itself (`variantgrid.settings.env.vgtest2`).

Secrets (database, broker, API tokens) come from `/etc/variantgrid/settings_config.json` via
`variantgrid/settings/components/secret_settings.py:get_secret`, path overridable with `SETTINGS_CONFIG`;
`config/settings_config.json` is the template and `config/ci/settings_config.json` the CI copy.

## Services, queues, logs

systemd units, all restarted by `scripts/restart_services.sh` (`stop_services.sh` / `start_services.sh` likewise):
`gunicorn` (:8000, behind nginx/https on deployments), `celeryd_beat`, and one `celeryd_<queue>` per queue in
`variantgrid/settings/components/celery_settings.py:CELERY_TASK_QUEUES`:

| Queue | Purpose |
|---|---|
| `db_workers` (default) | database-only tasks; could run on another machine |
| `web_workers` | tasks that need the web server filesystem (uploads, images, reports) |
| `analysis_workers` | node updates and analysis runs |
| `annotation_workers` | VEP dumps / runs / uploads |
| `variant_id_single_worker` | one process: inserts new Loci/Variants so there is exactly one row per coordinate |
| `scheduling_single_worker` | one process: schedulers that must not race (annotation, node tasks) |

Route a task with `@celery.shared_task(queue=...)` or an entry in `CELERY_TASK_ROUTES`; `claude/maps/tasks.md` lists every
task with its queue and enqueuers. The broker is RabbitMQ (`CELERY.broker_url` secret), the result backend and cache Redis.
`vg status` shows per-queue depth and consumer count from a passive declare.

Logs: `/var/log/variantgrid/` (`gunicorn.log`, `celery*.log`, `celeryd_beat.log`). Application events are also rows in
`eventlog/models.py:Event` (severity, app, name, details) - `vg status` lists the last ten ERRORs; `ViewEvent` records page views
(the `claude_agent` user's included). Uncaught exceptions go to Rollbar when `ROLLBAR.access_token` is configured, and Slack
notifications through `library/log_utils.py:AdminNotificationBuilder` when `SLACK.enabled`.

## Deploy and upgrade

`scripts/upgrade.sh <target>` on a deployment: `install_requirements.sh` (uv-managed `.venv` from `requirements.txt`), then
`scripts/migrator/migrator.py`, which pulls, runs `manage.py migrate`, and reads `manage.py manual_outstanding` (JSON) to
surface the deploy-time steps migrations registered with `manual/operations/manual_operations.py:ManualOperation`.
`manage` steps auto-run once their `requires` gates (`manual/gates.py`) clear and the command still exists; `other` steps are
printed for a human. Completion is a `ManualMigrationAttempt`. `scripts/deployed.sh` records the deploy in Rollbar and runs
`manage.py deployed`. `scripts/restart_services.sh` finishes it.

Pushed migrations are frozen (`CLAUDE.md`); a data fix that must run on every deployment is a `ManualOperation` in a
migration, not a note in a PR. Annotation upgrades (new VEP, new columns) are their own procedure: a new
`VariantAnnotationVersion` per build via `create_new_variant_annotation_version`, then re-annotation of every variant in
`AnnotationRun` batches - hours on the big deployments.

Caches: `CACHE_VERSION` in `default_settings.py` is the Redis key version - bump it after a change that makes cached values
wrong. Process-level caches (`library/cache.py:timed_cache`, the caching model managers) only clear on restart.

Static files: `manage.py collectstatic` writes `variantgrid/sitestatic/`; storage is Django's plain `StaticFilesStorage`
(no manifest hashing), which is why CI runs tests without collectstatic and `library/django_utils/unittest_utils.py`
overrides `STORAGES` to the plain backend.

## Data roots

| Setting | vg-test2 | Holds |
|---|---|---|
| `PRIVATE_DATA_ROOT` | `/opt/variantgrid/data` | uploads (`UPLOAD_DIR`), per-deployment private files; gitignored |
| `ANNOTATION_BASE_DIR` | `/data/annotation` | VEP code + caches (`ANNOTATION_VEP_BASE_DIR`), fastas, chain files, scratch, somalier |
| `MEDIA_ROOT` | `/opt/variantgrid/media_root` | user media |
| partition dumps | `/data/database/partition_dumps/<db>/` | pg_dump archives written before a partition is dropped (#1537) |

`vg status` reports free space on each. The root filesystem on vg-test2 is small (about 30 GB) and usually above 85 %;
large scratch output belongs under `/data`.

## Scale <a id="scale"></a>

The numbers that decide whether a design works. vg-test2 today (from `vg status`): snpdb 175 GB, 16.7 M variants, 15.8 M
loci, 53 k alleles, 211 samples, 726 classifications, 46 analyses, 200 k variant tags. Production deployments are larger;
SA Pathology holds on the order of 400,000 `VariantTag` rows (a known artefact is re-tagged in every analysis it appears in)
and many more variants and samples.

- **Whole-database operations are batched by default.** An unbatched liftover run held every Allele, VariantCoordinate and
  AlleleLiftover in memory and took the machine down. Page the queryset by pk range (not `.iterator()`, which holds a cursor
  open for hours), size the batch from a settings constant, and fan batches out as separate celery tasks so they run in
  parallel and one failure only loses its own batch: `snpdb/liftover.py:create_liftover_pipelines` and
  `snpdb/tasks/liftover_tasks.py:liftover_allele_batch` (`settings.LIFTOVER_BATCH_SIZE`) are the pattern. Two reasons, both
  weighed by the owner: cap resource usage, and limit how much one error destroys.
- **Aggregate in SQL, not Python.** Grouping, deduplicating or merging tags, genotypes or annotations belongs in
  `annotate()` / `Count` / `Min` + `update()` / `delete()` on a queryset. Streaming rows into a `set`, `dict` or `Counter` is a
  design error at these sizes, in management commands as much as in views.
- **Restrict variants to a build with an IN list on contig ids** (`Variant.get_contigs_q`); joining through GenomeBuildContig
  wrecks the planner's row estimate (#1720).
- **Partitioned tables** - one child table per collection or version, dropped rather than deleted from:
  `CohortGenotypeCollection`, `VariantCollection`, `VariantZygosityCountCollection`
  (`library/django_utils/django_partition.py:RelatedModelsPartitionModel`) and every `SubVersionPartition` annotation table
  (`annotation_variantannotation` per VariantAnnotationVersion). Before a drop, the archive pipeline writes a pg_dump -
  restore per `claude/runbooks/restore_partition_archive.md`. Bracket any query against a partition with `temporary_db_table`.
- **Long requests** go through `library/django_utils/major_operation.py:major_operation` (per-user concurrency cap, lower
  `statement_timeout`); grids inherit `MajorOperationViewMixin`. Postgres JIT is off (`-c jit=off`).
- **`EXPLAIN` on real data is what this box uniquely offers**: `manage.py profile_analysis_nodes --analysis <id> --rerun --explain`
  for analysis nodes; for anything else, `queryset.explain()` in `manage.py shell` inside a read-only transaction.

## Testing pipeline

CI (`.github/workflows/django-tests.yml`) runs the suite with `--parallel 4 --keepdb` against Postgres 16, Redis and
RabbitMQ service containers under `variantgrid/settings/env/github_actions.py`, skipping pushes that only touch `*.md` or `claude/**`;
`.github/workflows/agent-maps.yml` runs `vg map --check` on every push without a database. Browser regression tests live in
the private [variantgrid_autotests](https://github.com/SACGF/variantgrid_autotests) repo (Selenium, run with
`run_tests.py <instance.ini> [test | +keyword | -keyword]` against a deployed instance) - the after-deploy net, not the edit
loop. GitHub issues are closed by a human after that pipeline, never by a commit keyword.

## Authentication surface

`global_login_required.GlobalLoginRequiredMiddleware` requires login everywhere except `PUBLIC_PATHS`
(`variantgrid/settings/components/default_settings.py`), which exempts the API prefixes (`/classification/api/`, `/patients/api/`,
`/seqauto/api/`, `/upload/api/`, `/mme/api/`, `/beacon/`) so DRF's `IsAuthenticated` can answer 401 instead. A plain Django
`View` mounted under an exempt prefix is reachable anonymously (it usually 500s on `AnonymousUser`); new endpoints there must
be `rest_framework.views.APIView` subclasses. Verify with an unauthenticated request - a 500 means the view is unprotected.
