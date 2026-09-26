# sync — research notes

Verified against a96540a68 on 2026-09-25

The sync app moves classifications between this VariantGrid and another one - in practice a lab's instance (SA
Pathology, variantgrid.com) pushing its shared records up to Shariant, and an instance pulling Shariant's public records
back down as external labs. A `sync/models/models.py:SyncDestination` is one configured direction to one remote, a
`sync/models/models.py:SyncRun` is one execution of it, and a
`sync/models/models_classification_sync.py:ClassificationModificationSyncRecord` says "this exact published version
reached that destination". There are no URLs of its own: runs happen hourly from celery beat or from admin actions, and
the only user-facing surface is the sync status panel on the classification page. Models, tasks and signals are in the
generated maps ([models](../maps/models.md#sync), [tasks](../maps/tasks.md#sync), [signals](../maps/signals.md));
`claude/domain.md` has the classification vocabulary (share level, lab, organisation, external lab).

## Flows

### Configuration and dispatch

A destination is a name, an `enabled` flag and a free-form `config` JSON edited in the admin. Credentials are never in
the row: `config["sync_details"]` is a key into `settings.SYNC_DETAILS`, which
`variantgrid/settings/components/default_settings.py:get_shariant_sync_secrets` builds from the `SYNC.<key>` secrets
(host, OAuth client, username/password). `sync/models/models.py:SyncDestination.sync_details` resolves it and
`library/oauth.py:ServerAuth.for_sync_details` turns it into an authenticated `requests` wrapper. The admin's
"Validate configuration" action (`sync/admin.py:SyncDestinationAdmin.validate_configuration`) checks the key resolves.

Which code handles a destination is decided by matching its config against a registry:
`sync/sync_runner.py:register_sync_runner` decorates a `sync/sync_runner.py:SyncRunner` subclass with required config
values (a set means "any of"), and `sync/sync_runner.py:sync_runner_for_destination` returns the first match. The
registry is filled by import side effect - `sync/apps.py:SyncConfig.ready` imports `sync.shariant`, whose `__init__`
star-imports the two runners - so a runner module that is not imported silently does not exist. Two runners are
registered, both for `type` in {shariant, variantgrid}: upload and download, selected by `direction`.

`sync/sync_run.py:run_sync` is the one entry point: it builds a `sync/sync_runner.py:SyncRunInstance` (destination,
`full_sync`, `max_rows`), calls the runner, and in a `finally` marks the SyncRun FAILED if the runner left it
IN_PROGRESS. The SyncRun row is created lazily by the `SyncRunInstance.sync_run` cached property, so a runner that
raises before touching it still gets a FAILED row from that `finally`. `sync/tasks/sync_tasks.py:sync_all` calls it for
every enabled destination; the beat entry in `variantgrid/celery.py` is only registered when some `SYNC_DETAILS` entry
is `enabled`. The admin runs the same function synchronously in the request (delta, delta with `max_rows=1` for a smoke
test, or full).

### Upload

`sync/shariant/variant_grid_upload.py:VariantGridUploadSyncer.records_to_sync` selects the
`is_last_published` ClassificationModifications at a share level in `ShareLevel.DISCORDANT_LEVEL_KEYS` whose lab
appears in `config["mapping"]["labs"]` (a lab not listed is never sent), narrows by the optional `config["filters"]`
through `sync/shariant/query_json_filter.py:QueryJsonFilter` (a small JSON-to-Q language over `published_evidence`,
with the legacy `somatic` key meaning an allele-origin test), holds back gene-level variants (fusions, copy number
events) unless `remote_gene_level` is set, and for a delta run drops every modification that already has a successful
sync record at this destination (`ClassificationModificationSyncRecord.filter_out_synced`). Because a new publish makes
a new last-published modification, "changed since last sync" falls out of that exclusion with no timestamps involved.

`VariantGridUploadSyncer.classification_to_json` turns each record into a v2 API upsert: the id is the mapped lab plus
lab record id (`classification/models/classification_ref.py:ClassificationRef.make_lab_id_str`), the share level and
owner are mapped through `config["mapping"]`, `SHARIANT_PRIVATE_FIELDS` (patient id, DOB, sample id, family id, REDCap
id ...) are removed, every other known evidence key is sent (as None if absent, so a cleared value clears remotely), and
`sync/shariant/historical_ekey_converter.py:HistoricalEKeyConverter` rewrites keys an older local instance still uses.
A withdrawn classification is sent with `delete: True`. Batches of 50 go to the remote's v2 classification record API
with an `import_id` of this host's name and `status: complete` on the final batch, so the remote's
ClassificationImportRun knows when the upload is finished. A connection error or timeout is retried in place (the
records are upserts, so re-sending a batch the remote already processed is harmless); an HTTP error fails the run with
earlier batches already recorded. Each result in the response, returned in request order, becomes a sync record whose
`meta` holds the remote JSON, including the remote pk that `ClassificationModificationSyncRecord.remote_url` links to.

### Download

`sync/shariant/variant_grid_download.py:VariantGridDownloadSyncer` fetches the remote's public export
(the classification export API as JSON, one genome build, minus `exclude_labs` / `exclude_orgs`) in one
response and streams it with ijson. A delta run passes `since` = the remote's `Last-Modified` header stored in the meta
of the last SUCCESS run (`SyncRunInstance.last_success_server_date`) - the remote's clock, not ours. Each record's data
loses evidence keys this instance does not know (noted on a `source_url` value pointing back at the remote record), and
records are upserted in batches of 50 through `classification/models/classification_inserter.py:BulkClassificationInserter`
as the admin bot with `force_publish`, sleeping 10 s between batches. A record for a lab that exists locally and is not
`external` is skipped - the guard against importing Shariant's copy over our own records if `exclude_labs` is wrong - and
an unknown lab is created as an external Lab (and Organization, country Australia) with a `report_message` so someone
notices. `max_rows` is rejected.

### Reporting

`sync/classification_sync_status.py:classification_sync_status` is what the classification page shows
(`classification/views/views.py` passes it as `sync_statuses`): for each enabled upload destination whose runner is a
`sync/sync_runner.py:ClassificationUploadSyncRunner`, it reports Uploaded / changes pending / pending / held locally,
using the sync records for what actually arrived and the runner's own `records_to_sync(full_sync=True)` and
`exclusion_reasons` for what would happen next, so the explanation cannot drift from the selection. The two signal
receivers, `sync/signals/sync_health_check.py:sync_health_check` and
`sync/signals/sync_integration_status.py:sync_integration_status`, put each destination's last success on the health
check and the server status integrations table with a one-day warning age; both count NO_RECORDS as a success.

## Why it is shaped this way

The config is JSON rather than columns because every deployment's mapping is different (which labs go up, legacy lab
names mapped onto current ones, share levels downgraded on the way out, local usernames renamed) and the set of keys has
grown feature flag by feature flag. `remote_lab_record_url` and `remote_gene_level` are compatibility gates for a remote
that has not been upgraded yet: VariantGrid deployments upgrade on their own schedules (SA Pathology runs validated
releases), so the sender must not assume the receiver understands a lab-record URL or a gene-level variant. The same
concern is behind `HistoricalEKeyConverter`, whose conversion table is now down to dropping `variant_type`; its
commented-out entries are the record of earlier evidence-key renames.

Sync records are only written once the connection has worked (see the model docstring): a bad password should produce
one FAILED run, not a failure row per classification. Tracking per modification rather than per classification is what
makes the delta query a single exclude and what lets the status panel tell "uploaded" from "uploaded, later changes
pending". The runner registry exists so a destination type can be added without touching dispatch - it was built when
Alissa was a second target.

## History

The app was created in 2020 (`sync/migrations/0001_initial.py`) to upload to Shariant, then gained download from
Shariant (2021) and Alissa upload/download. Credentials moved from a single `SYNC` secret to keyed `SYNC_DETAILS` so one
instance can talk to several remotes (`sync/migrations/0006_one_off_manual_modify_config.py` is the deploy step, and
`get_shariant_sync_secrets` still refuses the old shape); SyncRun gained `full_sync` / `max_rows` in 2023. In 2026 the
upload runner moved to the v2 record API with the lab-record URL gate (SACGF/variantgrid_sapath#427), the
classification-page status panel landed (#1347), sync-created external labs started notifying, uploads got a longer
timeout and connection retries, and gene-level variants got the `remote_gene_level` gate (#1506, #1836). Alissa and the
MVL export it used were removed in September 2026 (variantgrid_private#3809); `sync/alissa/` may survive on disk as an
untracked `__pycache__` only.

## Traps

A per-record failure on the remote is recorded as a success. The v2 record endpoint answers a record it could not
process with `{"internal_error": ...}` (`classification/models/classification_utils.py:ClassificationPatchResponse.to_json`)
inside an HTTP 200, and `VariantGridUploadSyncer.sync` creates every sync record with the default `success=True`
without looking at the result, so that modification is excluded from every later delta run and the status panel says
"Uploaded". Only a full sync or a new publish re-sends it.

Withdrawal does not reach the remote on a delta run. `classification/models/classification.py:Classification.set_withdrawn`
flips a flag on the Classification and creates no new modification, so the last-published modification that was
already synced stays excluded by `filter_out_synced` and the `delete: True` branch of `classification_to_json` only fires
on a full sync or after a later publish. Un-withdrawing is never propagated (the code says so).

The April 2026 hardening of this app (53648bdbe, variantgrid_private#3831: filter key validation in `QueryJsonFilter`,
lab group name / `exclude_labs` validation before download creates Labs, `remote_pk` validation in `remote_url`, no
config in the `sync_runner_for_destination` error, URL scheme checks in `library/oauth.py`) is not in the current code:
the merge of master into its branch (4beb94b43, 2026-08-07) resolved every sync and oauth conflict in master's favour.

`sync/tasks/sync_tasks.py:sync_all` runs destinations one after another in a single task on the default `db_workers`
queue (no `queue=`), and the admin actions run a sync inside the web request, so a large full sync belongs in a shell or
celery, not the admin. The download holds the
whole export in memory (`stream=False`, then ijson over `response.content`).

The download's `since` comes only from the last SUCCESS run; a NO_RECORDS run stores `server_date` too but is ignored,
so quiet periods make each delta re-fetch from further back (harmless, as inserts are upserts). Conversely
`sync/models/models.py:SyncDestination.report` treats only SUCCESS as success, unlike the health and integration
receivers; it and `get_reports` have no callers. `SyncRunner.report_on` returns None for both remaining runners, so the
SyncRun admin's "Download report" action does nothing now that the Alissa report is gone.

Deleting a SyncRun cascades to its sync records, so the next delta run re-sends those records (upserts, so safe).
