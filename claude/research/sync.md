# VariantGrid Sync App — Reference Document

## Purpose and Overview

The `sync` app manages bidirectional data synchronization between VariantGrid and external systems (Shariant and Alissa). It handles classification upload/download operations, tracks sync operations, and maintains audit records.

**Key Design:** Admin-configurable sync destinations; delta and full sync modes; per-classification sync records; health monitoring integration.

---

## Models

### SyncDestination
Configured sync connection to an external system.
- Fields: name (unique text), config (JSONField), enabled (bool, default True)
- Config JSON keys:
  - `type`: "shariant", "variantgrid", or "alissa"
  - `direction`: "upload" or "download"
  - `sync_details`: Key into settings.SYNC_DETAILS for credentials
  - `mapping`: Lab, share level, and user mappings (upload only)
  - `filters`: Optional JSON filters for classifications
  - `genome_build`: Target genome build
  - `exclude_labs` / `exclude_orgs`: Labs/orgs to exclude
  - `include_sources` / `exclude_sources`: Sources to include/exclude (Alissa)
- Methods: `run(full_sync=False, max_rows=None)`, `sync_details` property, `report()`, `get_reports()` (static)

### SyncRun
Single execution of a sync operation.
- Fields: destination FK, status (SyncStatus), full_sync (bool), max_rows (int, nullable), meta (JSON, nullable)
- TimeStampedModel (created, modified)

### ClassificationModificationSyncRecord
Tracks individual classification modification sync operations.
- Fields: run FK, classification_modification FK, success (bool), meta (JSON, nullable — includes remote ID)
- Methods: `filter_out_synced(qs, destination)` (static), `remote_url` property

### SyncStatus Enum
- `IN_PROGRESS = 'P'`, `SUCCESS = 'S'`, `NO_RECORDS = 'N'`, `FAILED = 'F'`

---

## Sync Runners

### Classification Upload (Shariant/VariantGrid)
1. Query last-published classifications at DISCORDANT_LEVEL_KEYS share level for configured labs
2. Apply optional JSON filters
3. Exclude already-synced records (unless full_sync=True)
4. Transform: convert evidence keys (HistoricalEKeyConverter), apply lab/share-level/user mappings, strip private fields (patient_id, dob, etc.)
5. Batch POST to remote system
6. Record sync records with remote IDs from response

### Classification Download (from Shariant/VariantGrid)
1. Query remote system for public classifications
2. Filter by genome_build, exclude_labs, exclude_orgs
3. Use `since` parameter for delta sync (last successful sync time)
4. Sanitize unknown evidence keys; convert to local format; build source_url reference
5. Bulk insert/update via BulkClassificationInserter

### Alissa Upload
1. Export classifications in MVL (MIAME-compliant) JSON format
2. Filter by genome_build and excluded/included sources
3. POST to Alissa system
4. Alissa handles import (CONTRIBUTE or MIRROR mode)

---

## Celery Tasks

**`sync_all()`** — Periodic task; iterates enabled SyncDestinations; calls destination.run(); exceptions logged but don't halt other syncs.

---

## Admin Interface

### SyncDestinationAdmin
Actions available:
- Run now (delta sync)
- Run now (delta sync) single row — for testing
- Run now (full sync)
- Validate configuration — verifies credentials in settings

### SyncRunAdmin
- List display: id, destination, created, status, full_sync, meta summary
- Filters by destination and status
- Download report action

### ClassificationModificationSyncRecordAdmin
- Searchable by classification lab_record_id
- Shows success/failure status and metadata

No direct URL patterns — admin-only operations.

---

## Health Check

`sync_health_check` signal receiver reports last successful sync for each destination. Warning threshold: 1 day since last success.

---

## Integration Points

| App | Integration |
|-----|-------------|
| Classification | Reads ClassificationModification; uses EvidenceKeyMap for field validation; uses BulkClassificationInserter for bulk creates |
| Library | ServerAuth for HTTP authentication; OAuth token management |
| Health Check | Reports last successful sync per destination |

---

## Settings

- `settings.SYNC_DETAILS` — Dict mapping sync_details keys to credential dicts
- Credential dicts consumed by `library.oauth.ServerAuth.for_sync_details()`
