# Restore a partition archive

The pre-drop archival pipeline (#1537) writes a `pg_dump --format=custom` of a
`RelatedModelsPartitionModel`'s child partition tables before dropping them.
This runbook restores one of those dumps.

The Django admin → snpdb → PartitionArchive list view shows every archive ever
written, its status, sha256, and `dump_path`.

## 1. Confirm dump integrity

```bash
sha256sum /data/database/partition_dumps/<db>/<archive_name>.dump
# compare to PartitionArchive.sha256 for archive_pk=<pk>
```

## 2. Restore

The dump is in custom format; restore with `pg_restore`:

```bash
pg_restore --dbname=variantgrid --jobs=4 \
    /data/database/partition_dumps/<db>/<archive_name>.dump
```

`pg_restore` and `pg_dump` must use the same Postgres major version (or one
major newer than the dump). Run `pg_dump --version` on the box that produced
the dump and compare to the restore host.

## 3. Re-attach restored child tables to the parent

The dump captures plain CREATE TABLE statements for each child partition.
After restore the children exist as standalone tables; re-attach them to the
inheritance parent (post-#1534 the project will move to native partitioning, in
which case use `ALTER TABLE ... ATTACH PARTITION` instead):

```bash
psql -d variantgrid -c '
    ALTER TABLE "<child_table>" INHERIT "<parent_table>";
'
```

The child / parent names live on the `PartitionArchive` row in
`source_table_names` (parents) and can be derived per-PK via
`source.get_partition_table(base_table_name=t)`.

## 4. Mark the archive restored in Django admin

Admin → snpdb → PartitionArchive → select rows → "Mark restored". This clears
`data_archived_*` on the source row so the model's read paths stop raising
`DataArchivedError`.

## Schema renames between archive and restore

If `RECORDS_BASE_TABLE_NAMES` for the source model has changed since the
archive was written, the restore lays down tables under the names captured in
the dump. Rename them with `ALTER TABLE ... RENAME TO` before re-attaching to
the (renamed) parent.
