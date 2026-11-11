# Migrations, manual operations and one-off commands

Rules for writing or editing a migration, a `ManualOperation` (a management command run on deploy) or a one-off
fix/backfill command, in any app. How the upgrader runs these steps is in `manual/__manual_readme.md`.

## Editing a pushed migration
By default, express a change of mind about a pushed migration (a field removed, a default changed) as a new migration on
top: renaming or editing one that a deployment has applied leaves `django_migrations` pointing at a name that no longer
exists (`InconsistentMigrationHistory`), repaired by hand on every database. A migration that only exists locally
(unpushed) can still be reshaped or regenerated - check with `git log origin/master -- <migration file>`. When I say a
pushed migration is safe to edit (it hasn't reached a deployment, or I'll fix the databases), edit it as asked.
Changing a library call inside a pushed migration (a renamed helper, an inlined one) is fine when the migration's behaviour
stays identical: migrations never freeze the names of the code they import.

## Management commands on deploy
A new management command that must run on existing deployments goes in a migration as a `ManualOperation`
(`manual/operations/manual_operations.py`, example the `ManualOperation.operation_manage` calls in
`snpdb/migrations/0002_seed_data.py`).

## One-off fixes and backfills
A one-off fix or backfill - something run once on existing data - is deleted at the next migration checkpoint
(`manual/__manual_readme.md#one-off-commands`), so make that obvious: name it `one_off_<what_it_fixes>`, declare
`category = "one-off"`, and keep the code that exists only for the fix (its queries, recounts, helpers) in the command
file and its test, so removing it leaves core code unchanged.

Calling existing core functions is fine. When the right fix needs something core doesn't have yet (a new parameter, an
importer that can write into an existing run), add it to core in the same PR: a capability core would keep is the
proper fix, not something to route around. Design the fix the way the existing code would: a one-off that needs the
data to hold a state it normally never does (overlapping ranges, a row whose type disagrees with its parent) is a
workaround - find the existing backfill or importer that does the same job and extend it.
