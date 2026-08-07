# Manual

The manual app is to provide support for more human interaction required upgrade steps, as opposed to automatic model updates.

This allows us to have upgrade steps that can't easily be automated, such as updating external files, packages etc.

## Manual migration tasks

Migrations register manual steps with `ManualOperation` (see `operations/manual_operations.py`):

- `operation_manage([...])` — a `python3 manage.py ...` command (category `manage`).
- `operation_other([...])` — a free-text human step (category `other`).

Each registers a `ManualMigrationTask` (PK = the command string) plus a `ManualMigrationRequired`
row. `manage.py manual_outstanding` reports everything still outstanding as JSON, which the upgrader
(`scripts/migrator/migrator.py`) surfaces. Completion is recorded as a `ManualMigrationAttempt`.

## Dependency gates (auto-running manage steps)

`manage` steps can be auto-run by the migrator, but some must not run until a prerequisite is met
(e.g. ontology/annotation upgraded, transcripts in current cdot format). A step declares this at the
call site:

```python
ManualOperation.operation_manage(["match_patient_phenotypes", "--clear"], requires=["ontology-imported"])
```

`requires` is persisted on `ManualMigrationTask.requires`. Gate *definitions* live in `gates.py`,
keyed by gate name (not command):

- **AutoGate** — a predicate over live models (e.g. `ontology-imported`, `cdot-current`).
- **ManualGate** — confirmed once by an operator via `manage.py manual_gate --satisfy <name>`
  (e.g. `variant-annotation-current`).
- **`after:<task_id>`** — depends on another manual task completing.

`manual_outstanding` emits `requires` / `blocked_by` / `runnable` / `command_exists` per task. The
migrator's auto-manage pass (`am` menu option or `migrator.py --auto-manage`) runs every unblocked
`manage` task, re-evaluating between passes so `after:` gates unblock as prerequisites finish, and
stopping on the first failure. Blocked steps, obsolete steps (command no longer exists), and `other`
steps are reported for a human rather than run.

Use `manage.py manual_gate` to list gate status.

## Obsolete tasks

A `manage` task whose command has since been deleted is never auto-run (it would just crash) — it's
reported as obsolete for review. `migrations/0004_complete_obsolete_manual_tasks.py` marks the known
obsolete task ids complete so they drop out of the upgrade flow.
