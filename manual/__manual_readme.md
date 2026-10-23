# Manual

The manual app is to provide support for more human interaction required upgrade steps, as opposed to automatic model updates.

This allows us to have upgrade steps that can't easily be automated, such as updating external files, packages etc.

## Manual migration tasks

Migrations register manual steps with `ManualOperation` (see `operations/manual_operations.py`):

- `operation_manage([...])` — a `python3 manage.py ...` command (category `manage`).
- `operation_other([...])` — a free-text human step (category `other`).

Each registers a `ManualMigrationTask` (PK = the command string) plus a `ManualMigrationRequired`
row. `manual/upgrader.py:outstanding_tasks_json` lists everything still outstanding (`manage.py manual_outstanding`
prints it), which the upgrader (`scripts/upgrade.sh` -> `manage.py upgrade`) surfaces. Completion is recorded as a
`ManualMigrationAttempt`.

## Dependency gates (auto-running manage steps)

`manage` steps can be auto-run by the upgrader, but some must not run until a prerequisite is met
(e.g. ontology/annotation upgraded, transcripts in current cdot format). A step declares this at the
call site:

```python
ManualOperation.operation_manage(["match_patient_phenotypes", "--stale"], requires=["ontology-imported"])
```

(`--stale` rematches only the phenotype sentences matched with an older `PHENOTYPE_MATCHER_VERSION` or ontology;
each matcher version bump re-registers it, and one successful run satisfies every registration before it.)

`requires` is persisted on `ManualMigrationTask.requires`. Gate *definitions* live in `gates.py`,
keyed by gate name (not command):

- **AutoGate** — a predicate over live models (e.g. `ontology-imported`, `cdot-current`).
- **ManualGate** — confirmed once by an operator via `manage.py manual_gate --satisfy <name>`
  (e.g. `variant-annotation-current`).
- **`after:<task_id>`** — depends on another manual task completing.

Each outstanding task carries `requires` / `blocked_by` / `runnable` / `command_exists`. The
upgrader's auto-manage pass (`am` menu option or `upgrade.sh --auto-manage`) runs every unblocked
`manage` task, re-evaluating between passes so `after:` gates unblock as prerequisites finish. It carries
on past a failure without retrying it: the failed task stays outstanding, so anything gated `after:` it
stays blocked. Blocked steps, obsolete steps (command no longer exists), and `other` steps are reported
for a human rather than run.

## The upgrader

`scripts/upgrade.sh` installs requirements, then execs `manage.py upgrade`
(`manual/upgrader.py:Upgrader`), which runs every step in that one process with `call_command` - a
Django start per step was most of the upgrade's time. The menu takes several keys or ranges at once
(`m c 2-5`, or `upgrade.sh --steps m,c,2-5`), runs them in order and carries on past a failure, skipping a
step gated `after:` one that failed. `--quick` runs the standard steps and stops at the first failure.
A git pull that moves HEAD installs the new requirements and re-execs the process with the rest of the
selection, so nothing after the pull runs on the old code.

Use `manage.py manual_gate` to list gate status.

## Obsolete tasks

A `manage` task whose command has since been deleted is never auto-run (it would just crash) — it's
reported as obsolete for review, where selecting it offers to mark it complete.

## The migration squash

Every app's migrations were squashed into fresh initials and a `seed_data` migration (#736); an older database
upgrades through branch `migration_checkpoint` first. `squashed_history.py` holds the guard and the reset - see the
paragraph in `claude/guides/operations.md#deploy-and-upgrade`.

When something breaks in a way the squash could explain - a seed row or index a deployment lacks, a data fix that
never ran there, a schema that differs from a fresh build - start from #736. Everything the squash removed is still
in git at tag `pre_squash_checkpoint` (the tip of branch `migration_checkpoint`, kept forever - a GitHub ruleset
blocks deleting or moving it): the old per-app migrations (with their `RunPython` data fixes and `RunSQL`),
the `ManualOperation`s they registered, and the one-off management commands those ran (deleted from master in
593edc560 and 652c2fc3d). Read one with `git show pre_squash_checkpoint:<path>`, list an app's old migrations with
`git ls-tree --name-only pre_squash_checkpoint snpdb/migrations/`, and find which migration registered a command with
`git grep <command> pre_squash_checkpoint -- '*/migrations/*'`. To re-run a deleted command, bring it back as a new
command (or a new migration's `ManualOperation`) on master rather than checking the tag out, since the
checkpoint's system check stops every command there once its last migration is applied.
