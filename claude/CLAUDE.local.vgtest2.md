# CLAUDE.local.md for vg-test2

Copy this file to the repo root as CLAUDE.local.md (gitignored) on vg-test2. It is the per-host section that used to
live in `CLAUDE.md`; the checked-in file is the same on every box, so anything true of one machine only lives here.

## This box

`vg-test2` (test.variantgrid.com) is a shared lab: gunicorn and the celery workers run against a 175 GB database that human
testers are also using, so what you change they see. `python3 manage.py vg status` is the first thing to run in a session.

Safe without asking: anything read-only (`git`, `vg *`, `gh issue/pr view`, `manage.py shell` that only reads, `EXPLAIN`),
rendering pages as `claude_agent` with `vg page`, and tests with `--keepdb` (they use `test_snpdb`).
Ask first: restarting or stopping services, `manage.py migrate`, creating or deleting annotation versions or running VEP,
liftover across the database, and any write to `snpdb_variant`, `snpdb_allele` or `annotation_variantannotation`. The
`.claude/hooks/pre_bash.py` hook turns those into a confirmation prompt with the reason.

Row counts, data roots, logs and the deploy procedure are in `claude/guides/operations.md`.
