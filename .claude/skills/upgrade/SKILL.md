---
name: upgrade
description: Pull and upgrade this VariantGrid deployment - "pull and upgrade", "pull and update", "update the server", "upgrade and restart". Runs scripts/upgrade.sh --quick, runs the outstanding manual operations through the upgrader so they are recorded as done, restarts the services the box allows, and checks the result with vg status.
---

# Upgrade this deployment

The whole sequence from `claude/guides/operations.md#deploy-and-upgrade`, so a "pull and upgrade" never stops at the
`git pull`: migrate, collectstatic and the deploy-time commands are the part a bare pull skips.

1. **Which box.** `hostname`. A production host (`AGENTS.md`, This box) is not upgraded by this skill unless the user
   asked for a production rollout; a shared test box (vg-test2, the Shariant test boxes) is what this is for.
   Record the commit before the pull: `git rev-parse --short HEAD`.
2. **Upgrade.** From the repo root:
   ```bash
   source .venv/bin/activate && scripts/upgrade.sh --quick < /dev/null
   ```
   Exit 0 means every standard step passed (pull, migrate, collectstatic, deployment_check). Exit 1 at the menu
   ("No input available for the menu") means a step failed or manual steps are outstanding: read the output.
3. **Outstanding manual operations**, through the upgrader rather than as bare `manage.py` commands, so each one is
   recorded as a `ManualMigrationAttempt` and the next upgrade doesn't ask for it again:
   ```bash
   scripts/upgrade.sh --auto-manage < /dev/null      # runs the unblocked 'manage' steps, carries on past a failure
   python3 manage.py manual_outstanding              # what is left: 'other' steps are for a human
   ```
   Anything still outstanding goes in the report-back with the step's note (`manual/__manual_readme.md` says how a
   human records an `other` step done). An annotation backfill or VEP run is never started by this skill: report it.
4. **Restart**, even when the pull brought no code: collectstatic rewrote the hashed static names and running gunicorn
   keeps the manifest it loaded at start. On vg-test2 sudoers allows exactly one unit per call, named without
   `.service`:
   ```bash
   for unit in gunicorn $(systemctl list-units 'celeryd_*' --no-legend --plain | awk '{print $1}' | sed 's/\.service$//'); do
     sudo -n systemctl restart "$unit"
   done
   ```
   The pre_bash hook lets these through on vg-test2 without a prompt; another shared box prompts, which is the
   AGENTS.md rule for it.
5. **Check.** `python3 manage.py vg status` (every service active, queues consuming, 0 manual tasks outstanding, no
   new ERROR events from the restart) and `python3 manage.py vg page /variantopedia/dashboard`.
6. **Report**: the commits pulled (`git log --oneline <before>..HEAD | head -20`), the manual operations run and their
   outcome, anything outstanding for a human with its note, and the `vg status` summary. One screen.
