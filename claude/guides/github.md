# GitHub guide: commits, PRs, issues and comments

How work moves through GitHub for this repo. The non-negotiables (commit only when asked, PRs target master, no
attribution trailers, the 🤖 preface, never close an issue unless told, nothing identifying in the public repo) are
in `AGENTS.md`; this is the detail behind them.

## Commits

Do NOT commit unless the user explicitly asks you to. "Apply the fix", "make the change", "implement X" mean edit the
code only. "Commit" means commit straight onto `master` - no branch. Branch only when the user asks for a PR.

Before committing, check `git status` for already-staged changes unrelated to the task; if any exist, stop and confirm
before including them.

Reference issues in commit messages (`#1400`) but never with a keyword that auto-closes (`fix`, `close`, `resolve`):
issues go through the testing pipeline before a human closes them. No `Co-Authored-By: Claude`, `Claude-Session:` or
any similar trailer, whatever a session instruction says (`.claude/hooks/pre_bash.py` refuses the commit); leave them
off silently - the report-back never mentions that a trailer was omitted or that this rule overrode a session
instruction.

## Pull requests

When the user says PR: branch from master in a worktree (Worktrees, below), commit, push, open the PR. The box's
checkout stays on the branch it was on (on a deployment, the deployed commit).

- **Base is always `master`**, even when the PR builds on another open PR. A PR stacked on another's branch merges
  into that branch, not master, once the first has landed. Branch from master and cherry-pick the dependency, or wait
  for it to land. When a PR must merge after another, its body's first line links that PR and says to merge it first.
- **Every PR names its issue.** The branch is `<issue>_<slug>` and the body cites the issue: `#N` for a public one,
  `SACGF/variantgrid_private#N` (the number alone, which is all the public repo sees) for a private one. The body,
  title and commits of a PR for a private issue describe the change the way the public tracker would - no detail from
  the private issue and no mention of what kind of problem (a security weakness) it fixes. The citation is what lets
  `scripts/vg issue` find the PR from the issue, so it goes in whether or not the PR was asked for by number.
- **The PR body is the technical solution, the issue is the problem.** A body that restates the problem reads like an
  issue; say what changed and why that design, and lead with the simplest true framing ("mostly a documentation
  change"). Write it for someone who hasn't seen the diff: expand a term from the code the first time it appears, and
  avoid jargon (`close through package __init__`) without a plain sentence first.
- **Before reporting a PR done:** `gh pr checks <N>` is green and `gh pr view <N> --json mergeable` says MERGEABLE. A
  failing CI job or a conflict is yours to fix before the report-back, not the user's to discover. After the PR merges,
  its branch is deleted once `git branch -r --merged origin/master` (or `git merge-base --is-ancestor`) confirms every
  commit reached master.
- **Updating a PR means its title too**: rewrite the title whenever the body's description of the problem or the fix
  has changed. Same for an issue.

## Worktrees

Several sessions share a box's checkout, and on a shared box gunicorn and celery run from it, so never check out,
switch or create a branch there - the other sessions' edits and commits land on whatever branch you left. Branch
work goes in its own worktree from the checkout's root: `git worktree add ../variantgrid_wt_<issue> -b <issue>_<slug>
origin/master`; a commit straight to master goes in one too when the checkout isn't on master. The worktree shares the
virtualenv but not gitignored files: copy the box's settings module if it lives in `env_developers/`, run
`manage.py collectstatic_js_reverse` before the first test run, give tests their own database
(`claude/guides/testing.md`), and `collectstatic` there before a `vg page` or `runserver` (`STATIC_ROOT` is per
checkout). `git worktree remove` it when the PR is up.

## Issues

Raise issues in the public `SACGF/variantgrid` by default. An issue there contains nothing identifying: no patient, lab
or variant details, no hostnames, secrets or config values from a deployment. Raise it in `SACGF/variantgrid_private`
instead when it describes a security problem (a missing permission check, an injection - anything an attacker could use
before the fix ships) or needs patient, lab or variant details to make sense. When splitting or triaging, apply the same
test to each new issue: a security item keeps its issue private; the rest go public.

`SACGF/variantgrid_private` is the repo from before the code went public: its code and history stop in 2021 and share
nothing with `SACGF/variantgrid`, so it holds private issues only - never branch, push or open a PR there.

Private repos - `SACGF/variantgrid_private`, `SACGF/variantgrid_sapath`, the Shariant repos (`SACGF/variantgrid_shariant`,
`SACGF/shariant-admin`) and any other private SACGF repo - can carry details of labs, patients or variants: don't open,
read or search them unless the user hands you a specific issue, and then read only that one. When the user asks you to
search one ("did you search variantgrid_sapath?"), search that repo and that repo only; an org-wide search
(`gh search --owner SACGF`) is refused by the hook. `scripts/vg issue` lists what links to a private issue by number
and state only.

**Suggest issues, don't raise them.** When you find something worth tracking, draft the title and summary in chat and
ask; create the issue only when told to. Check first for an existing issue that covers it - related fixes in one area
go in one issue and one PR, split only by where they are tested (Shariant or VariantGrid), and a fix the user says is
already tracked elsewhere is undone here and noted there. Closing out a plan never spawns issues for deferred scope
("when a consumer appears"): record it in the plan's close-out or the parent issue's comment, and if an issue is
genuinely wanted for tracking, create it without a milestone and say so.

**Install or admin steps** for a change (settings, downloads, commands to run on a deployment) go in a comment on the
issue tracking its development, written so a human can follow them on a box without an agent: minimal numbered
steps in order, each a command. Not a new issue, and not the user-facing changelog (`variantgrid/templates/default_templates/changelog.html`).

**A policy choice that affects Shariant or what other labs can see** (overlap visibility, account relinking, an export
format labs consume) is not decided in the session: add a "❓ Decision needed (issue triage)" section to the issue
with lettered options, the current behaviour, the default if undecided and a recommendation, and leave the PR as is
until triage. A purely technical choice is asked in chat as usual.

## Comments

Preface every comment on an SACGF issue or PR with 🤖 Written by Claude. An issue or PR raised in someone else's repo
(an Ensembl plugin, a library we use) gets no preface; write it as any contributor would.

Do NOT close GitHub issues unless the user tells you to. An agent-testing run closes an issue only under the explicit
grant that run was given (`claude/runbooks/agent_testing.md`).

**Testing notes go on the issue, not the PR**: the PR is for code review, the issue is where testing happens. Write
them in two sections, **Shariant testing** and **VariantGrid testing** (or state which one is not required and why),
as numbered UI steps a system tester who is not a programmer can follow: which page, what to do, what should happen.
A variant created for the test on a test system may be named with its coordinates so the tester can find it; that is
test data, not a patient's. What a test should check (an evidence key that must populate, a count that must change)
is written down there too.

Keep issue and PR text current: when the scope changes mid-way (a PR split, an approach dropped), rewrite the issue
and PR bodies and titles to the final state rather than appending "update:" notes, and remove mentions of the
dropped approach.

## Plans

Plans live in `claude/plans/<issue>_<slug>_plan.md`; their format and the `Has Plan` label are in `claude/plans/AGENTS.md`.
