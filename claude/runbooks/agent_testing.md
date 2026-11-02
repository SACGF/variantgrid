# Agent testing: working a ZenHub test pipeline on a shared test box

The user asks for this as "agent-test the issues in the VariantGrid test column" or hands over a list of issues. The
box is a shared test deployment (vg-test2 unless said otherwise); the rules for it are `AGENTS.md`, This box, and
`claude/runbooks/zenhub_issue_triage.md` has the pipelines, the ZenHub queries and where an issue moves next.

## Scope

- **The issues named, or the ones currently in the pipeline named** - not every issue touched recently. Issues already
  tested and closed are not re-opened; an issue already agent-tested (a previous 🤖 comment with a result) is
  re-tested only if a PR or commit citing it has landed since that comment (`scripts/vg issue N` lists them with dates).
- A private issue (`SACGF/variantgrid_private`, `SACGF/variantgrid_sapath`): read only that issue, and keep its
  details out of anything public - a public comment on it is fine, a comment on a public issue that repeats it is not.

## What the user grants for the run

Two things need an explicit grant in the session, and the grant covers only the issues it names:

- **Closing issues.** Humans close issues; an agent closes one only when the user said so for this run ("I'm
  authorising you to close issues in this session"). Without that, the agent comments and leaves it open with a
  recommendation.
- **The writes AGENTS.md asks about on this box**: VCF uploads and the variant / allele / annotation writes and VEP
  runs they trigger, creating a variant from HGVS, zygosity-count changes, subdividing an annotation run,
  `vg page --js-writes-to-database`. When the pre_bash hook prompts for one that the grant covers, go ahead; for one
  it doesn't, ask.

Fake and test data is fine without a grant and stays afterwards.

## Sharing the box

Another agent or a human tester may be on the box at the same time.

- Name everything you create `agent_test_<issue>_...` (users, labs, VCFs, analyses, gene lists, tags) so runs stay
  independent, and leave it in place: the data is the evidence a tester can re-check.
- Make a fake user who owns a fake lab for the uploads, analyses and tags; use a superuser only where a step needs one.
- Judge a step by the change it causes in your own data (counts before and after on *your* VCF), not against a snapshot
  of the box, so other imports running at the same time can't affect the result.
- Order the tests so later ones build on earlier ones (upload first, then the variant created from HGVS, then
  analyses on the uploaded VCFs, then anything that subdivides or reruns annotation last).
- Benchmarks and bulk loads stay small (thousands of rows, minutes): the box has four cores and human testers.

## Per issue

1. `scripts/vg issue N`: the issue, every comment and the PRs. Read the merged PR's diff. The latest "Moved to
   VariantGrid test" / testing comment holds the test steps; a checklist in the body is the test plan.
2. Confirm the fix is on this box: `git merge-base --is-ancestor <merge sha> HEAD` in the deployed checkout. If it
   isn't, stop on that issue and say so.
3. Test it the way a user would: the page (`vg page <url> --queries`, and `--js-writes-to-database --click '#tab'`
   under the grant for JS), the API, or the management command, plus targeted `--keepdb` tests where they add
   evidence. For markup-escaping checks, put harmless markup in your own objects' names (`<i>agent_test</i> & <x>`)
   and check every listed grid shows it literally, once.
4. Comment on the issue, starting with `🤖 Written by Claude`: what changed (PR number), how it was tested and the
   result, and what still needs doing. For Shariant testing, numbered UI steps a tester can follow.
5. Outcome:
   - **Close** only if it passed, nothing further is needed, and the run has the closing grant.
   - **Leave open** (and move it per the triage runbook) when it needs Shariant testing (the comments say whether it
     is required), a cherry-pick to a release or production branch, a manual deploy step (`manage.py
     manual_outstanding`) or install steps, or when testing failed or was inconclusive. Say which in the comment.
   - A failure outside what the issue's commits touch is a new issue to draft for the user, not a fail of this one.

## Report back

One table: issue, PR, test result, closed / left open (why), pipeline it moved to. Then what the user has to do
(decisions, human-eyes checks, Shariant tests), and the names of the data you left on the box.
