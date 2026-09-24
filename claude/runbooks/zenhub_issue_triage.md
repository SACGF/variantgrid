# ZenHub issue queries and release triage

Issues for VariantGrid are spread over several GitHub repos (`SACGF/variantgrid`, `SACGF/variantgrid_private`,
`SACGF/variantgrid_shariant`, `SACGF/variantgrid_sapath`) but tracked on one ZenHub board. GitHub has the issue text and
comments (`gh issue view`); ZenHub has the pipeline (board column) each issue is in, which GitHub cannot see.

## 1. Get a token

1. Go to https://app.zenhub.com/settings/tokens and create a **GraphQL Personal API Key** (starts with `zh_`).
2. Put it in your shell profile, e.g. `export ZENHUB_TOKEN=zh_...` in `~/.bashrc`, so agents see it as `$ZENHUB_TOKEN`.
   Keep it out of the repo and out of any command output.

## 2. The API

One GraphQL endpoint: `POST https://api.zenhub.com/public/graphql` with `Authorization: Bearer $ZENHUB_TOKEN`.
The interactive explorer is at https://developers.zenhub.com/explore.

Issues are looked up by the GitHub repo's numeric id plus the issue number. Ids (from `gh api repos/SACGF/<repo> --jq .id`):

| Repo | GitHub id |
|---|---|
| `SACGF/variantgrid` | `299486514` |

Look the others up the same way when needed.

Which pipeline is an issue in:

```bash
curl -s https://api.zenhub.com/public/graphql \
  -H "Authorization: Bearer $ZENHUB_TOKEN" -H "Content-Type: application/json" \
  -d '{"query":"query { issueByInfo(repositoryGhId: 299486514, issueNumber: 1477) { number title state pipelineIssues { nodes { pipeline { name } workspace { name } } } } }"}'
```

Every issue in a pipeline (all repos), paged 100 at a time:

```graphql
query($id: ID!, $after: String) {
  searchIssuesByPipeline(pipelineId: $id, filters: {}, first: 100, after: $after) {
    totalCount
    pageInfo { hasNextPage endCursor }
    nodes { number title state repository { name ownerName } labels { nodes { name } } }
  }
}
```

The pipelines and their ids come from the workspace, e.g. via
`issueByInfo(...) { pipelineIssues { nodes { workspace { id name pipelinesConnection(first: 50) { nodes { id name } } } } } }`.

## 3. The board

Workspace **Everything Space** (`5bb3158c4b5806bc2beae448`). Pipelines, in board order:

| Pipeline | Id |
|---|---|
| New Issues | `Z2lkOi8vcmFwdG9yL1BpcGVsaW5lLzE1MTIzMzE` |
| Icebox | `Z2lkOi8vcmFwdG9yL1BpcGVsaW5lLzE1MTI0MTM` |
| Backlog | `Z2lkOi8vcmFwdG9yL1BpcGVsaW5lLzE1MTI1NDU` |
| Work on Next | `Z2lkOi8vcmFwdG9yL1BpcGVsaW5lLzIxNzQ4OTE` |
| In Progress | `Z2lkOi8vcmFwdG9yL1BpcGVsaW5lLzE1MTI2NDI` |
| Coding Complete | `Z2lkOi8vcmFwdG9yL1BpcGVsaW5lLzE1MTI3MTM` |
| Review/QA (Shariant) | `Z2lkOi8vcmFwdG9yL1BpcGVsaW5lLzE1MTI3Nzk` |
| Post QA Discussion | `Z2lkOi8vcmFwdG9yL1BpcGVsaW5lLzE1MTI4MzU` |
| VariantGrid test | `Z2lkOi8vcmFwdG9yL1BpcGVsaW5lLzIyMzAxNTk` |
| SA Path test | `Z2lkOi8vcmFwdG9yL1BpcGVsaW5lLzIzMzEyNDU` |
| Pending Shariant Release | `Z2lkOi8vcmFwdG9yL1BpcGVsaW5lLzE1MTI4ODY` |
| Pending SA Path release | `Z2lkOi8vcmFwdG9yL1BpcGVsaW5lLzI1MzY2OTI` |
| Pending variantgrid.com release | `Z2lkOi8vcmFwdG9yL1BpcGVsaW5lLzI5ODQyNjM` |

## 4. After a release: which issues went out

The usual job is: for each issue in Coding Complete, Review/QA (Shariant) and Post QA Discussion, did it make the release
branch (e.g. `origin/shariant_prod_2025_04`), and if so can it be closed.

1. `git fetch origin`, then list the pipelines' issues with the query above.
2. Match commits to issues by their message. References come in several shapes: `SACGF/variantgrid_private#3761`,
   `variantgrid#1222`, `issue #1206`, `SACGF/variantgrid##1136` (double hash) and a title glued to the number
   (`badly#1312`). A plain `#N` means `SACGF/variantgrid`, or `variantgrid_private` when N is above about 2000.
   Shariant and SA Path issue numbers clash with old public ones (`#187` in 2021 is a public issue), so for those,
   check the code rather than trusting a plain ref.
3. A release branch is master at a branch point plus cherry-picks, so test each commit with
   `git merge-base --is-ancestor <sha> <branch>` and also match on subject (a cherry-pick has a new sha).
   An issue with some commits in the branch and later ones only on master is only partly released.
4. For issues with no commit refs, check the change itself against the branch: `git grep` for the fix in the branch versus
   master, or compare library versions in `uicore/templates/uicore/page/base.html`. Issues reported against Shariant
   test after the branch's last commit are for the next release.
5. Read each candidate's comments (`gh issue view N -R SACGF/<repo> --json body,comments`). Close-ready means the code is in
   the branch *and* testing passed or wasn't needed. Leave for a person: an open test request, a pending decision, a
   `ManualOperation` / data fix that has to run on prod, or process-only issues (data imports, mappings, ClinVar clean-up).
6. For close-ready issues: comment (prefixed `🤖 Written by Claude`) saying it went into the named prod branch and has
   been released, with the commits or the test pass as evidence, then `gh issue edit N -R SACGF/<repo> --add-label
   "can we close this?"`. Humans close issues; agents never do.

## 5. Triage Coding Complete into test pipelines

Move an issue with the `id` field from `searchIssuesByPipeline` (a ZenHub id, not the GitHub number):

```graphql
mutation($i: MoveIssueInput!) { moveIssue(input: $i) { issue { number } } }
# variables: {"i": {"issueId": "<zenhub id>", "pipelineId": "<target pipeline id>", "position": 0}}
```

For each issue, find its commits on `origin/master` (section 4, step 2) and send it to one of:

- **Review/QA (Shariant)**: a user-visible change in anything Shariant has: classifications, allele/variant pages, search,
  gene pages, conditions, discordance, ClinVar export, lab/user admin, Keycloak, and all shared UI/JS/CSS/grid code.
  Include checks that a feature turned off for Shariant (`shariantcommon.py`) stays off.
- **VariantGrid test**: only in features Shariant lacks (analysis, VCF/sample upload, samples, patients, seqauto, the
  annotation upgrade, since Shariant pins `columns_version` 3). Comment `Shariant testing not required: <reason>`.
- **Stay in Coding Complete**, with a comment on why no manual test is needed: dead code, a refactor with no behaviour
  change, docs, deploy tooling, or a rare-crash fix covered by a unit test. List these for the user to review.
- **Stay, no comment, raise with the user**: no commits on master. The work may be on an unmerged branch
  (`git branch -r --contains <sha>`), or the issue may still have an open decision.
