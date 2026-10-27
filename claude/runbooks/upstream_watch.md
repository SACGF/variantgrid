# Upstream watch: things we carry until upstream catches up

Each row is a fork, patched copy or workaround that exists only because an upstream fix has not shipped. The
SessionStart hook (`.claude/hooks/session_start.sh`) prints a line when any `Next check` date has passed; the
agent that sees it runs the row's check, records the date and outcome, and either pushes `Next check` out a month or
carries out `When resolved` and deletes the row. Dates are ISO (`YYYY-MM-DD`) and the hook reads the fourth and fifth
columns, so keep the table shape.

| Item | Waiting on | Last checked | Next check |
|---|---|---|---|
| hgvs fork | biocommons/hgvs release with #899, #901, #902, #911 | 2026-10-27 | 2026-11-27 |
| UTRAnnotator copy | Ensembl/VEP_plugins#869 on the `release/<VEP>` branch | 2026-10-27 | 2026-11-27 |
| VEP plugin per-variant cache | Ensembl/ensembl-variation#1220, Ensembl/VEP_plugins#871 | 2026-10-27 | 2026-11-27 |
| CSI index docs | Ensembl/ensembl-vep#2068 | 2026-10-27 | 2027-01-27 |
| uv cooldown exceptions | the 7-day cooldown passing the pinned cutoffs | 2026-10-27 | 2026-10-30 |

## uv cooldown exceptions

`uv.toml` names packages whose security fix was less than seven days old when it was needed, each with a fixed
cutoff date that admits it. The entry freezes that package at the cutoff, so it comes out once the default cooldown
has passed the date: certifi and pillow after 2026-10-28, pypdf after 2026-10-29.

Check: `grep exclude-newer-package uv.toml` against today's date. When resolved: delete the expired entries, recompile
both requirements files (root `AGENTS.md`, Python packages), and delete this row when none are left.

## hgvs fork

`requirements.in` pins `hgvs` to a commit on davmlaw/hgvs: biocommons main plus four PRs that are merged upstream but
not in a release, plus a dup/inv guard near the start of a contig that has no upstream PR yet (#2103).

Check:

```bash
gh release list --repo biocommons/hgvs --limit 3
for n in 899 901 902 911; do gh pr view $n --repo biocommons/hgvs --json state,mergedAt -q '"\(.state) \(.mergedAt)"'; done
```

When resolved (a release tag newer than the merge dates): the guard still has to go upstream or stay as a fork on top of
the release. Either way, repoint `requirements.in` at the release or the new fork commit, recompile both requirements
files (root `AGENTS.md`, Python packages) and run `genes.tests.test_hgvs` and `genes.tests.test_hgvs_corpus`.

## UTRAnnotator copy

`annotation/annotation_data/generate_annotation/UTRAnnotator.pm` is Ensembl's plugin with Ensembl/VEP_plugins#869
applied (every 5' UTR effect reported in a fixed order; stock output differs between runs, Ensembl/VEP_plugins#870).
`download_annotation_data` installs it over Ensembl's because UTRAnnotator is in `_VARIANTGRID_VEP_PLUGINS`
(`annotation/vep_annotation.py:is_variantgrid_vep_plugin`).

Check, with `<VEP>` the `ANNOTATION_VEP_VERSION` deployments install:

```bash
gh pr view 869 --repo Ensembl/VEP_plugins --json state,mergedAt -q '"\(.state) \(.mergedAt)"'
diff <(curl -sL https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/<VEP>/UTRAnnotator.pm) \
     <(tail -n +10 annotation/annotation_data/generate_annotation/UTRAnnotator.pm) && echo "release branch has it"
```

When resolved (the diff is empty, or differs only in a way we do not need): delete the file, take UTRAnnotator out of
`_VARIANTGRID_VEP_PLUGINS`, drop its paragraph from `claude/research/annotation.md` (data files step) and this row. If
the PR is closed unmerged, keep the copy and note why here.

## VEP plugin per-variant cache

Ensembl/ensembl-vep#2069: a base-class helper (ensembl-variation#1220) and six plugins converted to it
(VEP_plugins#871) cut VEP CPU by a quarter on coding-heavy input. We ship nothing for this; the patched files and the
benchmark harness are in the scratch dir named in the session memory for #2499.

Check:

```bash
gh pr view 1220 --repo Ensembl/ensembl-variation --json state,mergedAt -q '"\(.state) \(.mergedAt)"'
gh pr view 871 --repo Ensembl/VEP_plugins --json state,mergedAt -q '"\(.state) \(.mergedAt)"'
```

When resolved (merged): nothing to do until deployments move to the VEP release that carries it; note the version on
#2499 and delete this row. If closed unmerged: decide on #2499 whether to carry the six plugins the way UTRAnnotator is
carried (the base-class half lives in ensembl-vep's vendored `Bio/` tree, which the plugins dir cannot override).

## CSI index docs

Ensembl/ensembl-vep#2068 suggests the VEP docs mention `tabix -C -m 12` for custom/plugin files. Informational only:
our indexes are built and distributed regardless (`annotation/annotation_data/generate_annotation/build_tabix_csi.sh`,
`download_annotation_data`).

Check: `gh issue view 2068 --repo Ensembl/ensembl-vep --json state -q .state`. When closed, delete this row.
