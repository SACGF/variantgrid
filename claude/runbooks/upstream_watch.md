# Upstream watch: things we carry until upstream catches up

Each row is a fork, patched copy or workaround that exists only because an upstream fix has not shipped. The
SessionStart hook (`.claude/hooks/session_start.sh`) prints a line when any `Next check` date has passed; the
agent that sees it runs the row's check, records the date and outcome, and either pushes `Next check` out a month or
carries out `When resolved` and deletes the row. Dates are ISO (`YYYY-MM-DD`) and the hook reads the `Next check`
column (the fifth `|`-separated field), so keep the table shape.

| Item | Waiting on | Last checked | Next check |
|---|---|---|---|
| hgvs fork | biocommons/hgvs release with #899, #901, #902, #911 | 2026-10-27 | 2026-11-27 |
| UTRAnnotator copy | Ensembl/VEP_plugins#869 on the `release/<VEP>` branch, run with `all=1` | 2026-11-11 | 2026-12-11 |
| VEP plugin per-variant cache | Ensembl/ensembl-variation#1220, Ensembl/VEP_plugins#871 | 2026-11-03 | 2026-12-03 |
| CSI index docs | Ensembl/ensembl-vep#2068 | 2026-11-03 | 2027-01-27 |
| Pre-112 VEP fasta | Shariant and runx1db2 moving to VEP >= 112 (Ensembl/ensembl-vep#1635 fixed in 112) | 2026-11-03 | 2027-02-03 |
| MaveDB stripped file | every deployment running MaveDB on a VEP_plugins `release/116` checkout from 2026-09-01 (Ensembl/VEP_plugins#853) | 2026-11-03 | 2027-02-03 |
| gnomAD SV overlap via two `--custom` | Ensembl/ensembl-vep#1989 (`coords=2`, `show_overlaps`) for #1645 | 2026-11-03 | 2027-02-03 |
| AnnotSV regulatory element sort | lgmgeo/AnnotSV#356 fixed in a released Annotations_Human bundle | 2026-11-10 | 2026-12-10 |
| chrX PAR table | a bioutils release with `bioutils.par` that installs on Python 3.14 (0.6.2 caps `<3.14`, biocommons/bioutils#97) | 2026-11-11 | 2026-12-11 |
| HPO phenotype synonyms | an HPO release after 2026-10-24 carrying obophenotype/human-phenotype-ontology#11777-#11799 | 2026-11-12 | 2026-12-12 |

## chrX PAR table

`snpdb/models/models_genome.py:CHRX_PSEUDOAUTOSOMAL_REGIONS` hardcodes chrX's pseudoautosomal regions for the sex guess
(#2745). `bioutils.par` (biocommons/bioutils#88) has the same regions for GRCh37, GRCh38 and T2T-CHM13v2.0, keyed by
our `GenomeBuild` names, but 0.6.2 declares `requires-python <3.14` (biocommons/bioutils#97; lifted on main by
biocommons/bioutils#98, so far only in the 0.6.3a0 prerelease).

Check:

```bash
gh release list --repo biocommons/bioutils --limit 3
curl -s https://pypi.org/pypi/bioutils/json | python3 -c "import json,sys; i=json.load(sys.stdin)['info']; print(i['version'], i['requires_python'])"
```

When resolved (a final release with `bioutils.par` and no `<3.14` cap): bump `bioutils` (it is transitive, so add it to
`requirements.in`) and recompile both requirements files (root `AGENTS.md`, Python packages). Build
`annotation/tasks/calculate_sample_stats.py:_pseudoautosomal_expression`'s ranges from
`bioutils.par.get_par_map(genome_build.name)["X"].values()` as `(start + 1, end)` (bioutils is interbase, the
`position__range` lookup 1-based inclusive), returning `Value(False)` on the `ValueError` for a build with no PAR data.
Delete `CHRX_PSEUDOAUTOSOMAL_REGIONS`, point `test_pseudoautosomal_regions_stay_out_of_the_chrx_sex_guess` at
`get_par_map`, and delete this row.

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

`annotation/annotation_data/generate_annotation/UTRAnnotator.pm` is Ensembl's plugin with the first revision of
Ensembl/VEP_plugins#869 (8891949) applied: every 5' UTR effect reported, `&`-joined in a fixed order (stock output
differs between runs, Ensembl/VEP_plugins#870). `download_annotation_data` installs it over Ensembl's because
UTRAnnotator is in `_VARIANTGRID_VEP_PLUGINS` (`annotation/vep_annotation.py:is_variantgrid_vep_plugin`).

The PR as revised for review (c27df0e4) leaves the default output one effect per transcript, now the first in the
order `uAUG_gained uSTOP_lost uAUG_lost uSTOP_gained uFrameshift`, and reports every effect with the plugin option
`all=1`. With `all=1` its output is our copy's plus a `consequence=<SO term>` field on each uORF entry. We want
`all=1`: the high-impact 5' UTR filter (`annotation/effect_filters.py:utr_high_impact_rules`) matches each effect's
own records, so a variant whose high-impact effect is not the first would drop out of it under the default. Our copy
ignores `all`, so the option can go on the `--plugin` string before the switch.

The plugin is installed per deployment, not per version, and a VEP upgrade brings a new VariantAnnotationVersion, so
the switch happens with the VEP upgrade to a release carrying #869. A pending `backfilled_utr_effects` pass needs
our copy installed (`annotation/vep_annotation.py:installed_vep_plugin_is_ours`), so it runs before that upgrade.

Check, with `<VEP>` the `ANNOTATION_VEP_VERSION` a new version would install:

```bash
gh pr view 869 --repo Ensembl/VEP_plugins --json state,mergedAt,commits \
    -q '"\(.state) \(.mergedAt) \(.commits[-1].oid[0:8])"'
curl -sL https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/<VEP>/UTRAnnotator.pm \
    | grep -c 'param_hash->{all}'
```

A new commit on the PR: diff it against c27df0e4 and update this section if the `all=1` output changes.

When resolved (#869 on `release/<VEP>`, every deployment upgraded to that VEP and none with `backfilled_utr_effects`
pending):
- add `,all=1` to the `VEPPlugin.UTR_ANNOTATOR` `--plugin` string (`annotation/vep_annotation.py:get_vep_command`),
  delete our `.pm` and take UTRAnnotator out of `_VARIANTGRID_VEP_PLUGINS`.
- drop the `installed_vep_plugin_is_ours` guard on `--mark-backfilled backfilled_utr_effects`
  (`annotation/management/commands/annotation_backfill_columns.py`).
- the `consequence` field: `annotation/effect_filters.py` still works out each record's effect from its keys, which
  older versions need, and its lookaheads skip the extra key. Put `consequence` first in
  `annotation/vep_plugin_display.py:_UTR_ANNOTATOR_KEY_ORDER`, where it heads each record on the variant page.
- update the UTRAnnotator paragraphs in `claude/research/annotation.md` (effects backfill, data files step, transcript
  rows) and delete this row.

If #869 is closed unmerged, new versions keep our copy: note why here.

## VEP plugin per-variant cache

Ensembl/ensembl-vep#2069: a base-class helper (ensembl-variation#1220) and six plugins converted to it
(VEP_plugins#871) cut VEP CPU by a quarter on coding-heavy input. We ship nothing for this; the patched files are the
upstream PRs' diffs and the measurements are in the comments on #2499.

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

## Pre-112 VEP fasta

`variantgrid/settings/components/annotation_settings.py:use_pre_vep112_fasta` points VEP at the chromosome-named
fasta, because VEP before 112 renamed NCBI-fasta contigs and silently dropped plugin/custom annotations
(Ensembl/ensembl-vep#1635). Shariant (`variantgrid/settings/env/shariantcommon.py`) and runx1db2
(`variantgrid/settings/env/runx1db2.py`) pin `ANNOTATION_VEP_VERSION = "110"` and call it.

Check: `grep -rn "use_pre_vep112_fasta\|ANNOTATION_VEP_VERSION *=" variantgrid/settings/env/`. When resolved (no env
file calls it): delete `use_pre_vep112_fasta`, the `vep_config["fasta"]` fallback in
`annotation/vep_annotation.py:_get_vep_fasta` stays only while a `vep_config` sets `fasta`, and delete this row.
`annotation/annotation_run_files.py:write_qs_to_vcf` writing contig names rather than accessions is independent of
this and stays.

## MaveDB stripped file

`annotation/annotation_data/generate_annotation/mavedb_strip.py` cuts the MaveDB download to the columns MaveDB.pm
emits, because the plugin hashed all ~250 columns of every record in a variant's span and a single long deletion
could OOM VEP. Ensembl/VEP_plugins#853 (on `release/116` from 2026-09-01) makes the plugin keep only those columns
and skip VRS-mapped records, with output unchanged, so the strip is redundant on that plugin - but a VEP 110
deployment, or a 116 plugins checkout made before then, still has the old MaveDB.pm.

Check: the pre-112 row above, plus `git -C <ANNOTATION_VEP_PLUGINS_DIR> log -1 --format=%cs -- MaveDB.pm` on each 116
deployment. When resolved: stop generating `.stripped` files (keep accepting the name in
`annotation/vep_annotation.py`, deployments already have them), delete `mavedb_strip.py` and this row.

## gnomAD SV overlap via two `--custom`

`annotation/vcf_files/bulk_vep_vcf_annotation_inserter.py:SVOverlapProcessor` runs the gnomAD SV file as two
`--custom` entries and stitches them back together, because `--custom` cannot return both the matched record's ID and
its coordinates (Ensembl/ensembl-vep#1645, open; the feature is Ensembl/ensembl-vep#1989). The StructuralVariantOverlap INV bug
(Ensembl/VEP_plugins#710) is fixed from VEP 114.

Check: `gh pr view 1989 --repo Ensembl/ensembl-vep --json state,mergedAt -q '"\(.state) \(.mergedAt)"'`. When resolved (in a release
deployments run): replace the two customs with one and drop the stitching in `SVOverlapProcessor`, then delete this
row.

## AnnotSV regulatory element sort

The `Annotations_Human` 3.5 bundle ships unsorted `*.sorted.bed` regulatory element files: ABC and EnhancerAtlas for GRCh38 and
CHM13, miRTargetLink for CHM13 and (ENSEMBL only) GRCh37. AnnotSV's `bedtools intersect -sorted` fails on them and it falls back, without a message,
to an intersect that loads the file: 1.4 GB per run instead of a few MB (#2749). Until a bundle ships them sorted,
each deployment sorts its own copy (`annotation/annotsv_annotation.py:get_annotsv_sort_regulatory_elements_script`: printed by a manual migration
step for existing deployments and by a `deployment_check` warning for new servers and bundles), per file, as the
files' owner:

```bash
LC_ALL=C sort -k1,1 -k2,2n <file> > <file>.tmp && mv <file>.tmp <file>
```

over `<ANNOTATION_ANNOTSV_ANNOTATIONS_DIR>/Annotations_Human/FtIncludedInSV/RegulatoryElements/*/*.sorted.bed`,
and only those: re-sorting the bundle's other bed files changes AnnotSV's output.
The output TSV is byte-identical either way.

Check:

```bash
gh issue view 356 --repo lgmgeo/AnnotSV --json state,closedAt -q '"\(.state) \(.closedAt)"'
```

and, when a newer bundle than 3.5 is out, on a fresh copy of it (prints nothing when sorted):

```bash
for f in Annotations_Human/FtIncludedInSV/RegulatoryElements/*/*.sorted.bed; do LC_ALL=C sort -c -s -k1,1 -k2,2n "$f"; done
```

When resolved (a released bundle passes that loop): take the sort out of the AnnotSV install paragraph in
`claude/guides/operations.md` and delete this row. The `deployment_check` warning stays - it is cheap and catches a
bundle that regresses.

## HPO phenotype synonyms

`patients/phenotype_matcher.py:PhenotypeMatcher._get_special_case_lookups` hardcodes free-text phrases that no HPO
term's name or synonym covers. obophenotype/human-phenotype-ontology#11778-#11799 added synonyms for 22 of them
(closed 2026-10-24, after the v2026-09-01 release), e.g. "Elevated renin", "Hot flushes", "Crowded dentition",
"Fetal hydrops". #11777 ("LQTS" on HP:0001657) was declined as a disease name on a finding; #11806, asking that the
term's existing "Long QT syndrome" synonyms go too, is open, so keep the `"LQTS"` entry.

Check:

```bash
gh release list --repo obophenotype/human-phenotype-ontology --limit 3
gh issue view 11806 --repo obophenotype/human-phenotype-ontology --json state -q .state
```

When resolved (a release dated after 2026-10-24): import its `hp.owl` (`manage.py ontology_import --hpo <file>`) on
a test deployment, then for each special-case key check whether the matcher's alias index (built from
`OntologyTerm.aliases` in `PhenotypeMatcher.__init__`) now resolves it to the same term. Delete those entries
(several keys differ from the synonym that was added - "increased renin" vs "Elevated renin", "recurrent urtis" vs
"Recurrent URTI" - so only the exact matches go), run `patients.tests.test_phenotype_matching`, and delete this row.
If #11806 is still open, carry it as its own row.
