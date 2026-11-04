# Upstream watch: things we carry until upstream catches up

Each row is a fork, patched copy or workaround that exists only because an upstream fix has not shipped. The
SessionStart hook (`.claude/hooks/session_start.sh`) prints a line when any `Next check` date has passed; the
agent that sees it runs the row's check, records the date and outcome, and either pushes `Next check` out a month or
carries out `When resolved` and deletes the row. Dates are ISO (`YYYY-MM-DD`) and the hook reads the `Next check`
column (the fifth `|`-separated field), so keep the table shape.

| Item | Waiting on | Last checked | Next check |
|---|---|---|---|
| hgvs fork | biocommons/hgvs release with #899, #901, #902, #911 | 2026-10-27 | 2026-11-27 |
| UTRAnnotator copy | Ensembl/VEP_plugins#869 on the `release/<VEP>` branch | 2026-11-03 | 2026-12-03 |
| VEP plugin per-variant cache | Ensembl/ensembl-variation#1220, Ensembl/VEP_plugins#871 | 2026-11-03 | 2026-12-03 |
| CSI index docs | Ensembl/ensembl-vep#2068 | 2026-11-03 | 2027-01-27 |
| Pre-112 VEP fasta | Shariant and runx1db2 moving to VEP >= 112 (Ensembl/ensembl-vep#1635 fixed in 112) | 2026-11-03 | 2027-02-03 |
| MaveDB stripped file | every deployment running MaveDB on a VEP_plugins `release/116` checkout from 2026-09-01 (Ensembl/VEP_plugins#853) | 2026-11-03 | 2027-02-03 |
| gnomAD SV overlap via two `--custom` | Ensembl/ensembl-vep#1989 (`coords=2`, `show_overlaps`) for #1645 | 2026-11-03 | 2027-02-03 |
| cdot thread-local FASTA swap | SACGF/cdot#234 in a cdot release | 2026-11-04 | 2026-12-04 |
| uv cooldown exceptions | the 7-day cooldown passing the pinned cutoffs in `uv.toml` | 2026-11-04 | 2026-11-06 |

## uv cooldown exceptions

`uv.toml` names packages whose security fix was less than seven days old when it was needed, each with a fixed
cutoff date that admits it. The entry freezes that package at the cutoff, so it comes out once the default cooldown
has passed the date: pypdf after 2026-11-05.

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
its coordinates (Ensembl/ensembl-vep#1645, open; the feature is Ensembl/ensembl-vep#1989). The StructuralVariantOverlap INV bug that first ruled the plugin out
(Ensembl/VEP_plugins#710) is fixed from VEP 114.

Check: `gh pr view 1989 --repo Ensembl/ensembl-vep --json state,mergedAt -q '"\(.state) \(.mergedAt)"'`. When resolved (in a release
deployments run): replace the two customs with one and drop the stitching in `SVOverlapProcessor`, then delete this
row.

## cdot thread-local FASTA swap

cdot's `GenomeFastaSeqFetcher` opens one pysam `FastaFile` per FASTA, and a handle read from two threads or from both
sides of a fork breaks for the life of the process (#2624). `HGVSMatcher.instance` and `get_genome_babelfish` cache
their seqfetchers per process, so `library/genomics/fasta_wrapper.py:make_contig_fastas_thread_local` swaps the values
of their `contig_fastas` for per-thread, per-pid `ThreadLocalFastaFile` proxies
(`genes/hgvs/biocommons_hgvs/data_provider.py:SingleBuildFastaSeqFetcher`, `snpdb/babelfish.py:GenomeFastaSequenceProvider`).
It leans on cdot only using those values through `.fetch()`, `.filename` and `in`.

Check: `gh issue view 234 --repo SACGF/cdot --json state -q .state` and the cdot release notes
(`gh release list --repo SACGF/cdot --limit 3`).

When resolved (a cdot release opens handles per thread and per process): bump cdot in `requirements.in`, recompile
both lock files, delete the two `make_contig_fastas_thread_local` calls, `ThreadLocalFastaFile`,
`make_contig_fastas_thread_local` and `library/tests/test_fasta_wrapper.py`'s tests of them, keep a threaded
regression test against the matcher, and delete this row.
