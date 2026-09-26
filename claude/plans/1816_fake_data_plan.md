# One place for fake data

Written by Claude Opus 5.5 (claude-opus-5-5), 2026-09-25 (revised to match the implementation)

Issue: https://github.com/SACGF/variantgrid/issues/1816 (Phase 2: `create_fake_data vcf | classifications | analysis`)
Status: in progress - implemented 2026-09-26, PR open; `create_fake_data all` run on vg-test2 (real annotation: uses real genes and variants, 1m41s first run, 6s re-run adding nothing)

## Why

Fake data was spread over a dozen modules with no common shape: the test annotation versions (refusing to run outside
tests), two helpers in the annotation app, three real transcripts under annotation's tests, the cohort / trio / quad /
duo / pedigree builders under snpdb's test utils, `make_cohort_genotype` in an analysis test mixin, and the two dev
subcommands (variant tags, reclassifications) behind `snpdb/management/commands/create_fake_data.py` - which imported
`analysis` and `classification` from `snpdb`, against the app dependency order.

Both dev subcommands also needed a database that already had annotated variants, so a new developer on an empty
database had no way to get a site with something on it.

Goal: one command, `manage.py create_fake_data`, that a person can run on an empty dev database and get users, labs,
genes, annotated variants, a trio with genotypes, an analysis, tags and classifications - and that tests build from the
same code. Each app owns the fake data for its own models; the command finds them without importing any app.

## Shape

### The registry (library, imports no app)

`library/fake_data.py`:

```python
class FakeData:
    name: str                    # the subcommand, e.g. "trio"
    help: str
    requires: tuple[str, ...]    # names of other FakeData this builds on - by name, never by import

@dataclass
class FakeDataContext:
    genome_build: GenomeBuild
    seed: int
    stdout: OutputWrapper
    genes: list[str]                              # filled by "genes"
    variant_ids_by_gene: dict[str, list[int]]     # filled by "variants"
    users: list[User]                             # filled by "people"
    labs: list[Lab]                               # filled by "people", germline first
    trio: Optional[Trio]                          # filled by "trio"
```

`FakeData` subclasses implement `add_arguments(parser)` (classmethod), `create(context, **options)` and
`delete(context, **options)`. `create` is idempotent: each step checks for its obviously fake rows and adds nothing when
they exist. `@register` adds a class to a module-level dict; `discover()` calls
`django.utils.module_loading.autodiscover_modules("fake_data")`, so every installed app's `<app>/fake_data.py` (or
`fake_data/__init__.py`, which imports its sibling step modules) registers itself on import. `zipf_weight` and
`in_preferred_order` (a step's preferred genes out of the context's, or all of them) live here.

Data passes between steps through the context, never through imports: `trio` in snpdb reads
`context.variant_ids_by_gene` that `variants` in annotation filled.

### The command

`snpdb/management/commands/create_fake_data.py` keeps its name and `category = "dev"`, imports `library.fake_data`
(and snpdb's own `GenomeBuild`), and builds one subparser per registered step plus `all` and `list`:

```
manage.py create_fake_data all [--genome-build GRCh37] [--seed N] [--delete] [any step's options]
manage.py create_fake_data <step> [--delete] [step options]     # runs its requirements first, on their defaults
manage.py create_fake_data list                                # steps, what each makes, requirements
```

Requirements run in dependency order before the named step, each step in its own transaction. `all` accepts every
step's options, and an option given applies to each step that has it (two steps may share a name, e.g. `--years`).
`--delete` runs the named step's delete only (for `all`, every step's delete in reverse order).

### The steps

| Step | Lives in | Requires | Creates |
|---|---|---|---|
| `people` | `snpdb/fake_data.py` | - | "Fake Health Network" org, a germline and a somatic lab, eight `fake_*` users (active, a random password printed when they are created) - the one set of labs/users every later step uses |
| `genes` | `genes/fake_data.py` | - | when the build has no real `GeneAnnotationRelease` (the test annotation's release counts as fake): the RUNX1 / GATA2 / PTEN transcripts; otherwise nothing, and later steps use `WELL_KNOWN_GENES` |
| `annotation` | `annotation/fake_data.py` | `genes` | when the build has no `VariantAnnotationVersion`: `create_fake_annotation_version`; otherwise nothing |
| `variants` | `annotation/fake_data.py` | `annotation` | on a fake annotation version: `--variants` (300) SNVs and 1bp indels at coding bases of the genes' exons (ref read from the build's fasta, whole exons at a time) with a `VariantAnnotation` each in the version's partition (gene, transcript, consequence, impact, variant class, gnomAD AF, g./c./p. HGVS), hanging off one `AnnotationRun` marked `create_fake_data variants`; on a real version: nothing created, `context.variant_ids_by_gene` read from the real annotation (variants no allele has claimed) |
| `trio` | `snpdb/fake_data.py` | `people`, `variants` | a VCF with proband / mother / father samples, cohort, `CohortGenotypeCollection` with genotype rows over up to `--max-genotypes` of the variants (inherited het from either parent, recessive, de novo, hom, missing), a `Trio` and a `Pedigree` |
| `analysis` | `analysis/fake_data/__init__.py` | `trio` | an `Analysis` over the trio: dominant and de novo TrioNodes, a population (<= 1%) and an impact (>= MODERATE) filter, a union venn - nodes saved and `update_analysis` launched |
| `tags` | `analysis/fake_data/variant_tags.py` | `people`, `variants` | the variant tags for the tag stats page, tagged by the shared users (per-user weights and active windows in `FAKE_TAGGERS`) on the context's genes |
| `classifications` | `classification/fake_data/__init__.py` | `people`, `variants` | `--alleles` (40) variants classified by one or both labs - agreeing, discordant across buckets, ~8% withdrawn - created, matched to pre-resolved allele infos and published like any record, so groupings and overlaps (discordance) are derived as usual |
| `reclassifications` | `classification/fake_data/reclassifications.py` | `people`, `variants` | the curation histories for the reclassification analytics page, curated by the shared users at the shared labs (per-lab curation behaviour in `CURATION_BY_LAB_FOCUS`) |

A version is fake when its release's `GeneAnnotationImport.url == "fake"` (what `create_fake_annotation_version`
creates); `annotation/fake_data.py:is_fake_annotation_version` says so in one place. `variants` writes annotation rows
only into a fake version, so a box with real annotation never gets fake `VariantAnnotation` rows.
`get_fake_annotation_version` keeps refusing outside tests; the `annotation` step calls `create_fake_annotation_version`
only after checking the build has no annotation version.

Alleles are shared: `snpdb/fake_data.py:fake_alleles_for_variants` makes one for each variant that has none, so tags and
classifications on a variant land on the one allele. `delete_unused_fake_alleles` deletes those no other fake data still
uses (VariantTag and Classification point at alleles with PROTECT). The classification steps share allele info creation
and deletion in `classification/fake_data/shared.py`.

`genes` and `annotation` leave their data on `--delete` (tests and everything annotated build on them); `variants`
deletes its annotation rows and keeps the `Variant` rows, which are only coordinates.

### Test builders move next to it

The builders tests use and the steps reuse live in the same per-app module, so there is one implementation of "a fake
trio":

| From | To |
|---|---|
| snpdb test utils `fake_cohort_data.py` (cohort, trio, quad, duo, pedigree) | `snpdb/fake_data.py` (names are arguments with the old values as defaults; pedigree takes an optional cohort) |
| `make_cohort_genotype` from the analysis inheritance test mixin | `snpdb/fake_data.py` |
| annotation's `tests/test_data_fake_genes.py` | `genes/fake_data.py` (`_create_fake_gene_version` / `_insert_transcript_data` became public: `create_fake_gene_version`, `insert_transcript_data`) |
| annotation's `fake_annotation.py` | `annotation/fake_data.py` |
| `get_variant_ids_by_gene` | stays in `annotation/fake_data.py` |

Every import was updated (tests, `variantgrid/test_runner.py`, `annotation/vep_annotation.py`) and the old modules
deleted. Test-only builders that nothing outside tests would use (`mme/tests/fakes.py`, `ConditionMock`,
`AnalysisSetupMixin`) stay where they are.

## Verification

- `snpdb/tests/test_create_fake_data.py`: `call_command("create_fake_data", "all", ...)` with small sizes on GRCh37;
  every step made rows; the analysis, trio, classification listing, allele, tag stats and reclassification analytics
  pages return 200; running `all` again changes no counts; `all --delete` leaves no fake rows.
- `scripts/vg imports cycles` and `lint-imports` stay clean (snpdb no longer imports analysis / classification).
- The whole suite passes with `--keepdb --parallel 4`.
- `claude/guides/testing.md` fixture index points at the new locations.
- The sparse CI fastas carry the three transcripts' exons, recorded with `FastaRecordingRunner`.
