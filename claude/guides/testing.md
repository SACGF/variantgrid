# Testing guide

Verified against fbff51279 on 2026-09-10

The fixture index for "how do I conjure a Variant / Sample / Trio / Classification / Analysis in three
lines of a test" (`claude/plans/agent_system.md` §4.3). Every builder line is `path:function(signature)`;
read the builder before relying on a detail not stated here.

## Running tests

```bash
python3 manage.py test --keepdb                                             # whole suite (~3 min)
python3 manage.py test --keepdb snpdb.tests.test_variant                    # one module
python3 manage.py test --keepdb snpdb.tests.test_variant.VariantTest.test_x # one class / method
python3 manage.py vg tests --changed                                        # only the modules a change puts at risk
```

- `--keepdb` skips re-running hundreds of migrations; always use it. Drop it only when a migration you wrote changes.
- `vg tests --changed` (contract in `agent_system.md` §4.2) maps changed files to test labels through the
  first-party import graph and prints `manage.py test` args: a changed module selects every test module that
  transitively imports it, templates / JS / `urls.py` select the app's `tests.test_urls`, a migration selects the
  app's whole `tests` package. The engine is `library/vg/test_selection.py:select_tests(changed=None, base=None)`.
  `scripts/vg tests --explain` shows which changed file selected each label without booting Django; `--run`
  executes them.
- `--parallel 4 --keepdb` runs the whole suite (3,212 tests) in about 65 s wall on vg-test2 (4 cores) with
  nothing failing (2026-09-10); serially it is about 3 minutes. Django clones `test_snpdb` into
  `test_snpdb_1..4`. `VariantGridTestRunner` drops the clones before each `--keepdb` run (about 3 s each to
  recreate) because Django would otherwise reuse them unmigrated. CI runs it the same way (`.github/workflows/django-tests.yml`). A test that fails only in parallel
  is sharing a file under `data/` or a fixed temp path - give it its own directory rather than marking it serial.
  `tblib` (in requirements) lets a worker ship a failing test's traceback back to the parent; without it any
  error under `--parallel` aborts the whole run with `cannot pickle 'traceback' object`.
- `VariantGridTestRunner.setup_databases` seeds the fake annotation versions for GRCh37 and GRCh38 into the
  main test database *before* Django clones it, so the per-class `get_fake_annotation_version` calls are a few
  `get_or_create` lookups rather than 240 ontology re-imports. Two consequences: a test that counts
  `AnnotationVersion` / `VariantAnnotationVersion` / `OntologyTerm` rows sees the seeded ones, and a test that
  needs its own ACTIVE `VariantAnnotationVersion` (`one_active_vav_per_build`) calls
  `annotation/fake_annotation.py:retire_seeded_annotation_version(genome_build)` first. Under `--keepdb` the
  seeded rows persist between runs; the seed is idempotent.
- Every app package needs an `__init__.py`: `manage.py test <app>.tests` fails at discovery with
  `expected str ... not NoneType` when the app is an implicit namespace package.
- `TEST_RUNNER` is `variantgrid/test_runner.py:VariantGridTestRunner` (see External services below).

## Fixture builders

Almost everything wants a `GenomeBuild` and a `User`; get them with
`GenomeBuild.get_name_or_alias("GRCh37")` (or `GenomeBuild.grch37()` / `.grch38()`) and
`User.objects.get_or_create(username=...)[0]`. Builders are `get_or_create` based and safe to call in `setUpTestData`.

### Users, labs, organisations
- `classification/tests/models/test_utils.py:ClassificationTestUtils.setUp()` → creates org `instx`, labs `instx/labby` and external `instx/ext`, users `joejoe` (labby) and `joejoe2` (both); read back with `.lab_and_user()` / `.external_lab_and_user()` (each returns `(Lab, User)`).
- No shared Lab/Org builder outside that: the pattern is `Organization.objects.get_or_create(name=, group_name=)`, `Country.objects.get_or_create(name=)`, `Lab.objects.get_or_create(name=, city=, country=, organization=, group_name="org/lab")` then `lab.group.user_set.add(user)` (copy from `classification/tests/utils/test_urls.py` `setUpTestData`).
- Object permissions: `library/guardian_utils.py:assign_permission_to_user_and_groups(user, obj)`; `admin_bot()` there is uncached under `UNIT_TEST` on purpose.

### Variants and alleles
- `snpdb/tests/utils/vcf_testing_utils.py:slowly_create_test_variant(chrom, position, ref, alt, genome_build)` → `Variant` (creates Sequence/Locus, normalises symbolic alts); needs the build's contigs (always loaded). Bypasses `VariantPKLookup`, test only.
- `snpdb/tests/utils/vcf_testing_utils.py:create_mock_allele(variant, genome_build)` → `Allele` linked via a `VariantAllele` (origin IMPORTED_TO_DATABASE, tool DBSNP).
- `snpdb/tests/utils/vcf_testing_utils.py:slowly_create_loci_and_variants_for_vcf(genome_build, vcf_filename, get_variant_id_from_info=False)` → None; inserts every record of a VCF as Loci/Variants (first ALT only, honours SVLEN); `get_variant_id_from_info` pins pks from `INFO/variant_id`.
- `annotation/fake_annotation.py:create_fake_variants(genome_build)` → None; loads `annotation/tests/test_data/test_columns_version1_<build>.vep_annotated.vcf` (4 records) via the function above, pks from INFO. Pick one with `Variant.objects.filter(Variant.get_no_reference_q()).first()`.
- `genes/tests/gene_fusion_test_utils.py:create_gene_fusion(gene_a, gene_b, directionality_known=True, resolver=None)` → `GeneFusion` (plus its gene-level `Variant`), built the way the VCF insert pipeline would; needs the gene symbols to resolve.

### VCF, samples, cohorts, trios
- `snpdb/tests/utils/fake_cohort_data.py:create_fake_cohort(user, genome_build)` → `Cohort` of 3 samples (`proband`, `mother`, `father`) on one `VCF` (`genotype_samples=1`, import SUCCESS, a `VCFFilter`), with a `CohortGenotypeCollection`; permissions assigned to `user`.
- `snpdb/tests/utils/fake_cohort_data.py:create_fake_trio(user, genome_build)` → `Trio` (mother affected, father not) over `create_fake_cohort`.
- `snpdb/tests/utils/fake_cohort_data.py:create_fake_quad(user, genome_build, sibling_affected=False)` → `Quad` over its own 4-sample cohort (`proband`, `mother`, `father`, `sibling`).
- `snpdb/tests/utils/fake_cohort_data.py:create_fake_duo(user, genome_build, relationship=DuoRelationship.MOTHER, parent_affected=False)` → `Duo` over its own 2-sample cohort (`proband`, `parent`).
- `snpdb/tests/utils/fake_cohort_data.py:create_fake_pedigree(user, genome_build)` → `Pedigree` (`PedFile` + `PedFileFamily`) over `create_fake_cohort`.
- `analysis/tests/inheritance_node_mixin.py:make_cohort_genotype(cgc, variant, samples_zygosity)` → None; writes one `CohortGenotype` row for `cgc` (`trio.cohort.cohort_genotype_collection`). Zygosity string per sample in cohort order: `E`=HET, `R`=HOM_REF, `O`=HOM_ALT, `U`=UNKNOWN, `.`=MISSING, so a trio `"ERR"` is a het de novo.
- Single-sample VCF + `CohortGenotypeCollection`: copy `_create_vcf_sample` from `analysis/tests/test_sample_node_levels.py` (private, but the canonical 12 lines).

### Genes, transcripts, annotation versions
- `annotation/tests/test_data_fake_genes.py:create_fake_transcript_version(genome_build, release=None)` → `TranscriptVersion` ENST00000300305.7 (RUNX1, Ensembl, minus strand, chr21) with real exon data for GRCh37 and GRCh38; passing a `GeneAnnotationRelease` also links the symbol into that release. `tv.gene_version.gene_symbol` is the usual next step.
- `annotation/tests/test_data_fake_genes.py:create_gata2_transcript_version(genome_build)` → NM_001145661.2 (GATA2, RefSeq, minus strand, chr3).
- `annotation/tests/test_data_fake_genes.py:create_pten_transcript_version(genome_build)` → NM_000314.8 (PTEN, RefSeq, plus strand, chr10) - use with GATA2 to cover both orientations.
- `mme/tests/fakes.py:make_gene_version(gene_id, symbol, annotation_consortium, version=1)` → `GeneVersion` on GRCh38; call twice with the same `gene_id` to model a symbol rename.
- `annotation/fake_annotation.py:get_fake_annotation_version(genome_build)` → `AnnotationVersion` (the thing `Analysis.set_defaults_and_save` and most views need). Creates a `GeneAnnotationRelease` "42.20240101", the test ontology, an ACTIVE `VariantAnnotationVersion` (columns_version 2, gnomAD 2.1.1/3.1, dbNSFP 4.0a), a `ClinVarVersion` and an HPA version. Raises unless `settings.UNIT_TEST`. Call it once per build in `setUpTestData`.
- `annotation/fake_annotation.py:retire_seeded_annotation_version(genome_build)` → None; demotes the runner-seeded ACTIVE `VariantAnnotationVersion` to HISTORICAL, for a test that creates its own ACTIVE one (only one per build is allowed).
- `annotation/fake_annotation.py:get_fake_vep_version(genome_build, annotation_consortium, columns_version)` → dict of `VariantAnnotationVersion` kwargs (not saved); for tests that need a VAV at another `columns_version` - pair it with `FIXTURE_VEP_VERSIONS` and the `test_columns_version<N>_<build>.vep_annotated.vcf` fixtures in `annotation/tests/test_data/`.
- `annotation/fake_annotation.py:create_fake_variant_annotation(variant, variant_annotation_version)` → `VariantAnnotation` (with its `AnnotationRangeLock` and `AnnotationRun`); only `hgvs_g` is filled, set other columns on the returned row.
- `annotation/fake_annotation.py:create_fake_clinvar_data(clinvar_version)` → None; a `ClinVar` row (Pathogenic, 42/42 ids) on the first fake variant plus a PubMed `ClinVarCitation`; calls `create_fake_variants` itself.
- `ontology/tests/test_data_ontology.py:create_ontology_test_data()` → None; imports `ontology/tests/test_data/biomart_omim.tsv` and `small.owl` (OMIM + a small HPO). `create_test_ontology_version()` → `OntologyVersion`. Both already run inside `get_fake_annotation_version`.
- `mme/tests/fakes.py:make_term(term_id, service, index, name)` → unsaved `OntologyTerm`.

### Classifications
- `classification/models/classification.py:Classification.create(user, lab, lab_record_id=None, data=None, save=True, source=SubmissionSource.VARIANT_GRID, make_fields_immutable=False, populate_with_defaults=False, **kwargs)` → `Classification` with its first `ClassificationModification`; `data` is `{evidence_key: {"value": ...}}`. Marked deprecated in favour of `create_with_response`, but it is what every test uses. Needs the lab/user from `ClassificationTestUtils`.
- `mme/tests/fakes.py:make_classification(lab, user, gene_symbol="BRCA1", share_level=ShareLevel.PUBLIC.value, clinical_significance=ClinicalSignificance.VUS, withdrawn=False, is_last_published=True)` → `Classification` with a deterministic latest published modification; defaults are MME-eligible, vary one argument to break one rule.
- `classification/fake_reclassifications.py:FakeReclassifications` is the `manage.py create_fake_data reclassifications` subcommand (labs, curators, allele infos, classifications with multi-year histories); dev-box data for pages rather than a unit-test builder. `--delete` removes it.

### Analyses, nodes, tags
- `analysis/tests/utils.py:AnalysisSetupMixin` → mix into a `TestCase`: its `setUpTestData` sets `cls.grch37`, calls `get_fake_annotation_version`, and saves `cls.analysis` for a per-class user. Extend, then add samples/nodes.
- By hand: `analysis = Analysis(genome_build=grch37); analysis.set_defaults_and_save(user)` (needs an `AnnotationVersion` for the build). Nodes are plain creates, `SampleNode.objects.create(analysis=analysis, sample=sample)`, `TrioNode.objects.create(analysis=..., trio=..., inheritance=...)`; chain with `child.add_parent(parent)` then `child.save()` (see `_create_child_node` in `analysis/tests/test_models.py`).
- `analysis/tests/inheritance_node_mixin.py:InheritanceNodeTestsMixin` → shared inheritance-mode assertions for `TrioNode` / `QuadNode`; subclass sets `INHERITANCE`, `NODE_CLASS`, `REQUIRE_PARENT_ZYGOSITY`, `SHARED_ZYGOSITIES` and calls `cls.create_shared_variants(cgc)`. `DuoNode` has its own `analysis/tests/test_duo_node.py` - one parent means a different mode set and a `parent`-keyed zygosity table.
- Tags: `Tag.objects.get_or_create(pk="artefact")[0]` then `VariantTag.objects.create(genome_build=, analysis=, variant=, tag=, user=)` (`analysis/tests/test_variant_tags.py`). `analysis/fake_variant_tags.py:FakeVariantTags` is the `create_fake_data tags` subcommand for dev boxes.
- `snpdb/tests/utils/tag_testing_utils.py:create_classify_queue_tag(tag_id="ToDo", bucket=AlleleOriginBucket.UNKNOWN)` → a live `Tag` with `requires_classification` set, for anything about the classify queue: the to-do lists, the "New Classification" buttons, tag resolution. The default bucket is "Both"; pass `AlleleOriginBucket.SOMATIC` for the somatic path. Name the tag for the property under test rather than reaching for the seeded `settings.TAG_REQUIRES_CLASSIFICATION` name.

### Mocks and recorded data
- `snpdb/tests/utils/mock_clingen_api.py:MockClinGenAlleleRegistryAPI` / `MockServerErrorClinGenAlleleRegistryAPI` (always 502); `genes/tests/utils/mock_transcript_sequence_retrieval.py:MockTranscriptSequenceFetcher`; `classification/tests/utils/data_utils.py:ConditionMock` (its `setUp` is flagged FIXME - build ontology data with `create_ontology_test_data` instead); `mme/tests/fakes.py` `Fake*` dataclasses for MME profiles without the DB.
- `genes/tests/utils/hgvs_corpus.py:load_hgvs_corpus()` / `all_hgvs()` → the HGVS strings in `genes/tests/test_data/hgvs_corpus.tsv`.

## Base classes and helpers

`library/django_utils/unittest_utils.py:URLTestCase` is the base for every `<app>/tests/test_urls.py`. A test is a
list of `(url_name, reverse kwargs, expected status)`; a `"GET_PARAMS"` key in the kwargs becomes the query string:

```python
class Test(URLTestCase):
    def testUrls(self):
        URL_NAMES_AND_KWARGS = [("cohorts", {}, 200), ("view_user", {"pk": self.user_owner.pk}, 200)]
        self._test_urls(URL_NAMES_AND_KWARGS, self.user_non_owner)
        self._test_urls(URL_NAMES_AND_KWARGS, expected_code_override=302)  # GlobalLoginRequiredMiddleware bounces anon
```

- `_test_datatable_urls(...)` hits each URL and again with `?dataTableDefinition=1`; `_test_autocomplete_urls(names_obj_kwargs, user, in_results)` and `_test_datatables_grid_urls_contains_objs(names_obj, user, in_results)` assert an object's pk is (or is not) in the JSON - the standard permission check for grids.
- `URLTestCase` overrides settings: plain `StaticFilesStorage` (manifest storage would demand a staticfiles manifest that collectstatic never writes), `CELERY_TASK_ALWAYS_EAGER=True` (tasks run inline), `ANNOTATION_CACHED_WEB_RESOURCES=[]` (nothing auto-loads), `LIFTOVER_CLASSIFICATIONS=False`, `GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID=None`, `LOG_PARTITION_WARNINGS=False`. A plain `TestCase` gets none of these; add `@override_settings(CELERY_TASK_ALWAYS_EAGER=True)` yourself when a test fires a task.
- `prevent_request_warnings` decorates a test that expects 404/405 so `django.request` stays quiet.
- Query budgets: wrap the call in `CaptureQueriesContext(connection)` and assert on `production_query_count(ctx.captured_queries)`, which drops savepoints and reads of `PRODUCTION_CACHED_TABLES` (GenomeBuild, GeneSymbol, FlagType, ResolvedVariantInfo, Allele, Organization, Lab - their managers cache in production but not under `UNIT_TEST`). Raw `assertNumQueries` over-counts by exactly those.
- `VG_QUERY_PROFILE=/path.jsonl manage.py test ...` swaps in `QueryProfilingClient`: one JSON line per GET with status, query count, ms and duplicated (N+1) SQL. Add `VG_QUERY_TRACE=<regex>` to log a stack trace to `<path>.trace` for each matching statement.

## External services in tests

`variantgrid/test_runner.py:VariantGridTestRunner.setup_test_environment` sets `override_class` on two clients, so
tests never reach the network:

- ClinGen Allele Registry → `MockClinGenAlleleRegistryAPI`. HGVS lookups come from `snpdb/tests/test_data/clingen_hgvs_responses.json` (verbatim API JSON keyed by HGVS string), canonical alleles from `CANONICAL_ALLELES` in the mock. A miss raises `ValueError` naming the file to add the recording to; record it against the real registry and paste the JSON in.
- RefSeq/Ensembl transcript sequences → `MockTranscriptSequenceFetcher`, served from `genes/tests/test_data/transcript_sequences_{refseq,ensembl}.fasta`. A miss raises `BadTranscript` naming the FASTA; fetch with `TranscriptVersionSequenceInfo.get(accession)` on a real deployment and append `>accession` + sequence.

Reference genome fastas are read for real (HGVS resolution, VCF normalisation). CI uses the sparse copies in
`variantgrid/data/reference/sparse_test_fastas/` (all `N` except recorded regions; `variantgrid/settings/env/github_actions.py`
points there). A test that passes locally but sees `N` bases on CI fetches a region they lack: rerun the suite with
`--testrunner=variantgrid.test_runner.FastaRecordingRunner` on a box with the real fastas, then
`scripts/generate_sparse_test_fastas.py /tmp/vg_fasta_regions.jsonl variantgrid/data/reference/sparse_test_fastas`.

## Slow tests / traps

- `UNIT_TEST = sys.argv[1:2] == ['test']` (`default_settings.py`). Under it: the cache is in-process `LocMemCache`, the celery broker is `memory://`, Postgres JIT is off, axes and rollbar are off, passwords hash with MD5 rather than PBKDF2 (~1 s a hash, and the suite creates ~75 users), `ObjectManagerCachingImmutable` / `ObjectManagerCachingRequest` stop caching, and `admin_bot()` is uncached. Nothing in the suite should touch the dev Redis or RabbitMQ; if it does, a setting override is leaking.
- `CELERY_TASK_ALWAYS_EAGER` is `False` in `celery_settings.py`; with the memory broker a `.delay()` in a plain `TestCase` is enqueued and never runs. Use `URLTestCase` or override it.
- `get_fake_annotation_version` is the expensive fixture (ontology import, several versions): call it in `setUpTestData`, once per build, never per test. The runner seeds GRCh37 and GRCh38 (see Running tests), so those two builds are cheap; a test running alone still pays for the first call.
- Bulk variant builders are named `slowly_*` for a reason: fine for a handful of records, wrong for thousands.
- `snpdb/tests/test_fasta_index.py` loads the genome fasta `.fai`; `upload/tests/vcf/test_vcf_preprocess.py` writes a placeholder fasta and only asserts the command line (the bcftools pipe itself is not run in tests).
- Skipped on purpose: `genes/tests/test_hgvs.py` "Needs Ensembl contigs", and four `classification/tests/views/test_classification_view.py` cases pending variantgrid_private#3740.
