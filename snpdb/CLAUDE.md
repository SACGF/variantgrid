# snpdb — agent notes
Owns: Variant/Allele/Locus/Sequence, GenomeBuild/Contig, VCF/Sample/Cohort/Trio, CohortGenotype packing, zygosity counts, liftover, ClinGen alleles, Lab/Organization, UserSettings, VariantGridColumn/CustomColumnsCollection, DataTables engine.
Start with:
- models/models_variant.py — Sequence, Locus, Variant, VariantCoordinate, Allele, VariantAllele, LiftoverRun, AlleleLiftover
- models/models_genome.py — GenomeBuild, Contig, GenomeBuildContig, GenomeFasta
- models/models_vcf.py — VCF, Sample, VCFFilter; models/models_cohort.py — Cohort, CohortSample, CohortGenotypeCollection, CohortGenotype, Duo, Trio, Quad
- liftover.py — allele liftover pipelines; clingen_allele.py — ClinGen Allele Registry linking
- views/datatable_view.py — DatatableConfig, RichColumn, DatabaseTableView; grids.py — AbstractVariantGrid and the list grids
- Import models from `snpdb.models` (the package re-exports every models_*.py); `snpdb/models.py` exists only for PyCharm.

Patterns here:
- Variant is build-specific, Allele is the build-independent hub. A Variant has no genome_build FK: its build comes from locus.contig via GenomeBuildContig (models/models_variant.py:Variant.genome_builds). Cross-build work goes through VariantAllele (models/models_variant.py:Allele.variant_for_build), never by comparing Variants.
- Restrict variant querysets to a build with models/models_variant.py:Variant.get_contigs_q (an IN list on contig ids); joining through GenomeBuildContig wrecks the planner's row estimate (#1720).
- Reference variants (alt == "=", models/models_variant.py:Variant.REFERENCE_ALT) exist on purpose; exclude them with models/models_variant.py:Variant.get_no_reference_q when you mean real calls.
- Canonicalise coordinates before lookup or insert with models/models_variant.py:VariantCoordinate.as_internal_canonical_form (alts >= settings.VARIANT_SYMBOLIC_ALT_SIZE become <DEL>/<DUP>/<INV> with svlen); models/models_variant.py:Variant.qs_from_variant_coordinate does this for you. Variant is unique on (locus, alt, svlen).
- Bulk-insert variants through variant_pk_lookup.py:VariantPKLookup (hash to pk, COPY of unknowns in batch_check); tests/utils/vcf_testing_utils.py:slowly_create_test_variant is the one-at-a-time test version.
- Always save Sequence via the model: models/models_variant.py:Sequence.save fills seq_sha256_hash, and the unique constraint is on the hash, not on seq.
- Three per-sample facts on models/models_vcf.py:VCF (delegated by Sample and Cohort): has_sample_columns (FORMAT + sample columns; drives the importer choice and grid sample columns), has_genotype (a GT field, so zygosity means something) and has_depth (AD/DP, so thresholds and allele frequency mean something). A depth-only caller such as TSO 500 splice variants has sample columns and depth but no genotype - every zygosity is stored as unknown, so nothing should filter on it.
- Genotypes are packed one row per variant per cohort (models/models_cohort.py:CohortGenotype samples_* arrays, indexed by CohortSample.cohort_genotype_packed_field_index). Query them with models/models_cohort.py:CohortGenotypeCollection.get_annotation_kwargs and get_zygosity_q.
- Each CohortGenotypeCollection (and VariantZygosityCountCollection, VariantCollection) is its own partition table via library/django_utils/django_partition.py:RelatedModelsPartitionModel — create_partition on save, delete_related_objects to drop.
- Resolve which build a request is for with genome_build_manager.py:GenomeBuildManager.get_current_genome_build (GET param, URL path, user default, first annotated build, in that order).
- `settings.TAG_REQUIRES_CLASSIFICATION` is seed data - the name a fresh install's classify queue tag gets. What behaves as one is `Tag.requires_classification`, set per tag on the tag settings page; ask `models/models.py:Tag.classify_queue_qs` (or `classify_queue_qs_for_bucket`), never a tag name.
- Read user preferences through models/models_user_settings.py:UserSettings.get_for_user — Global, Organization, Lab then User overrides, later wins.
- Search handlers register with search.py:search_receiver (see signals/variant_search.py); every other receiver is connected in apps.py:SnpdbConfig.ready, not at import.

Gotchas:
- Some contigs are shared between builds (MT, unplaced scaffolds), so one Variant can carry a VariantAllele per build; filter the variantallele join by genome_build or rows duplicate (grids.py:AbstractVariantGrid.get_initial_queryset, #1626).
- A Variant has at most one Allele per build - unique_together on (variant, genome_build) since #1361. Two Alleles wanting the same variant/build is resolved with models/models_variant.py:Allele.merge, which also sends models/models_variant.py:allele_merged_signal so classification re-homes the clinical contexts and groupings it moved.
- Gene-level events (fusions) are Variants on a fake contig with a gene id for a position; guard coordinate code with models/models_variant.py:Variant.get_gene_level_q and read gene_level_variants.py first.
- GenomeBuild, Allele, Lab and Organization managers cache lookups in production only (library/django_utils/django_object_managers.py:ObjectManagerCachingImmutable / ObjectManagerCachingRequest, off under settings.UNIT_TEST); expect stale instances, and count queries in tests with library/django_utils/unittest_utils.py:production_query_count.
- Allele.grch37 / grch38 / variants are cached_property; refetch the Allele after a liftover or merge. models/models_variant.py:Allele.merge refuses when both sides already have a ClinGenAllele.
- Only a VCF's own cohort has the FORMAT/INFO JSON: tasks/cohort_genotype_tasks.py packs a custom cohort's
  CohortGenotype rows from the source collections and writes `format`/`info` as empty, so anything read out
  of them (grid_columns/grid_sample_columns.py:get_copy_number_annotation) is blank for a custom cohort.
- Adding or removing a CohortSample bumps models/models_cohort.py:Cohort.increment_version, which renumbers packed indexes for custom cohorts and orphans old CohortGenotypeCollections; go through the model, never bulk-update around it. For a whole new membership use models/models_cohort.py:Cohort.set_samples - it applies the diff under one version bump (and so one CohortGenotype rebuild) instead of one per row.
- Permissions cascade from the VCF: models/models_cohort.py:Cohort.can_view defers to the VCF, and deleting a VCF deletes its cohort (models/models_vcf.py:vcf_pre_delete_handler).
- Deleting a VCF or Sample from a view is a soft delete: tasks/soft_delete_tasks.py:soft_delete_vcfs sets ImportStatus.MARKED_FOR_DELETION and a celery task removes it; filter_for_user hides those rows in the meantime.
- models/models_vcf.py:VCF.delete_internal_data keeps the VCF and Sample rows and drops or recreates the partitions (recreate_partitions=False is the archive path in archive.py).
- Whole-table work on snpdb_variant or snpdb_allele is millions of rows in prod: page by pk range and fan out celery tasks (tasks/liftover_tasks.py:liftover_allele_batch, settings.LIFTOVER_BATCH_SIZE) rather than iterating one queryset.
- Liftover is per Allele, not per Variant: liftover.py:create_liftover_pipelines batches AlleleLiftover records and liftover.py:allele_can_attempt_liftover decides eligibility. Builds sharing a contig link with AlleleConversionTool.SAME_CONTIG and no external call.
- Liftover has left an Allele with 2 VariantAlleles in the destination build - the real match plus a
  "build difference" variant that other Alleles at the locus also picked up. That's what a variant with 2
  ClinGen Alleles is; `one_off_dedupe_variant_alleles` drops the artefact link rather than merging them.
- A tool that has ever errored on an allele/build is skipped forever (models/models_variant.py:AlleleLiftover.get_failed_conversion_tools),
  which is why re-clicking liftover does nothing. Pass `retry_conversion_tools=` to liftover.py:create_liftover_pipelines to
  override it for chosen tools - the liftover page's per-tool retry buttons, "Create Variant" and the admin action all do (#1273).
- ClinGen Allele Registry calls are network I/O (clingen_allele.py:populate_clingen_alleles_for_variants); models/models_variant.py:Variant.can_have_clingen_allele bounds what may be sent.
- Somalier decides how to genotype from the header of the VCF we hand it: a FORMAT `AD` line means it re-genotypes every sample from the depths and applies its own QC at relate time (`--min-depth` 7, `--min-ab` 0.3), no `AD` line means it trusts `GT`. variants_to_vcf.py:vcf_export_to_file picks one per VCF - declaring `AD` we can't fill in zeroes out every sample (#183).
- somalier keeps each site's two alleles alphabetically and reads the genotype against that pair rather than against the record's `REF`/`ALT` (brentp/somalier#163), so variants_to_vcf.py:vcf_export_to_file writes the `AD` pair in the site's order - `REF`, `ALT` and `GT` stay as called, and only a VCF with no depths flips `GT` instead. Getting it wrong is invisible in a jointly called VCF - it cancels out pairwise - and shows up as inflated relatedness once `--unknown` is in play (#183). `settings.SOMALIER["compensate_allele_order"]` turns it off for a somalier that reads the record's alleles; `deployment_check`'s `somalier_allele_order` genotypes both allele orders through the installed binary and fails saying which way to set it, so nobody has to remember (variantgrid/deployment_validation/somalier_check.py). Changing it means rebuilding the extracts.
- `CohortGenotype` stores "no value" as `-1` (models/models_cohort.py:CohortGenotype.MISSING_NUMBER_VALUE), so `is not None` is not enough when reading `samples_allele_depth` / `samples_read_depth` / `samples_allele_frequency`.
- Views are split by topic (views/views_data.py, views_cohort.py, views_lab.py, views_user_settings.py, views_liftover.py, …); views/views.py holds only index, wiki and genome build/contig pages.

Tests:
- Fixture builders: tests/utils/fake_cohort_data.py:create_fake_cohort / create_fake_trio / create_fake_quad / create_fake_pedigree build VCF + samples + cohort in one call.
- Variants: tests/utils/vcf_testing_utils.py:slowly_create_test_variant and create_mock_allele; annotation/fake_annotation.py:get_fake_annotation_version for anything that touches annotation.
- Never hit the real registry: inject tests/utils/mock_clingen_api.py:MockClinGenAlleleRegistryAPI (or its ServerError sibling) via clingen_api=, as test_liftover.py and test_clingen_allele.py do.
- URL coverage: tests/test_urls.py:Test (URLTestCase); add new views to URL_NAMES_AND_KWARGS or PRIVATE_OBJECT_URL_NAMES_AND_KWARGS, datatable endpoints to testDataGridUrls.
- Query-count guards: tests/test_query_counts.py (view_sample must stay flat as trios grow).
- Slow: test_partition_archive_task.py and test_admin_archive_action.py run pg_dump into a temp dir; run them alone with --keepdb.

Deep reference: __snpdb_readme.md · claude/research/snpdb.md · claude/maps/models.md#snpdb
