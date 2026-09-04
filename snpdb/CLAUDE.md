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
- Genotypes are packed one row per variant per cohort (models/models_cohort.py:CohortGenotype samples_* arrays, indexed by CohortSample.cohort_genotype_packed_field_index). Query them with models/models_cohort.py:CohortGenotypeCollection.get_annotation_kwargs and get_zygosity_q.
- Each CohortGenotypeCollection (and VariantZygosityCountCollection, VariantCollection) is its own partition table via library/django_utils/django_partition.py:RelatedModelsPartitionModel — create_partition on save, delete_related_objects to drop.
- Resolve which build a request is for with genome_build_manager.py:GenomeBuildManager.get_current_genome_build (GET param, URL path, user default, first annotated build, in that order).
- Read user preferences through models/models_user_settings.py:UserSettings.get_for_user — Global, Organization, Lab then User overrides, later wins.
- Search handlers register with search.py:search_receiver (see signals/variant_search.py); every other receiver is connected in apps.py:SnpdbConfig.ready, not at import.

Gotchas:
- Some contigs are shared between builds (MT, unplaced scaffolds), so one Variant can carry a VariantAllele per build; filter the variantallele join by genome_build or rows duplicate (grids.py:AbstractVariantGrid.get_initial_queryset, #1626).
- Gene-level events (fusions) are Variants on a fake contig with a gene id for a position; guard coordinate code with models/models_variant.py:Variant.get_gene_level_q and read gene_level_variants.py first.
- GenomeBuild, Allele, Lab and Organization managers cache lookups in production only (library/django_utils/django_object_managers.py:ObjectManagerCachingImmutable / ObjectManagerCachingRequest, off under settings.UNIT_TEST); expect stale instances, and count queries in tests with library/django_utils/unittest_utils.py:production_query_count.
- Allele.grch37 / grch38 / variants are cached_property; refetch the Allele after a liftover or merge. models/models_variant.py:Allele.merge refuses when both sides already have a ClinGenAllele.
- Adding or removing a CohortSample bumps models/models_cohort.py:Cohort.increment_version, which renumbers packed indexes for custom cohorts and orphans old CohortGenotypeCollections; go through the model, never bulk-update around it.
- Permissions cascade from the VCF: models/models_cohort.py:Cohort.can_view defers to the VCF, and deleting a VCF deletes its cohort (models/models_vcf.py:vcf_pre_delete_handler).
- Deleting a VCF or Sample from a view is a soft delete: tasks/soft_delete_tasks.py:soft_delete_vcfs sets ImportStatus.MARKED_FOR_DELETION and a celery task removes it; filter_for_user hides those rows in the meantime.
- models/models_vcf.py:VCF.delete_internal_data keeps the VCF and Sample rows and drops or recreates the partitions (recreate_partitions=False is the archive path in archive.py).
- Whole-table work on snpdb_variant or snpdb_allele is millions of rows in prod: page by pk range and fan out celery tasks (tasks/liftover_tasks.py:liftover_allele_batch, settings.LIFTOVER_BATCH_SIZE) rather than iterating one queryset.
- Liftover is per Allele, not per Variant: liftover.py:create_liftover_pipelines batches AlleleLiftover records and liftover.py:allele_can_attempt_liftover decides eligibility. Builds sharing a contig link with AlleleConversionTool.SAME_CONTIG and no external call.
- ClinGen Allele Registry calls are network I/O (clingen_allele.py:populate_clingen_alleles_for_variants); models/models_variant.py:Variant.can_have_clingen_allele bounds what may be sent.
- Views are split by topic (views/views_data.py, views_cohort.py, views_lab.py, views_user_settings.py, views_liftover.py, …); views/views.py holds only index, wiki and genome build/contig pages.

Tests:
- Fixture builders: tests/utils/fake_cohort_data.py:create_fake_cohort / create_fake_trio / create_fake_quad / create_fake_pedigree build VCF + samples + cohort in one call.
- Variants: tests/utils/vcf_testing_utils.py:slowly_create_test_variant and create_mock_allele; annotation/fake_annotation.py:get_fake_annotation_version for anything that touches annotation.
- Never hit the real registry: inject tests/utils/mock_clingen_api.py:MockClinGenAlleleRegistryAPI (or its ServerError sibling) via clingen_api=, as test_liftover.py and test_clingen_allele.py do.
- URL coverage: tests/test_urls.py:Test (URLTestCase); add new views to URL_NAMES_AND_KWARGS or PRIVATE_OBJECT_URL_NAMES_AND_KWARGS, datatable endpoints to testDataGridUrls.
- Query-count guards: tests/test_query_counts.py (view_sample must stay flat as trios grow).
- Slow: test_partition_archive_task.py and test_admin_archive_action.py run pg_dump into a temp dir; run them alone with --keepdb.

Deep reference: __snpdb_readme.md · claude/research/snpdb.md · claude/maps/models.md#snpdb
