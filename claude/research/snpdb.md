# snpdb — research notes

Verified against 7c4408c62 on 2026-09-06

snpdb is the core app: the genome (`GenomeBuild`, `Contig`), the variant identity model (`Sequence`, `Locus`,
`Variant`, `Allele`), samples and their packed genotypes (`VCF`, `Sample`, `Cohort`, `CohortGenotypeCollection`),
liftover between builds, ClinGen Allele Registry linking, labs and organisations, user settings, the grid column
catalogue and the DataTables engine every list page uses. The rules are in `snpdb/CLAUDE.md`; the models are listed
in `claude/maps/models.md#snpdb`; this is the story behind them. The name is historical ("SNP database" predates
VariantGrid); nothing here is SNP-specific.

## Flows

### Inserting a variant

Every VCF line becomes one `Locus` (contig, position, ref) and one `Variant` per alt. Because `snpdb_variant` is
tens of millions of rows, insertion goes through `snpdb/variant_pk_lookup.py:VariantPKLookup`, which hashes each
coordinate to a tuple `(contig_id, position, ref_id, alt_id, svlen)` in `VariantPKLookup._get_variant_hash` and
resolves hashes to pks per contig in SQL (`VariantPKLookup._get_variant_ids`, a `Concat` annotation looked up with
`locus__position__in`). `null` svlen is encoded as `''` because 0 is a real value gene-level variants rely on.
`VariantPKLookup.add` canonicalises first (`VariantCoordinate.as_internal_canonical_form`) and short-circuits a
coordinate whose ref or alt sequence is not yet in `sequence_pk_by_seq` straight to the unknown list without a query;
`VariantPKLookup.batch_check` resolves the rest and, with `insert_unknown=True`, calls `VariantPKLookup._insert_unknown`,
which inserts sequences (`bulk_create(ignore_conflicts=True)` then a re-read, since bulk_create returns no pks), then
loci and variants by Postgres `COPY` (`upload/vcf/sql_copy_files.py:loci_sql_copy_csv`, `variants_sql_copy_csv`).
`batch_check` must stay on the calling thread: the final `insert_all=True` pass depends on its side effects.

The insert runs on the `variant_id_single_worker` queue (`variantgrid/settings/components/celery_settings.py`) in
`upload/tasks/vcf/unknown_variants_task.py:InsertUnknownVariantsTask`, one process so there is exactly one row per
coordinate; a belt-and-braces cache lock inside `InsertUnknownVariantsTask.process_items` covers the case where two
copies of the task are ever started. The producer is `upload/tasks/vcf/unknown_variants_task.py:SeparateUnknownVariantsTask`,
run per split VCF file. `Sequence.save` (`snpdb/models/models_variant.py:Sequence.save`) fills `seq_sha256_hash`;
the unique constraint is on the hash because a text index row is capped at 8 kB and 10 kB substitutions exist. The
bulk path sets the hash by hand for the same reason.

### Coordinates, symbolic alts and canonical form

`snpdb/models/models_variant.py:VariantCoordinate` is a pydantic value object with two exit doors: `as_external_explicit`
for VCF and HGVS, `as_internal_symbolic` for the database. `VariantCoordinate.as_internal_symbolic` turns an alt at or
beyond `settings.VARIANT_SYMBOLIC_ALT_SIZE` (1000 in `variantgrid/settings/components/default_settings.py`) into
`<DUP>` (single-base ref, HGVS says dup), `<DEL>` (single-base alt equal to ref[0]) or `<INV>` (equal lengths, reverse
complement) with an `svlen`. `VariantCoordinate.as_internal_canonical_form` is the rule "one representation per
variant": symbolic only when `abs(svlen)` clears the threshold, otherwise explicit if `VariantCoordinate.can_be_made_explicit`
(`<CNV>` and `<INS>` cannot be), and an alt equal to the ref becomes `Variant.REFERENCE_ALT`. `Variant.qs_from_variant_coordinate`
applies this and resolves the contig up front from `GenomeBuild.chrom_contig_mappings` so the filter lands on the
leading edge of the `(contig, position, ref)` unique index instead of joining GenomeBuildContig (#1720). `Variant` is
unique on `(locus, alt, svlen)`: two CNVs can share an alt and differ only in length. `Variant.end` is stored, not
computed, for overlap queries.

### Linking a Variant to an Allele

`Allele` is the build-independent hub; `VariantAllele` is the per-build link with an `origin`
(`snpdb/models/models_enums.py:AlleleOrigin`: imported, imported-normalised, liftover) and the `allele_linking_tool`
that made it (`AlleleConversionTool`). `snpdb/clingen_allele.py:populate_clingen_alleles_for_variants` is the batch
path: variants that fail `Variant.can_have_clingen_allele` (gene-level, an alt that cannot form g.HGVS, or larger than
`ClinGenAllele.CLINGEN_ALLELE_MAX_ALLELE_SIZE`, with `<DUP>` counted double for a registry bug) get a bare `Allele`
with no network call, once - a variant that can never be registered counts as done on every later call, which is what
stopped it collecting a fresh `Allele` per call (#1361 / #1844). The rest are sent as one `hgvs_put` batch, and the
response creates `ClinGenAllele` rows or stores the error on `VariantAllele.clingen_error`. Alleles are created with
`ignore_conflicts=True` and re-read because the Allele may already exist from another build; any empty Allele whose
`VariantAllele` insert was dropped as a conflict is deleted at the end of the call. The single-variant path is `snpdb/clingen_allele.py:get_variant_allele_for_variant`,
and `snpdb/clingen_allele.py:variant_allele_clingen` handles the registry saying "this coordinate is CA123" when CA123
already belongs to a different Allele by calling `Allele.merge`.

`snpdb/models/models_variant.py:Allele.merge` moves the ClinGen id, flags, clinical contexts, classifications,
variant tags, `ImportedAlleleInfo`s, `AlleleLiftover`s, `ClinVarRecordCollection`s and every `VariantAllele` onto the
survivor (deleting a link that would violate the unique constraint) and refuses when both sides already carry a
ClinGenAllele, logging the attempt in `AlleleMergeLog`. It ends by sending
`snpdb/models/models_variant.py:allele_merged_signal`, which
`classification/signals/classification_hooks_allele_merge.py:allele_merged_handler` uses to put the moved
classifications' clinical contexts and groupings back under the surviving Allele - both are derived records that are
otherwise only re-homed when a classification is next published.
`ClinVarAllele` and `AlleleGrouping` are left alone: both are rebuilt from classifications and both carry uniqueness a
blind update would violate. `VariantAllele.needs_clingen_call` only retries a
stored error when it was a server error, and only when `settings.CLINGEN_ALLELE_REGISTRY_LOGIN` is set.

### Liftover

Liftover is per Allele, not per Variant: the goal is a Variant in every annotated build for each Allele.
`snpdb/liftover.py:create_liftover_pipelines` pages the allele queryset by pk (`snpdb/liftover.py:_batch_alleles`, so no
cursor is held open for the hours a big run takes and alleles lifted over by earlier batches drop out), holds one batch
of `settings.LIFTOVER_BATCH_SIZE` in memory and writes one VCF per (batch, tool, destination build). The celery entry
`snpdb/tasks/liftover_tasks.py:liftover_alleles` fans out `snpdb/tasks/liftover_tasks.py:liftover_allele_batch` per pk
range, and each batch re-queries `Allele.missing_variants_for_build` so a failure only loses its own batch.

`snpdb/liftover.py:_get_build_liftover_dicts` chooses the method per allele, stopping at the first that works:
`snpdb/liftover.py:_liftover_using_existing_contig` when the contig is shared between builds (MT), which just creates the
`VariantAllele` with `AlleleConversionTool.SAME_CONTIG` and no VCF at all
(`snpdb/liftover.py:_run_liftover_using_same_contig` - if the destination build's variant already has an Allele of its
own it merges the two and records the `AlleleLiftover` as `SKIPPED` instead); `snpdb/liftover.py:_liftover_using_dest_variant_coordinate`, which
asks the ClinGen record for the destination g.HGVS (`ClinGenAllele.get_g_hgvs`) and, if `settings.LIFTOVER_DBSNP_ENABLED`
(off by default), dbSNP; and `snpdb/liftover.py:_liftover_using_source_variant_coordinate`, whose only option is
`AlleleConversionTool.BCFTOOLS_LIFTOVER` (`snpdb/bcftools_liftover.py:bcftools_liftover`). `PICARD` and `CROSSMAP` exist
in the enum but have no code path. The VCF ID column is the Allele pk, only standard contigs are written
(`snpdb/liftover.py:_non_standard_contig_error` records a per-allele error instead of failing the run, #1197), and the
file goes through the normal upload pipeline as `upload/import_task_factories/import_task_factories.py:LiftoverImportFactory`.
Results land in `upload/vcf/bulk_allele_linking_vcf_processor.py:BulkAlleleLinkingVCFProcessor.batch_handle_variant_ids`,
which sets each `AlleleLiftover.status` and merges alleles into the lowest pk when the destination variant is already
linked elsewhere (`BulkAlleleLinkingVCFProcessor.merge_alleles`), so `merge(a, b)` and `merge(b, a)` racing agree.

`LiftoverRun` has no status of its own; status is per allele on `snpdb/models/models_variant.py:AlleleLiftover` using
`snpdb/models/models_enums.py:ProcessingStatus`. `AlleleLiftover.get_failed_conversion_tools` is read once per batch
so a tool that already failed an allele is not retried, and `LiftoverRun.get_clingen_auto_fail_liftover_run` collects the
"known to fail before we start" cases into one run per build so the failure is still recorded.

### Cohorts and packed genotypes

Genotypes are not stored per sample. `snpdb/models/models_cohort.py:CohortGenotype` is one row per (collection,
variant) with the samples packed into parallel arrays and a one-character-per-sample `samples_zygosity` string, indexed
by `CohortSample.cohort_genotype_packed_field_index`. Every VCF gets an automatic Cohort; custom cohorts pick samples
across VCFs. `snpdb/tasks/cohort_genotype_tasks.py:cohort_genotype_task` builds the rows with raw SQL
(`snpdb/tasks/cohort_genotype_tasks.py:_get_insert_cohort_genotype_sql`): for each sample it left-joins that sample's own
VCF-cohort partition and concatenates `coalesce(partition.col[i], empty)` per packed column. A custom cohort whose
samples all sit inside one VCF cohort becomes a sub-cohort sharing the parent's packing instead
(`snpdb/tasks/cohort_genotype_tasks.py:create_cohort_genotype_and_launch_task`).

Membership is versioned. `snpdb/models/models_cohort.py:Cohort.increment_version` bumps `version`, renumbers the packed
indexes of a custom cohort and marks older `CohortGenotypeCollection`s for deletion; `CohortSample.save` and
`CohortSample.delete` each call it, so per-row edits rebuild the packing per row. `Cohort.set_samples` (2026-09-05) applies
a whole new membership under one bump. `snpdb/models/models_cohort.py:CohortVersion` (#1551) gives version-specific
caches something to cascade from, the way `NodeVersion` does in analysis.

Each collection is its own partition table: `snpdb/models/models_cohort.py:CohortGenotypeCollection` is a
`library/django_utils/django_partition.py:RelatedModelsPartitionModel`, whose `save` runs
`CREATE TABLE … INHERITS` with a check constraint on the collection id and whose `delete_related_objects` drops it.
Since 2024-11 a VCF's genotypes are split in two: the UNCOMMON collection and a COMMON one holding variants above
`CohortGenotypeCommonFilterVersion.gnomad_af_min` in every listed gnomAD version, so a rare-variant filter can skip the
common partition entirely. `CohortGenotypeCollection.get_annotation_kwargs` only skips it when the analysis's gnomAD
version is one the partition was built for (#1119, #1582); `snpdb/common_variants.py:get_classified_high_frequency_variants_qs`
keeps classified common variants out of the common side.

Querying is a `FilteredRelation` on `cohortgenotype` plus a regex on `samples_zygosity`
(`CohortGenotypeCollection.get_zygosity_q`), inverted with a negative lookahead rather than `~Q()`. The sample-level
path uses `Substr` instead (`snpdb/models/models_vcf.py:Sample.get_annotation_kwargs`) because benchmarking (#1494) found
substring-plus-IN faster for one sample and regex faster for wide cohorts. `Cohort.get_any_sample_called_variant_collection`
(#1551) pre-computes a `VariantCollection` for sub-cohorts so the analysis EXCLUDE filter becomes a hash join.

### VCF and Sample lifecycle

`snpdb/models/models_vcf.py:VCF.import_status` and `Sample.import_status` share `snpdb/models/models_enums.py:ImportStatus`.
Permissions cascade from the VCF: `Sample.can_view` is the VCF's answer or a sample-level guardian grant, and
`Sample.filter_for_user` ORs both; `VCF.filter_for_user` also hides `ImportStatus.DELETION_STATES` and, by default,
archived VCFs. Deleting from the UI is `snpdb/tasks/soft_delete_tasks.py:soft_delete_vcfs`, which marks rows and queues
`snpdb/tasks/soft_delete_tasks.py:remove_soft_deleted_vcfs_task` on the single scheduling worker, one VCF at a time because
larger deletes blew out the commit log. `VCF.delete_internal_data` keeps the VCF and Sample rows and either recreates the
partitions (reload in place) or drops them (`snpdb/archive.py:archive_vcf`, #1536), stamping
`library/django_utils/data_archive_mixin.py:DataArchiveMixin`. Before any partition drop the archive pipeline writes a
pg_dump (`snpdb/partition_archive.py:archive_partitioned_model`, `snpdb/tasks/partition_archive_tasks.py:perform_partition_archive`,
#1537); `RelatedModelsPartitionModel._warn_if_no_archive` only warns, it does not block.

Per-deployment zygosity counts (`snpdb/models/models_zygosity_counts.py:VariantZygosityCountCollection`) are another
partition, updated by raw SQL in `snpdb/variant_zygosity_count.py:update_all_variant_zygosity_counts_for_vcf` on the
single worker, with `VariantZygosityCountForVCF` and `VariantZygosityCountForSample` as the audit rows that let a delete
subtract exactly what an import added. `data_version` on the collection is bumped in the same transaction as the last
write so analysis nodes can pin what they read. Per-sample stats moved to `snpdb/models/models_cohort_stats.py:CohortGenotypeStats`,
whose chrX het/hom ratio is the "detected sex" (`CohortGenotypeStats.chrx_sex_guess`).

`snpdb/models/models_vcf.py:VCFSourceSettings` rewrites a VCF's sample-field bindings by regex on its `source` header,
because callers reuse standard FORMAT ids for other meanings. `snpdb/signals/signal_handlers.py:backend_vcf_import_success_handler`
creates BED intersections for sequencing samples, and `trio_post_save_handler` launches karyomapping for a new Trio only
(Duo and Quad do not).

### Labs, users and the current build

`snpdb/models/models.py:Lab` and `Organization` each own a Django Group named by `group_name`, materialised on access
(`Lab.group`, `Organization.group`) and pre-created in `Lab.save`; lab membership *is* group membership, so
`Lab.add_member` also joins the organisation group and `Lab.remove_member` leaves it. `Organization.is_member` is
defined through `Lab.valid_labs_qs`. Membership changes are logged and notified (`Lab._log_membership_change`), and a lab
whose membership is decided outside VariantGrid says so through `Lab.external_membership`. `LabHead` is a separate row.

`snpdb/models/models_user_settings.py:UserSettings` is a dataclass, not a model: `UserSettings.get_for_user` layers
`GlobalSettings`, `OrganizationUserSettingsOverride`, `LabUserSettingsOverride` and `UserSettingsOverride`, later wins,
deriving the lab and organisation from the user's default lab when not given. `snpdb/genome_build_manager.py:GenomeBuildManager.get_current_genome_build`
resolves the build for a request in a fixed order (GET parameter, build name in the URL path, the user's default, the
first annotated build) and caches it on the request threadlocal; the URL regex only knows GRCh37 and GRCh38.
`snpdb/apps.py:SnpdbConfig.ready` connects `user_post_save_handler` only outside `UNIT_TEST`, imports every
`snpdb/signals/` module for its `@search_receiver` side effect, and connects the trio and VCF-import handlers.

### Grids

`snpdb/views/datatable_view.py:DatatableConfig` defines both the client table and the server response; every list page
subclasses it. `DatatableConfig.apply_filters` is the one chokepoint (config filter, then the search box via
`DatatableConfig.power_search`, then the client's column rules via `library/django_utils/filter_rules.py:rules_to_q`),
and `DatatableConfig.ordering` always appends `F("pk").desc()` so paging is stable. `DatatableConfig.known_count` and
`approximate_count` (#1700) let a page skip an exact count when a node count or a planner estimate is at hand.
`snpdb/views/datatable_view.py:DatabaseTableView` serves the JSON with `json_allow_nan = False` so a NaN annotation value
renders blank instead of breaking `JSON.parse`, and inherits `library/django_utils/major_operation.py:MajorOperationViewMixin`.

`snpdb/grids.py:AbstractVariantGrid` is the variant table (the analysis node grid and the standalone variant lists).
Columns are built per user from `UserSettings.columns` through `snpdb/grid_columns/custom_columns.py:get_variant_grid_columns`
from the catalogue in `snpdb/models/models_columns.py:VariantGridColumn` (composite columns since 2026-09-02).
`AbstractVariantGrid.get_initial_queryset` restricts the `variantallele` join to the grid's build (#1626, shared MT contig)
and applies `_get_q` once, because applying column rules there as well would add a second join per filtered relation.
`AbstractVariantGrid._genomic_order_by` orders contigs with a `CASE` over `GenomeBuild.standard_contigs` for the same
shared-contig reason.

### Search

`snpdb/search.py:SearchInput.search` sends `search_signal` and every `@search_receiver` (`snpdb/search.py:search_receiver`)
answers with a `SearchResponse`. The decorator does the shared work: admin-only gating, `preview_enabled`, the regex
`pattern`, result caps (variants uncapped because they merge into alleles), exception capture into a
`SearchMessageOverall`, and under `settings.PREFER_ALLELE_LINKS` the conversion of variant hits into allele hits.
`SearchInput.get_visible_variants` restricts to a build's contigs and, on Shariant, to classified variants.
`SearchResultMatchStrength` decides whether the UI jumps straight to a single result. The variant receivers live in
`snpdb/signals/variant_search.py` (locus, dbSNP, gnomAD, HGVS via `snpdb/signals/variant_search.py:_search_hgvs`, ids,
fusions as lookup only); the rest of the app registers cohort, trio, duo, quad, sample, VCF, lab, organisation and user
searches from `snpdb/signals/`.

## Why it is shaped this way

- **No genome_build on Variant.** A Contig can belong to more than one build (MT is identical in GRCh37 and GRCh38),
  so the build is a property of the locus's contig (`Variant.genome_builds`). The price is paid in three places: the
  `variantallele` filter in the grids, the CASE ordering, and `Variant.get_contigs_q` as an IN list because the
  GenomeBuildContig join collapses the planner's estimate to one row (#1720).
- **Allele as the hub.** Classifications, ClinGen ids, flags and liftover hang off the Allele so that the same change in
  two builds is one thing. Merging alleles is therefore a data-migration-grade operation with an audit log.
- **Packed genotypes.** One row per (cohort, variant) instead of one per (sample, variant) keeps the biggest tables an
  order of magnitude smaller and makes "any sample called" a regex or substring on one column. The cost is that any
  membership change is a full rebuild, which is why versions, sub-cohorts and `set_samples` exist.
- **Partition tables, dropped not deleted.** Deleting millions of rows is slower and bloats the WAL; dropping a child
  table is instant. Everything version- or collection-scoped is a partition, and archives exist so a drop is reversible.
- **Gene-level variants are a declared hack.** Fusions with no coordinate are Variants on one fake contig
  (`snpdb/gene_level_variants.py`, #1506) with the gene id as position, `svlen = 0` rather than null so the unique
  constraint still bites. `Variant.get_gene_level_q` is the one predicate; VCF writers are safe because they enumerate
  `GenomeBuild.standard_contigs`.

## History

- 2024-05: BCFtools liftover replaced NCBI Remap (variantgrid_private#2647); 2024-07 liftover code moved out of `Allele`
  into `snpdb/liftover.py` (#1014) and tuples became `VariantCoordinate`.
- 2024-11: common/rare split of CohortGenotypeCollection (variantgrid_private#3704); 2024-12 T2T-CHM13v2 (#814) and the
  symbolic `<DEL>` representation (#1214); 2025-04 `as_internal_canonical_form` (#1275).
- 2026-04/06: the gnomAD-version guard on the common partition (#1119, #1582); trio and quad search (#1623);
  build-restricted grids (#1626); 2026-05 data archive (#1536) and the pre-drop pg_dump pipeline (#1537); profiling that
  settled regex versus substring (#1494); sub-cohort VariantCollections and CohortVersion (#1551).
- 2026-08: per-allele liftover errors for non-standard contigs (#1197); known counts for grids (#1700); gene fusions as
  variants (#1506); contig IN lists (#1720); retired rather than deleted tags (#1751); jqGrid removed, variant grids on
  the native DataTables engine (#1785, #1815); symbolic CNV HGVS from coordinates (#1571).
- 2026-09: user awards (#1819); proband sex from chrX genotypes (variantgrid_sapath#439); `snpdb/views/views.py` split
  into topic modules; composite grid columns; Duo analyses (#1829); `Cohort.set_samples` and the cohort stats tab.

## Traps

- `Variant.is_symbolic` and `Variant.is_reference` are properties; `VariantCoordinate` methods are methods.
- `Allele.grch37` / `grch38` / `variants` are cached properties; refetch after liftover or merge.
- `GenomeBuild`, `Allele`, `Lab` and `Organization` managers cache in production only
  (`library/django_utils/django_object_managers.py:ObjectManagerCachingImmutable` / `ObjectManagerCachingRequest`).
  `GenomeBuild.get_name_or_alias` and `builds_with_annotation_cached` are also 60-second `timed_cache`s, and
  `SiteMessage.get_site_messages` is a 30-second Django cache entry.
- `GenomeFasta.get_for_genome_build` uses `get_or_create`, so an empty row created once is returned forever with its
  index never loaded.
- `AlleleLiftover.BULK_UPDATE_BATCH_SIZE` is 1000 because `bulk_update` writes one CASE arm per row.
- Search must never mint a fusion identity; `snpdb/signals/variant_search.py:search_variant_gene_fusion` is lookup only.
- Client-side grid renderers add links; they never convert data (`snpdb/grids.py:get_standard_overrides`).
- `snpdb/forms.py:VCFForm` hides unused sample-field columns, so not every model field is on the form.
- `ManualVariantEntry` lives in annotation (`annotation/models/models.py:ManualVariantEntryCollection`); snpdb only
  supplies its form, view and grid.
