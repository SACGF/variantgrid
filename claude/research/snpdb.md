# VariantGrid SNPDB App — Reference Document

## Purpose and Overview

The **snpdb** app is the core application of VariantGrid. Originally named for "SNP database" before the VariantGrid name was adopted, it contains the foundational models for:
- Genomic variant and allele management (build-specific and build-independent)
- VCF import and processing workflows
- Sample management with genotype data
- Cohort analysis (groupings of samples)
- Liftover between genome builds
- User and organizational management
- Permissions and access control
- Zygosity tracking and cohort genotype collections
- Custom column configurations and user settings
- Integration with ClinGen Allele Registry for variant normalization

---

## Key Models

### Organization and Lab

**Organization** — Company/institution level grouping
- Fields: name, short_name, group_name, classification_config, active
- Methods: is_member(), can_write(), get_classifying_labs(), get_sharing_labs()

**Lab** — Individual laboratory within organization
- Fields: name, city, country, url, contact info, organization FK, clinvar_key FK, group_name
- Methods: is_member(), can_write(), valid_labs_qs(), classifications properties
- Related to Classifications, ClinVarKey, LabProject

**LabHead** — Maps users to labs with lab head role (unique_together: lab, user)

**LabProject** — Lab-specific projects with leader, members, families, dates

### Genome Models

**GenomeBuild** — Reference genome builds (GRCh37, GRCh38, T2T-CHM13v2.0)
- Fields: name (PK), accession, alias (unique), slug, enabled, igv_genome
- Methods: grch37(), grch38(), t2tv2() class methods; get_name_or_alias() (fuzzy, cached); detect_from_filename(); builds_with_annotation()
- Caching: ObjectManagerCachingImmutable, timed_cache on get_name_or_alias()

**Contig** — Chromosome/scaffold representation
- Fields: name, refseq_accession, genbank_accession, ucsc_name, role (SequenceRole), molecule_type (AssemblyMoleculeType)

**GenomeBuildContig** — Junction table with ordering (genome_build FK, contig FK, order)

### Variant Models

**Sequence** — DNA sequences (reference and alternate alleles)
- Fields: val (sequence string)

**Locus** — Chromosomal location (contig FK, position, ref FK)
- Unique on (contig, position, ref)

**Variant** — Build-specific variant (locus FK, alt FK, genome_build FK)
- Unique on (locus, alt, genome_build)
- Partitioned via RelatedModelsPartitionModel

**VariantCoordinate** — Pydantic model for variant coords with validation
- Fields: chrom, position, ref, alt, svlen (optional for symbolic variants)

**Allele** — Build-independent variant representation
- Optional link to ClinGen Allele Registry via ClinGenAllele OneToOne
- Methods: variant_alleles(), grch37/grch38 (cached), variant_for_build()
- Implements FlagsMixin

**VariantAllele** — Junction: Allele ↔ Variant ↔ GenomeBuild
- Fields: variant, allele, genome_build FKs; origin (IMPORTED_TO_DATABASE, IMPORTED_NORMALIZED, LIFTOVER, LIFTOVER_NORMALIZED); allele_linking_tool, clingen_error

**AlleleMergeLog** — Audit trail for Allele.merge() operations

**LiftoverRun** — Tracks liftover operations (PENDING/PROCESSING/SUCCESS/ERROR)

**AlleleLiftover** — Individual allele liftover results with status and error FK

### VCF and Sample Models

**VCF** — Represents imported VCF files
- Fields: name, date, genome_build FK, project FK, user FK, genotype_samples count, import_status, header, source
- Field mappings: allele_depth_field, allele_frequency_field, ref/alt/read_depth fields, genotype_field, etc.
- Methods: filter_for_user(), get_sample_ids(), get_variant_qs(), delete_internal_data()
- Pre-delete signal cascades to cohort deletion

**Sample** — Individual sample from VCF
- Fields: vcf FK, vcf_sample_name, name, no_dna_control, research_consent, patient FK, specimen FK, import_status, variants_type (UNKNOWN/GERMLINE/MIXED/SOMATIC_ONLY)
- Methods: get_genotype(variant), get_variant_qs(), filter_for_user(), get_annotation_kwargs()
- Preview system with fa-microscope icon

**SampleFilePath** — Paths to BAM/CRAM/BED/VCF files (sample FK, file_path, file_type)

**VCFFilter** — VCF filter codes (PASS, LowQual, etc.)

**SampleStats, SampleStatsPassingFilter** — Variant counts per sample by annotation type

### Cohort Models

**Cohort** — Collection of samples for grouped analysis
- Fields: name, user FK, version, import_status, genome_build FK, vcf OneToOne FK (null for custom cohorts), parent_cohort FK (self-referential), sample_count
- Methods: add_sample(), get_samples(), create_sub_cohort(), increment_version()
- GuardianPermissionsAutoInitialSaveMixin (inherits permissions from VCF)

**CohortSample** — Junction between Cohort and Sample with packed field indexing
- Fields: cohort FK, sample FK, cohort_genotype_packed_field_index (0-based), sort_order
- Save/delete signals trigger cohort.increment_version()

**CohortGenotypeCollection** — Partitioned table for genotype data
- Packs genotype info into single row for fast multi-sample queries (similar to Gemini)
- Methods: get_annotation_kwargs(), get_zygosity_q(), get_packed_index_by_sample_id()
- Types: UNCOMMON (rare variants) and COMMON (high AF variants)

**CohortGenotype** — Actual genotype data (packed format per row)
- Fields: variant FK, collection FK, samples_zygosity (string), samples_allele_depth/read_depth/genotype_quality/phred_likelihood (arrays), samples_filters (array), vc_zygosity_count

**Trio** — Pedigree for Mendelian inheritance analysis
- Fields: name, user FK, cohort FK, mother/father/proband FKs to CohortSample, mother_affected, father_affected
- Post-save signal creates karyomapping analysis via Celery

**SampleGenotype** — Non-model class for genotype information extracted from packed fields

### Zygosity Count Models

**VariantZygosityCountCollection** — Partitioned container for zygosity counts across VCFs/samples

**VariantZygosityCount, VariantZygosityCountForVCF, VariantZygosityCountForSample** — Tracks count of each zygosity (HOM_REF/HET/HOM_ALT/MISSING) per variant

### User Settings Models

**UserSettings** — Aggregator (non-model) for cascading settings
- Cascade: GlobalSettings → OrganizationUserSettingsOverride → LabUserSettingsOverride → UserSettingsOverride
- Key properties: default_genome_build, default_lab, email_weekly_updates, columns, tag_colors, igv_port, timezone, allele_origin_focus

**GlobalSettings, OrganizationUserSettingsOverride, LabUserSettingsOverride, UserSettingsOverride** — Concrete settings levels

**TagColorsCollection** — Collection of color assignments for tags (FK to User)

**UserGridConfig** — Per-user grid configuration (rows, filters, visibility)

**UserPageAck** — Acknowledgment tracking (unique on user, page_id)

**UserDataPrefix** — Maps server paths to local paths (for BAM file access)

**AvatarDetails** — Avatar system integration with SpaceThemedAvatarProvider

### Custom Columns Models

**VariantGridColumn** — Built-in columns for data grids
- Fields: column_name, group_name, description, missing_value, source_field, custom_type

**CustomColumnsCollection** — User's selected columns for analysis display (FK to User)

**CustomColumn** — User-defined column combining built-in and custom logic

### Other Models

**Tag** — Simple tag system (id as primary key text)

**ClinVarKey** — API credentials for ClinVar submissions
- Fields: id (PK), api_key, org_id, name, last_full_run, assertion_method_lookup JSON, citations_mode

**CachedGeneratedFile** — UUID-based storage for generated files (graphs, downloads)
- Methods: get_or_create_and_launch(), check_ready()
- Pre-delete signal removes physical files

**SiteMessage** — System-wide messages displayed to all users (cached in Django cache)

**Wiki, ImportedWiki** — Markdown documentation system with InheritanceManager

**ClinGenAllele** — Bridges VariantGrid alleles to ClinGen Allele Registry

**GenomicIntervalsCollection / GenomicInterval** — Region-based analysis definitions

---

## Key Enums

| Enum | Values |
|------|--------|
| ImportStatus | CREATED, IMPORTING, REQUIRES_USER_INPUT, ERROR, SUCCESS, MARKED_FOR_DELETION, DELETING |
| VariantsType | UNKNOWN, GERMLINE, MIXED, SOMATIC_ONLY |
| SampleFileType | BAM, CRAM, BED, VCF |
| AlleleOrigin | IMPORTED_TO_DATABASE, IMPORTED_NORMALIZED, LIFTOVER, LIFTOVER_NORMALIZED |
| AlleleConversionTool | SAME_CONTIG, CLINGEN_ALLELE_REGISTRY, DBSNP, NCBI_REMAP, PICARD, CROSSMAP, BCFTOOLS_LIFTOVER |
| SequenceRole | ASSEMBLED_MOLECULE, UNLOCALIZED_SCAFFOLD, UNPLACED_SCAFFOLD, ALT_SCAFFOLD, FIX_PATCH, NOVEL_PATCH |
| CohortGenotypeCollectionType | UNCOMMON, COMMON |

---

## Key Views and URL Patterns

### Core Views (views/views.py)
- `index()` — Homepage
- `data()`, `vcfs()`, `samples()` — Data management pages
- `view_vcf(vcf_id)` — VCF details with statistics
- `view_sample(sample_id)` — Sample details with genotypes
- `cohorts()`, `view_cohort(cohort_id)` — Cohort management
- `trios()`, `view_trio(pk)` — Trio analysis
- `group_permissions()`, `bulk_group_permissions()` — Permission assignment
- `view_user_settings()` — User preference management
- `custom_columns()`, `tag_settings()`, `igv_integration()` — Custom settings
- `chrom_density_graph()`, `homozygosity_graph()` — Sample visualization
- `manual_variant_entry()` — Manual variant creation
- `liftover_runs()` — Manage liftover operations

### JSON Views (views/views_json.py)
- `job_status(job_id)` — Celery task status
- `create_cohort_genotype()` — Generate counts
- `vcf_populate_clingen_alleles()` — Batch ClinGen lookups
- `clone_custom_columns()`, `set_tag_color()`

### REST API (views/views_rest.py)
- `GET /api/sample_variant_zygosity/<sample_id>/<variant_id>/`
- `GET /api/trio/<pk>/`
- `GET /api/variant_allele_for_variant/<variant_id>/<genome_build_name>/`
- `POST /api/project/create`
- `GET /docs/` — OpenAPI 4.0.0 schema

### Autocomplete Views
- CohortAutocompleteView, SampleAutocompleteView, TrioAutocompleteView
- ProjectAutocompleteView, UserAutocompleteView, LabAutocompleteView
- VCFAutocompleteView, TagAutocompleteView

---

## Management Commands

| Command | Purpose |
|---------|---------|
| import_lab_info | Import lab metadata |
| replace_vcf | Replace VCF data in-place |
| create_zygosity_counts_for_existing_vcfs | Bulk recount zygosity |
| liftover_alleles | Run liftover between builds |
| delete_unused_variants | Garbage collection |
| clean_orphan_obj_perms_safe | Clean dangling permissions |
| user_emails | Bulk email operations |
| drf_user_api_token | DRF token management |
| site_messages | Manage system messages |
| expire_all_sessions | Force user re-login |

---

## Signals and Key Business Logic

### Signal Handlers

**user_post_save_handler()** — On new user creation:
- Adds user to PUBLIC_GROUP_NAME and LOGGED_IN_USERS_GROUP_NAME
- Creates UserDataPrefix if settings.INITIAL_USER_DATA_PREFIX_KWARGS set
- Auto-creates org/lab if settings.USER_CREATE_ORG_LABS configured
- Sends admin notification

**backend_vcf_import_success_handler()** — On VCF import success:
- Triggers BED file intersection calculations via Celery

**trio_post_save_handler()** — On trio creation:
- Creates genome karyomapping analysis via Celery task

**cgf_pre_delete_handler()** — On CachedGeneratedFile deletion:
- Removes physical file from disk

**vcf_pre_delete_handler()** — Cascades VCF deletion to cohort

### Key Data Flows

**VCF Import:**
1. VCF created (import_status=CREATED)
2. File parsed, samples created
3. Backend annotation triggered via Celery
4. On success: import_status=SUCCESS
5. CohortGenotypeCollection(s) created with packed genotype data
6. backend_vcf_import_success_signal fires → BED intersection calculated
7. Zygosity counts calculated

**Variant Normalization:**
1. Variant imported → Locus/Sequence/Variant created
2. ClinGen Allele Registry lookup (if settings allow)
3. If found → Allele created with clingen_allele link
4. VariantAllele bridges variant to build-independent allele
5. If liftover needed → AlleleLiftover created, run async

**Multi-sample Query:**
1. Annotate Variant queryset with CohortGenotypeCollection (get_annotation_kwargs)
2. Filter by packed zygosity string (get_zygosity_q)
3. Access sample-specific fields via array indices

---

## Performance Considerations

### Caching Strategies
- **Request-level caching**: Organization, Lab managers (ObjectManagerCachingRequest)
- **Immutable caching**: GenomeBuild (ObjectManagerCachingImmutable)
- **Django cache**: Site messages (30s TTL)
- **Timed cache**: GenomeBuild.get_name_or_alias() (60s)

### Data Partitioning
- **CohortGenotypeCollection**: Separate table per partition (variant count optimization)
- **Packed fields**: Samples stored as packed strings/arrays in single row (optimized for multi-sample queries, not normalized)

### Liftover Tools Supported
ClinGen Allele Registry (preferred), dbSNP API, NCBI Remap (obsolete), Picard LiftoverVCF, CrossMap, BCFtools/liftover

---

## Permission Model

- **Guardian**: Object-level permissions on VCF, Sample, Cohort, Trio
- **Group-based**: Organizations/Labs create Django groups for bulk permission management
- **Cascading**: Cohort inherits VCF permissions
- **User Settings**: Per-user, per-lab, per-org overrides for display preferences
