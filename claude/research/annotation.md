# Annotation App - VariantGrid Research Document

## Overview

The `annotation` Django app manages variant annotations from VEP (Variant Effect Predictor) and other external sources. It handles the full annotation pipeline: VEP execution, result storage in partitioned PostgreSQL tables, ClinVar integration, gene-level annotations, protein data, and citation management. It is a core app that most analysis features depend on.

---

## Versioning Architecture

The annotation system uses a hierarchical versioning model where a top-level `AnnotationVersion` groups together sub-versions for different annotation data sources.

### AnnotationVersion

Top-level grouping of all annotation sub-versions.

- **Fields:**
  - `variant_annotation_version` (FK to VariantAnnotationVersion)
  - `gene_annotation_version` (FK to GeneAnnotationVersion, nullable)
  - `clinvar_version` (FK to ClinVarVersion, nullable)
  - `human_protein_atlas_version` (FK to HumanProteinAtlasAnnotationVersion)
  - `ontology_version` (FK)
- **Key Methods:**
  - `latest(genome_build)`: Returns the most recent AnnotationVersion for a build
  - `validate()`: Checks sub-version consistency - raises `InvalidAnnotationVersionError` if inconsistent (checks that all sub-versions are present, that the same GeneAnnotationRelease is referenced across variant and gene versions, and that OntologyVersion matches)
  - `new_sub_version(genome_build)`: Creates a new sub-version for incremental updates
  - `sql_partition_transformer()`: Modifies generated SQL to query the correct PostgreSQL partition

`SubVersionPartition` is the abstract base class for all sub-version types.

---

### VariantAnnotationVersion

VEP-specific version record.

- **Fields:**
  - `genome_build`
  - `annotation_consortium`: `REFSEQ` or `ENSEMBL`
  - `gene_annotation_release` (FK, nullable)
  - `vep` (int): VEP version number
  - `vep_cache`: VEP cache identifier string
  - `columns_version`: `1`, `2`, or `3` - tracks schema evolution of annotation columns
  - Ensembl/RefSeq version numbers
  - Data source versions: `thousand_genomes`, `cosmic`, `hgmd`, `dbsnp`, `gencode`, `sift`, `dbnsfp`, `gnomad`, `refseq`
  - `distance` (default 5000): Distance in bases for upstream/downstream gene overlap
  - `active` (bool)
- **Key Methods:**
  - `latest(genome_build, active=True)`: Returns the most recent active version for a build
  - `get_pathogenic_prediction_funcs()`: Returns callables for pathogenicity prediction scoring
  - `damage_predictions_description`: Human-readable description of active damage predictors
  - `cdot_gene_release_filename`: Filename for the cdot gene release associated with this version

---

## Core Annotation Models

### AbstractVariantAnnotation

Abstract base class shared by `VariantAnnotation` and `VariantTranscriptAnnotation`. Contains the fields present for every transcript annotation.

- **Fields:**
  - `version` (FK to VariantAnnotationVersion)
  - `variant` (FK)
  - `annotation_run` (FK to AnnotationRun)
  - `gene`, `transcript`, `transcript_version` (FKs, nullable)
  - `amino_acids`
  - `cadd_phred`
  - `canonical` (bool): Whether this is the canonical transcript for the gene
  - `codons`
  - `consequence`: VEP consequence term(s)
  - `distance`: Distance to nearest gene (for intergenic variants)
  - `domains`
  - `ensembl_protein`
  - `exon`: Exon number/total (e.g., `"3/10"`)
  - `intron`: Intron number/total
  - `fathmm_pred_most_damaging`
  - `flags`
  - `gerp_pp_rs`
  - `grantham`: Grantham score for amino acid substitution
  - `hgvs_c`: HGVS coding notation
  - `hgvs_p`: HGVS protein notation
  - `impact`: VEP impact (HIGH/MODERATE/LOW/MODIFIER)
  - `interpro_domain`
  - `maxentscan_alt`, `maxentscan_diff`, `maxentscan_ref`: MaxEntScan splice scores
  - `mutation_assessor_pred`
  - `polyphen2_hvar_pred`
  - `protein_position`
  - `revel_score`
  - `sift`
  - `splice_region`
  - `symbol`: Gene symbol from VEP
  - `mavedb_score`, `mavedb_urn`: MaveDB functional score data

---

### VariantAnnotation

The representative (canonical) transcript annotation for a variant. One record per variant per version. Extends `AbstractVariantAnnotation` with population frequency, clinical, and genomic fields.

- **Additional Fields:**
  - `hgvs_g`: HGVS genomic notation
  - **Population Frequencies:**
    - `af_1kg`: 1000 Genomes allele frequency
    - `af_uk10k`: UK10K allele frequency
    - `topmed_af`: TOPMed allele frequency
    - `gnomad_af`: gnomAD overall allele frequency
    - `gnomad2_liftover_af`: gnomAD v2 liftover frequency
    - `gnomad_ac`, `gnomad_an`, `gnomad_hom_alt`: gnomAD allele count/number/homozygotes
  - **gnomAD by Population:**
    - `gnomad_afr_af`, `gnomad_amr_af`, `gnomad_asj_af`, `gnomad_eas_af`, `gnomad_fin_af`, `gnomad_mid_af`, `gnomad_nfe_af`, `gnomad_oth_af`, `gnomad_sas_af`
  - **gnomAD v4 Filtering Allele Frequencies:**
    - `gnomad_faf95`, `gnomad_faf99`
    - `gnomad_fafmax_*` fields
  - **gnomAD Sex Chromosome:**
    - `gnomad_xy_af`, `gnomad_xy_ac`, `gnomad_xy_an`
    - `gnomad_hemi_count`
    - `gnomad_non_par`: Whether variant is outside pseudoautosomal region
  - **gnomAD Popmax:**
    - `gnomad_popmax_af`, `gnomad_popmax_ac`, `gnomad_popmax_an`, `gnomad_popmax_hom_alt`, `gnomad_popmax`: population code
  - **gnomAD Structural Variants:**
    - `gnomad_sv_overlap_af`, `gnomad_sv_overlap_percent`, `gnomad_sv_overlap_name`
  - **Database IDs:**
    - `dbsnp_rs_id`
    - `cosmic_id` (COSV format), `cosmic_count`
    - `mastermind_count_1`, `mastermind_count_2`, `mastermind_count_3`
  - **Splice Predictors:**
    - `dbscsnv_ada_score`, `dbscsnv_rf_score`
    - `spliceai_pred_dp_ag`, `spliceai_pred_dp_al`, `spliceai_pred_dp_dg`, `spliceai_pred_dp_dl`
    - `spliceai_pred_ds_ag`, `spliceai_pred_ds_al`, `spliceai_pred_ds_dg`, `spliceai_pred_ds_dl`
    - `spliceai_gene_symbol`
  - **Rank Scores:**
    - `cadd_raw_rankscore`, `revel_rankscore`, `bayesdel_noaf_rankscore`, `clinpred_rankscore`, `vest4_rankscore`, `metalr_rankscore`, `alphamissense_rankscore`
  - **Aloft:**
    - `aloft_prob_tolerant`, `aloft_prob_recessive`, `aloft_prob_dominant`, `aloft_pred`
  - **Conservation Scores:**
    - `phylop_30way_mammalian`, `phylop_46way_primate`, `phylop_100way_vertebrate`
    - `phastcons_30way_mammalian`, `phastcons_46way_primate`, `phastcons_100way_vertebrate`
  - **Prediction Summaries:**
    - `predictions_num_pathogenic`, `predictions_num_benign`
  - **Other:**
    - `somatic` (bool)
    - `variant_class`
    - `overlapping_symbols`: Pipe-separated list of overlapping gene symbols
    - `pubmed`: PubMed IDs from VEP
- **Properties:** `gnomad_url`, `has_gnomad`, `has_cosmic`, `has_splicing`, `has_conservation`, and other boolean presence checks

---

### VariantTranscriptAnnotation

Stores annotation for every transcript overlapping a variant. One record per transcript per variant per version. Partitioned by VariantAnnotationVersion.

---

### VariantGeneOverlap

Stores gene overlap records for variants that are upstream or downstream of a gene. Partitioned by VariantAnnotationVersion. Used for distance-based gene filtering.

---

## ClinVar Models

### ClinVarVersion

Represents a ClinVar VCF import.

- **Fields:**
  - `filename`
  - `sha256_hash`: Deduplication of imports
  - `genome_build`
- Pre-delete signal cascades deletion to dependent ClinVar records

---

### ClinVar

Clinical significance classifications per variant.

- **Fields:**
  - `variant` (FK, on_delete=PROTECT)
  - `version` (FK to ClinVarVersion)
  - `clinvar_variation_id`
  - `clinvar_allele_id`
- **Germline fields:**
  - `preferred_disease_name`
  - `disease_database_name`
  - `review_status`
  - `clinical_significance`
  - `conflicting_clinical_significance`
  - `highest_pathogenicity` (int): Numeric pathogenicity rank
  - `clinical_sources` (pipe-separated)
  - `origin` (bitfield): Allele origin flags
  - `suspect_reason_code`
- **Oncogenic fields:**
  - `oncogenic_review_status`
  - `oncogenic_classification`
- **Somatic fields:**
  - `somatic_review_status`
  - `somatic_clinical_significance`
- **Properties:**
  - `stars`: Review status star rating (0-4)
  - `is_expert_panel_or_greater`: Whether reviewed by an expert panel
  - `allele_origin_bucket`: Simplified origin category
  - `germline_disease_database_terms`: Parsed disease database terms
  - `citation_ids`: Associated citation IDs
  - `json_summary()`: JSON representation of key fields

---

### ClinVarRecordCollection

Caches individual ClinVar record data retrieved from ClinVar XML (per variation ID, not per VCF import).

- **Fields:**
  - `clinvar_variation_id` (PK)
  - `allele` (FK, nullable)
  - `urls` (ArrayField): Source XML URLs
  - `last_loaded`: Timestamp of last fetch
  - `parser_version`: Tracks parser version for cache invalidation
  - `max_stars`: Highest star rating among contained records
  - `expert_panel` (OneToOne FK): Expert panel record if present
- **Key Methods:**
  - `update_with_records_and_save()`: Updates collection from fresh XML data

---

### ClinVarRecord

An individual ClinVar submission (SCV accession).

- **Fields:**
  - `record_id` (PK, SCV accession)
  - `clinvar_record_collection` (FK)
  - `org_id`
  - `submitter`
  - `genome_build`
  - `review_status`
  - `stars` (0-4)
  - `submitter_date`
  - `date_last_evaluated`
  - `c_hgvs`: HGVS coding notation
  - `variant_coordinate`
  - `condition`
  - `clinical_significance`
  - `somatic_clinical_significance`
  - `gene_symbol`
  - `interpretation_summary`
  - `allele_origin_bucket`

---

## Gene-Level Annotation Models

### GeneAnnotationVersion

Sub-version for gene-level annotations. Extends `SubVersionPartition`.

- **Fields:**
  - `gene_annotation_release` (FK)
  - `ontology_version` (FK, nullable)
  - `dbnsfp_gene_version` (FK, nullable)

---

### GeneAnnotation

Per-gene annotation record.

- **Fields:**
  - `gene` (FK)
  - `version` (FK to GeneAnnotationVersion)
  - `dbnsfp_gene` (FK, nullable)
  - `hpo_terms`: Pipe-separated HPO term IDs
  - `omim_terms`: Pipe-separated OMIM term IDs
  - `mondo_terms`: Pipe-separated MONDO term IDs
  - `gene_disease_moderate_or_above`: Bool flag for moderate+ gene-disease associations
  - `gene_disease_supportive_or_below`: Bool flag for supportive/below gene-disease associations
  - `gnomad_oe_lof`: gnomAD observed/expected LoF ratio

---

### DBNSFPGeneAnnotation

Gene-level scores from the dbNSFP database. Partitioned by version.

- **Fields:**
  - `gene_symbol`
  - `refseq_transcript`, `ensembl_transcript` (FKs)
  - **Constraint/Essentiality scores:**
    - `gene_damage_index_score`, `gene_damage_index_phred`, `gene_damage_index_phi`
    - `ghis`, `prec`, `hipred_score`
    - `gnomad_pli`, `gnomad_prec`, `gnomad_pnull`
    - `loftool`
    - `gene_indispensability_score`, `gene_indispensability_pred`
  - **Pathway:**
    - `pathway_biocarta`
    - `pathway_kegg`
  - **Disease/Trait:**
    - `gwas_trait_association`
  - **Gene Ontology:**
    - `gene_ontology_biological_process`
    - `gene_ontology_cellular_component`
    - `gene_ontology_molecular_function`
  - **Interactions:**
    - `interactions_biogrid`
  - **Expression:**
    - `expression_egenetics`
    - `expression_gnf_atlas`
  - **Essential Gene Screens:**
    - `essential_gene_crispr`
    - `essential_gene_crispr2`
    - `essential_gene_gene_trap`

---

## Protein and Tissue Data

### HumanProteinAtlasAnnotationVersion

Version record for Human Protein Atlas imports.

---

### HumanProteinAtlasTissueSample

Reference tissue/cell type records from the Human Protein Atlas.

---

### HumanProteinAtlasAnnotation

Tissue-level protein abundance data per gene.

- **Fields:**
  - `gene_symbol` (FK)
  - `tissue_sample` (FK to HumanProteinAtlasTissueSample)
  - `value` (float): Abundance value

---

## VEP Field Management

### ColumnVEPField

Maps VEP output fields to `VariantGridColumn` database columns. Controls which VEP fields are available for which pipeline versions.

- **Fields:**
  - `column` (unique): The database column name
  - `variant_grid_column` (FK to VariantGridColumn)
  - `genome_build` (nullable = applies to all builds)
  - `pipeline_type` (nullable = applies to all pipelines)
  - `category`: Field category for grouping
  - `source_field`: Raw VEP output field name
  - `vep_plugin`, `vep_custom`: Whether sourced from a VEP plugin or custom annotation
  - `min_columns_version`, `max_columns_version`: Schema version range
  - `min_vep_version`, `max_vep_version`: VEP version range
  - `summary_stats`: Whether to compute summary statistics
- **Key Methods:**
  - `filter_for_build()`: Returns fields applicable to a given genome build
  - `get_source_fields()`: Returns the VEP source field names to extract

---

## Citations

### Citation

Literature reference record.

- **Fields:**
  - `id` (PK): Formatted as `"PMID:1234567"`, `"PMCID:..."`, etc.
  - `title`
  - `abstract`
  - `source`: Source database

---

### CitationFetchRequest / CitationFetchResponse

Async fetch infrastructure for retrieving citation data from external databases. Requests are queued and responses cached.

---

### ClinVarCitationsCollection / ClinVarCitation

Links ClinVar records to Citation records. Enables display of supporting literature alongside ClinVar classifications.

---

## Supporting Models

### VCFAnnotationStats

Per-VCF, per-annotation-version statistics about annotated variants. Used for QC reporting.

---

### AnnotationRun

Tracks annotation processing jobs.

- **Fields:**
  - `status`: Job status
  - `pipeline_type`: `STANDARD` or `CNV`
  - `annotation_version_id`
  - `dump_count`: Number of variants dumped for annotation

---

### ManualVariantEntry

Stores user-entered variant assertions outside the normal VCF import pipeline.

---

### GenePubMedCount

Caches publication count per gene. Used for prioritizing genes by literature evidence.

---

### CachedWebResource

Dynamic external data loader. Loads and caches data from external sources including:
- HGNC gene nomenclature
- MANE transcript designations
- Pfam protein domains
- Panel App gene panels
- gnomAD gene constraint scores

---

## VEP Integration Flow

1. Variants are collected and written to a temporary VCF
2. VEP is configured per genome build (consortium, cache, plugins, custom annotations)
3. VEP executes and produces annotated VCF output
4. `BulkVEPVCFAnnotationInserter` parses VEP output and stores results in partitioned tables:
   - `VariantAnnotation`: One record per variant (representative/canonical transcript)
   - `VariantTranscriptAnnotation`: One record per variant per transcript
5. `ColumnVEPField` mapping controls which VEP output columns are written to which database columns

---

## Partition Strategy

Each `VariantAnnotationVersion` gets its own PostgreSQL table partition for `VariantAnnotation`, `VariantTranscriptAnnotation`, and `VariantGeneOverlap`.

- `sql_partition_transformer()` on `AnnotationVersion` modifies ORM-generated SQL to target the correct partition at query time
- `ColumnVEPField` tracks field availability per version via `min/max_columns_version` and `min/max_vep_version`
- This design allows old annotation data to be retained while new versions are added without schema migration

---

## AnnotationVersion Consistency Validation

`AnnotationVersion.validate()` enforces consistency across sub-versions:

1. All required sub-versions must be present (non-null)
2. The `GeneAnnotationRelease` referenced by `VariantAnnotationVersion` must match the one referenced by `GeneAnnotationVersion`
3. The `OntologyVersion` must be consistent

Raises `InvalidAnnotationVersionError` with details if any check fails.

---

## Admin Interfaces

| Admin Class | Notable Features |
|---|---|
| ClinVarRecordCollection admin | `refresh_old()` action to re-fetch stale records; `refresh_force_vcv()` to force re-fetch by VCV ID |
| ClinVar admin | Cache filter: `Fresh` vs `Stale` to identify records needing refresh |
| Citation admin | Searchable by `id` and `source` fields |

---

## Celery Tasks

| Task | Purpose |
|---|---|
| `annotate_variants` | Drives VEP execution and result insertion |
| `annotation_scheduler_task` | Schedules pending annotation jobs |
| `cached_web_resource_tasks` | Updates HGNC, MANE, Pfam, Panel App cached data |
| `calculate_sample_stats` | Computes per-sample annotation statistics |
| `cohort_sample_gene_damage_counts` | Aggregates gene damage counts across cohort samples |
| `import_clinvar_vcf_task` | Imports a ClinVar VCF file |
| `process_manual_variants_task` | Processes manually entered variants |

---

## Data Flow Summary

1. **VEP annotation:** Variants → VEP execution → `BulkVEPVCFAnnotationInserter` → partitioned `VariantAnnotation` / `VariantTranscriptAnnotation` tables
2. **ClinVar:** ClinVar VCF import → `ClinVar` records per variant; XML fetch → `ClinVarRecordCollection` + `ClinVarRecord` for individual submissions
3. **Gene annotation:** `GeneAnnotationVersion` bundles gene-level annotations from dbNSFP, ontology terms, and disease associations into `GeneAnnotation`
4. **Protein data:** Human Protein Atlas imports populate `HumanProteinAtlasAnnotation` for tissue expression data
5. **Version control:** All annotation data is scoped to versioned objects; `AnnotationVersion.validate()` ensures cross-source consistency before use
