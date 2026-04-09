# VariantGrid Variantopedia App — Reference Document

## Purpose and Overview

The `variantopedia` app is the variant encyclopedia and metadata search module. It provides the primary interface for searching, viewing, and managing variant data, plus server health monitoring.

**Key Design:** No new models defined; orchestrates data from snpdb, annotation, classification, genes, and analysis apps into comprehensive variant views.

---

## Key Views and URL Patterns

### Core Variant Views
- `variants(request, genome_build_name=None)` — `/variants` — Entry point for variant browsing
- `view_variant(request, variant_id, genome_build_name=None)` — `/view_variant/<variant_id>` — Main variant detail page; handles genome build selection; includes IGV data
- `variant_details_annotation_version(request, variant_id, annotation_version_id, template, extra_context)` — `/view/<variant_id>/<annotation_version_id>` — Core variant details with specific annotation version; generates transcript selections, canonical transcripts, annotation descriptions
- `view_variant_annotation_history(request, variant_id)` — `/view_variant_annotation_history/<variant_id>` — Historical annotation versions

### Search and Discovery
- `search(request)` — `/search` — Multi-type search; supports SearchForm and SearchAndClassifyForm; auto-redirects to single preferred result

### Allele Views
- `view_allele(request, allele_id)` — `/view_allele/<allele_id>` or `/a<allele_id>` — Allele-level data across all genome builds; AlleleOriginGrouping descriptions with overlap status; classification filtering per user; ClinGen Allele Registry link
- `view_allele_from_variant(request, variant_id)` — Redirects variant to linked allele (respects PREFER_ALLELE_LINKS setting)
- `export_classifications_allele(request, allele_id)` — CSV export of classifications for an allele
- `create_variant_for_allele(request, allele_id, genome_build_name)` — POST; creates liftover pipelines across genome builds

### Variant Tagging
- `variant_tags(request, genome_build_name=None)` — `/variant_tags/` — Lists all variant tags with counts
- `variant_tag_detail(request, variant_id, tag)` — Details for a specific tag on a variant
- `variant_tags_export(request, genome_build_name)` — CSV export of tags grid
- `tagged_variant_export(request, genome_build_name)` — CSV export of tagged variants grid

### Server Administration (superuser only)
- `server_status(request)` — `/server_status` — Comprehensive health monitoring: Celery workers, long-running SQL queries, reference FASTA accessibility, disk usage, variant annotation status; actions: Test Slack, Health Check, Test Rollbar, Kill PID
- `server_status_activity(request, days_ago)` — `/server_status_activity/detail/<int:days_ago>` — Historical activity with dashboard notices
- `server_status_settings(request)` — HGVS matcher status, import run overview
- `health_check_details(request)` — `/health_check_details` — Detailed health check results from all registered signals
- `database_statistics(request)` — `/database_statistics/detail` — Variant upload statistics; cumulative samples/variants/genotypes

### Analysis Views
- `nearby_variants_tab(request, variant_id, annotation_version_id)` — Nearby variants (tab format)
- `nearby_variants(request, variant_id, annotation_version_id)` — Full nearby variants page
- `variant_sample_information(request, variant_id, genome_build_name)` — Sample-level info; multiallelic handling
- `gene_coverage(request, gene_symbol_id)` — Gene-level 20x coverage statistics
- `variant_wiki(request, genome_build_name=None)` — Variant wiki with markdown editing
- `dashboard(request)` — User-facing dashboard; latest sequencing VCFs, sample enrichment, case existence check

---

## Business Logic

### Nearby Variant Discovery
Multiple search strategies based on variant location:
- **Codon**: Variants in same codon
- **Exon**: Variants in same exon
- **Domain**: Variants in overlapping protein domains
- **Gene**: Variants in overlapping gene(s)
- **Range**: Variants within configurable distance (e.g., 50bp)
Each returns filtered queryset plus summary text; supports clinical significance filtering.

### Allele Origin Grouping Analysis (AlleleOriginGroupingDescription)
Categorizes classifications by AlleleOriginBucket (germline, somatic, etc.) and analyzes overlap status:
- AGREEMENT, CONFIDENCE, DISCORDANCE, MEDICALLY_SIGNIFICANT_DISCORDANCE
- SINGLE_SUBMITTER, NO_SHARED_RECORDS, NOT_COMPARABLE

### Server Health Monitoring
- Celery worker inspection via Celery control API
- Long-running query detection from PostgreSQL pg_stat_activity (`RunningQuery` dataclass)
- Reference genome accessibility verification
- Disk space monitoring with thresholds

---

## Grids

- **AllVariantsGrid** — All variants in a genome build; customizable columns; sortable by zygosity counts
- **NearbyVariantsGrid** — Variants near a target; filter by region type (codon, exon, domain, range, genes)
- **VariantTagsGrid** — Tag-centric (one row per tag-variant pair); filter by analysis/gene/tag
- **TaggedVariantGrid** — Variant-centric (variants that have tags); filter by specific tag
- **VariantTagCountsColumns** — Tag counts per variant; expandable detail rows
- **VariantWikiColumns** — Variant wiki entries; filter by genome build

---

## Forms
- **SearchForm** — search (CharField), mode (optional "preview")
- **SearchAndClassifyForm** — extends SearchForm; adds classify (BooleanField, hidden)

---

## Template Tags

**nearby_variants_tags.py:**
- `@nearby_variants(context, variant, annotation_version, clinical_significance=True)` — inclusion_tag; displays nearby variants within variant details context

---

## Key Data Structures

**AlleleOriginGroupingDescription** (dataclass):
- Properties: allele_origin_grouping, discordance_report, overlap_status, shared_counts, unshared_counts, should_show_diffs
- Static method: `describe()` — performs overlap analysis

**ShareLevelRecordCounts** (dataclass): Classification share level statistics with lab_count.

---

## Integration Points

| App | Integration |
|-----|-------------|
| snpdb | Variant, Allele, Sample, VCF, GenomeBuild models; variant search; IGV data; Tag/UserSettings |
| annotation | VariantAnnotation, AnnotationVersion, transcript annotations, canonical transcripts |
| classification | Classification, AlleleOriginGrouping, DiscordanceReport, ClassificationGrouping, ClinicalSignificance |
| genes | GeneSymbol, CanonicalTranscript, transcript selection logic |
| analysis | VariantTag model for user-created variant tags |
| seqauto | VCFFromSequencingRun, sample enrichment kits, 20x gene coverage stats |
| eventlog | Event logging for search actions |
| pathtests | Cases for user functionality |
| patients | Clinician model for role-based access |
| library | Health check signals, health_check_signal aggregation |
