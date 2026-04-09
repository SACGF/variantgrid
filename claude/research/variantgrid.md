# VariantGrid: High-Level Architecture & Codebase Overview

## What Is VariantGrid?

VariantGrid is a Django-based web application for storing, analyzing, and classifying genetic variants. It is used by:
- Research institutions (public instance at VariantGrid.com)
- Australian clinical genomics via **Shariant** — a national variant classification sharing platform
- SA Pathology (South Australian public pathology service)
- Disease-specific databases (e.g., RUNX1Db for rare disease variant interpretation)

Its core function is to allow bioinformaticians and clinical geneticists to:
1. Upload VCF files and automatically annotate variants
2. Filter and analyze variant data across samples, cohorts, and trios
3. Classify variants using ACMG/AMP guidelines and custom evidence schemas
4. Share classifications across laboratories, detect discordances, and submit to ClinVar
5. Match variants to patient phenotypes and disease ontologies

**Technology stack**: Python 3, Django 5, PostgreSQL, Redis, Celery, Gunicorn

---

## Top-Level Directory Structure

```
variantgrid/
├── variantgrid/          # Main Django project (settings, urls, celery, wsgi)
├── snpdb/                # Core variant/sample database
├── classification/       # Variant classification system (largest app)
├── analysis/             # Interactive variant analysis pipeline
├── annotation/           # VEP integration and variant annotation
├── genes/                # Gene symbols, transcripts, gene lists, panels
├── patients/             # Patient demographics and sample linkage
├── pathtests/            # Pathology test definitions
├── pedigree/             # Pedigree file management
├── ontology/             # Disease/phenotype ontology (HPO, MONDO, OMIM)
├── upload/               # VCF and file import pipelines
├── seqauto/              # Sequencing run automation and QC
├── flags/                # Flagging system for variants/records
├── review/               # Review workflows
├── sync/                 # Data synchronization (Shariant, Alissa)
├── eventlog/             # Audit logging
├── email_manager/        # Email templating and delivery
├── expression/           # Placeholder (unused, migrations only)
├── variantopedia/        # Variant encyclopedia and metadata search
├── sapath/               # Symlink to SA Pathology custom app
├── oidc_auth/            # OpenID Connect authentication
├── vcauth/               # User authentication and authorization
├── uicore/               # UI utilities, templatetags, JSON views
├── library/              # Shared utilities (caching, permissions, graphs)
├── scripts/              # Deployment, startup, backup scripts
├── data/                 # Static reference data
├── config/               # Per-environment configuration overrides
├── manual/               # Manual deployment tasks
└── manage.py             # Django management entry point
```

---

## Django Apps Overview

| App | Purpose |
|-----|---------|
| **snpdb** | Core variant database — Variant, Allele, Sample, VCF, Project, Cohort, Trio, Lab, liftover, zygosity counts |
| **classification** | ACMG variant classification, discordance tracking, ClinVar export, condition matching |
| **analysis** | Interactive analysis pipeline with filterable graph nodes; drag-and-drop UI |
| **annotation** | VEP integration, variant annotation versions, citation management, ClinVar data |
| **genes** | HGNC gene symbols, transcripts (cdot), gene lists, Panel App integration, GenCC |
| **patients** | Patient demographics, phenotypes, external system identifiers |
| **pathtests** | Pathology test catalog |
| **pedigree** | Pedigree (PED) files and family member records |
| **ontology** | HPO, MONDO, OMIM, Orphanet ontology terms and relationships |
| **upload** | File upload pipeline — VCF import, variant insertion, annotation triggering |
| **seqauto** | Sequencing job scanning, QC tracking, run management |
| **flags** | Lightweight flagging of any domain object (variants, classifications, etc.) |
| **review** | Formal review workflows with questions and responses |
| **sync** | Bidirectional sync to Shariant and Alissa |
| **eventlog** | Audit trail — user actions, HTTP requests, severity levels |
| **email_manager** | Templated emails (HTML + plain text) |
| **variantopedia** | Variant detail pages, search, metadata aggregation |
| **oidc_auth** | OIDC/Keycloak SSO integration |
| **vcauth** | Fine-grained object-level permissions via django-guardian |
| **uicore** | Shared UI components, datatables, JSON views, template tags |
| **library** | Utilities: caching helpers, permission decorators, graphing, Django ORM utils |

---

## Core Data Model (snpdb)

The `snpdb` app is the foundation. All other apps reference its core entities:

- **GenomeBuild** — GRCh37, GRCh38, T2T-CHM13v2.0
- **Contig / GenomeBuildContig** — chromosome/contig mappings per build
- **Locus** — chromosome + position
- **Variant** — Locus + ref allele + alt allele (per genome build)
- **Allele** — Cross-build abstraction (one allele can have variants in multiple builds)
- **VariantAllele** — Junction: Allele ↔ Variant ↔ GenomeBuild
- **VCF** — Uploaded VCF file
- **Sample** — A sample within a VCF
- **Project / LabProject** — Organizational groupings
- **Lab / Organization / Company** — Hierarchy for labs and institutions
- **Cohort / CohortSample** — Groups of samples
- **Trio** — Parent–proband relationship for trio analysis
- **LiftoverRun / AlleleLiftover** — Cross-build variant coordinate conversion
- **UserSettings / UserGridConfig** — User preferences and column layouts

---

## Classification App (Most Feature-Rich)

The `classification` app is the largest and most complex. It handles:
- ACMG/AMP-based variant curation with flexible evidence keys
- Multi-lab classification sharing (share levels: User → Lab → Institution → Public)
- Discordance detection and resolution workflows
- ClinVar submission pipeline
- Condition/disease ontology matching
- Import/export of classification data in multiple formats

See `claude/research/classifications.md` for detailed coverage.

---

## Analysis Pipeline (analysis app)

The analysis app provides an interactive, graph-based variant filtering system:

- **Analysis** — A named analysis session belonging to a user/project
- **AnalysisTemplate** — Reusable analysis configurations
- **AnalysisNode** — A configurable filter/calculation step in a graph
  - Node types: AllVariants, CohortNode, ClassificationsNode, FilterNode, GeneListNode, VennNode, TagNode, TissueNode, etc.
- **NodeVersion** — Versioned node computation result
- **NodeCache** — Cached query results per node
- **VariantTag** — User-applied tags to variants within an analysis

Users build analyses visually by connecting nodes into a directed graph. Each node filters or computes on variants, with results cached and displayed in a JQGrid datatable.

---

## Variant Annotation Pipeline (annotation app)

Variants are annotated automatically using Ensembl VEP:

- **AnnotationVersion** — Tracks VEP version + database versions used
- **VariantAnnotation** — Per-transcript annotation (consequence, SIFT, PolyPhen, etc.)
- **CachedWebResource** — External resource cache (e.g., dbSNP, ClinVar data)
- **Citation** — Literature references

VEP is configured with numerous plugins:
- gnomAD v2.1.1, v4.0, v4.1 (population allele frequencies)
- SpliceAI (splice variant prediction)
- dbNSFP (in-silico predictions: REVEL, CADD, etc.)
- dbscSNV (splicing predictions)
- MaxEntScan
- PhastCons / PhyloP (conservation)
- COSMIC (somatic mutation database)

Genome builds supported: GRCh37, GRCh38, T2T-CHM13v2.0.

---

## Ontology App

Provides standardized disease/phenotype terms used in condition matching:

- **OntologyTerm** — HPO, MONDO, OMIM, Orphanet terms
- **OntologyRelation** — Parent/child/associated relationships
- **OntologyImport** — Import tracking for ontology files
- Fuzzy-matching logic to map free-text condition strings to standard terms

---

## Authentication & Authorization

- **Standard Django auth** — Username/password
- **OIDC / Keycloak** (`oidc_auth`) — SSO for institutional deployments
- **django-guardian** (`vcauth`) — Object-level permissions (per-analysis, per-sample access)
- **Share levels** on classifications (User / Lab / Institution / Public)
- **Lab membership** — Users belong to labs; labs belong to organizations

---

## Asynchronous Task Processing (Celery)

Celery handles long-running work. Key tasks include:

| Task | Frequency | Purpose |
|------|-----------|---------|
| `sync_data` | Hourly | Bidirectional sync to Shariant |
| `seqauto-scan` | 06:00 and 19:00 | Scan for new sequencing jobs |
| `discordance-emails-weekly` | Monday morning | Send discordance summary emails |
| `notify-server-status` | Daily 19:00 | System health notification |
| `warn-low-disk-space` | Hourly | Disk space monitoring |
| `heartbeat` | Every 30 min | Uptime monitoring ping |
| VEP annotation | On-demand | Batch variant annotation |
| VCF import | On-demand | Parse, normalize, insert variants |
| Liftover | On-demand | Cross-build coordinate conversion |
| Classification import | On-demand | Bulk classification processing |

---

## Settings Architecture

Settings are split by concern and loaded based on hostname:

```
variantgrid/settings/
├── __init__.py                    # Hostname-based settings loading
└── components/
    ├── default_settings.py        # Core Django settings, feature flags
    ├── annotation_settings.py     # VEP configuration, reference data paths
    ├── celery_settings.py         # Celery broker and beat schedule
    ├── seqauto_settings.py        # Sequencing automation config
    └── secret_settings.py         # Credentials from env/files
```

Environment-specific overrides live in `config/` (git-ignored) and `env_developers/`.

**Key feature flags in settings:**
- `UPLOAD_ENABLED` — Allow VCF uploads
- `SEQAUTO_ENABLED` — Sequencing automation
- `DISCORDANCE_ENABLED` — Inter-lab discordance tracking
- `CLINVAR_EXPORT_ENABLED` — ClinVar submission
- `USE_OIDC` — OIDC authentication
- `SHARIANT_SYNC_ENABLED` — Shariant data sync

---

## External Integrations

| Service | Purpose |
|---------|---------|
| **Ensembl VEP** | Variant annotation (v110) |
| **ClinGen Allele Registry** | Cross-variant allele identification |
| **ClinVar** | Submit classifications; import existing records |
| **Shariant** | Australian national variant classification sharing |
| **MONDO / HPO / OMIM / Orphanet** | Disease/phenotype ontologies |
| **Panel App Australia** | Gene-disease panel data |
| **GenCC** | Gene–condition relationship curation |
| **Keycloak** | OIDC SSO for institutions |
| **AWS S3** | File storage (if configured) |
| **AWS SES** | Email delivery (if configured) |
| **Rollbar** | Exception tracking |
| **Slack** | Admin notifications |
| **Heartbeat service** | Uptime monitoring |

---

## Key Dependencies

```
Django 5.x                  # Web framework
celery 5.6.0                # Async task queue
PostgreSQL + psycopg2       # Primary database
django-psqlextra            # Materialized views, JSON ops, partial indexes
redis                       # Cache + Celery broker
djangorestframework         # REST API
django-guardian             # Object-level permissions
django-crispy-forms         # Bootstrap form rendering
biopython 1.86              # Bioinformatics utilities
cdot 0.2.26                 # Canonical transcript library
hgvs                        # HGVS nomenclature parsing and normalization
cyvcf2                      # VCF file parsing
pandas / numpy / scipy      # Data analysis
reportlab                   # PDF generation
mozilla-django-oidc         # OIDC authentication
boto3 / django-storages     # AWS S3 integration
```

---

## URL Structure (Top-Level)

Each app registers its own URL namespace. The root router includes:

```
/                     → Index / redirect
/admin/               → Django admin
/accounts/            → Registration and login
/oidc/                → OIDC flows
/analysis/            → Analysis pipeline
/annotation/          → Annotation management
/classification/      → Classification system (100+ URLs)
/genes/               → Gene data management
/patients/            → Patient records
/snpdb/               → Variant database views
/upload/              → File upload
/variantopedia/       → Variant encyclopedia
/api/                 → REST API endpoints
```

---

## Management Commands

Custom management commands (via `manage.py`):

- `import_cdot_latest` — Update transcript annotations from cdot
- `import_hgnc` — Refresh gene symbol data from HGNC
- `import_canonical_transcript` — Update canonical transcript designations
- `ontology_import` — Import HPO/MONDO/OMIM data files
- `scan_run_jobs` — Detect new sequencing jobs (seqauto)
- `deployment_check` — System health validation
- `rematch_unmatched_gene_list_symbols` — Re-resolve gene symbols

---

## Deployment & Scripts

```
scripts/
├── install/              # Ubuntu/RHEL installation guides
├── backup/               # PostgreSQL backup scripts
├── dbscripts/            # Database utilities
├── start_services.sh     # Start all services
├── stop_services.sh      # Stop all services
├── restart_services.sh   # Restart services
├── celery_start.sh       # Start Celery workers
├── gunicorn_start.sh     # Start Gunicorn web server
└── upgrade.sh            # Version upgrade automation
```

Production runs: Gunicorn + Nginx (external) + Celery workers + Redis + PostgreSQL.
