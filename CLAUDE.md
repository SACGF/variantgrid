# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

VariantGrid is a Django/PostgreSQL web application for storing, annotating, and classifying genomic variants. It supports multiple genome builds (GRCh37, GRCh38), integrates with Ensembl VEP for annotation, and manages ACMG-based variant classifications. Key deployments include Shariant (Australian variant sharing), SA Pathology clinical use, and variantgrid.com.

## Commands

### Running Tests

There are a lot of migrations in the project, so when running tests, you can save a lot of time
by using --keepdb ie:


```bash
# Run all tests
python3 manage.py test --keepdb

# Run a specific test module
python3 manage.py test --keepdb snpdb.tests.test_variant

# Run a specific test class or method
python3 manage.py test --keepdb snpdb.tests.test_variant.VariantTest.test_something
```

### Linting
```bash
# Run pylint (output to lint.txt)
./scripts/linting/run_pylint.sh

# Auto-fix whitespace/formatting issues with autopep8
./scripts/linting/format_code.sh
```

### Django Management
```bash
python3 manage.py runserver
python3 manage.py migrate
python3 manage.py shell
```

## Settings Architecture

Settings use a **hostname-based split-settings pattern**. `variantgrid/settings/__init__.py` auto-detects the hostname and loads the matching settings file:

1. First checks `variantgrid/settings/env_developers/<hostname>.py`
2. Falls back to `variantgrid/settings/env/<hostname>.py`

Each environment file imports from `variantgrid/settings/components/`:
- `default_settings.py` — all base Django settings
- `celery_settings.py` — RabbitMQ/Celery queues
- `annotation_settings.py` — VEP paths and annotation config
- `seqauto_settings.py` — sequencing automation

Secrets (DB credentials, API keys) are loaded from `/etc/variantgrid/settings_config.json` via `get_secret()`. See `config/settings_config.json` for a template.

The `DJANGO_SETTINGS_MODULE` defaults to `variantgrid.settings` which triggers hostname detection. Set it explicitly (e.g., `variantgrid.settings.env.vgtest`) for specific environments.

## App Architecture

There are per-app research documents generated in claude/research

### Core genetic apps (dependency order)
- **`library/`** — Shared utilities (not a Django app). No Django models; provides: permissions (`guardian_utils`), notifications (`log_utils.NotificationBuilder`), preview system, jQGrid/DataTables base classes, caching, and 20+ utility modules in `library/utils/`.
- **`snpdb/`** — Foundation genetic models: `Variant`, `Allele`, `Locus`, `GenomeBuild`, `VCF`, `Sample`, `Cohort`, `Lab`, `Organization`. Also called "SNPDB" (the original project name). Most other apps depend on this.
- **`genes/`** — Genes, transcripts, HGVS resolution, canonical transcript management. Contains `hgvs/` subdirectory with biocommons/pyhgvs converters.
- **`annotation/`** — VEP annotation integration, ClinVar, variant annotation versions.
- **`analysis/`** — Interactive DAG-based variant filtering pipeline. Nodes produce Django Q objects composed into querysets. Supports sample/trio/cohort/pedigree analysis modes.
- **`classification/`** — The largest app. Full ACMG classification workflow: evidence keys, versioned records, discordance detection, ClinVar export, condition text matching, multi-lab sharing.

### Supporting apps
- **`uicore/`** — Shared UI components: template tags, ValidatedJson, DataTables integration, Bootstrap/FontAwesome patterns.
- **`upload/`** — VCF upload and import processing pipeline.
- **`patients/`** — Patient and phenotype management.
- **`ontology/`** — HPO/OMIM/MONDO term management and matching.
- **`flags/`** — Flexible flagging system for any model.
- **`seqauto/`** — Sequencing automation (SeqAuto) workflows.
- **`sync/`** — Import/export syncing with external systems (Alissa, other VariantGrid instances).
- **`variantopedia/`** — Variant detail pages ("Variantopedia" wiki-style pages).
- **`vcauth/`** / **`oidc_auth/`** — Authentication/OIDC.

## Security

### Authentication
The project uses `global_login_required.GlobalLoginRequiredMiddleware`, which enforces login on **all** views globally. Individual views do **not** need `@login_required` decorators — their absence is intentional, not a security gap. Do not flag missing `@login_required` as a security issue during audits.

DRF is configured with `DEFAULT_PERMISSION_CLASSES = [IsAuthenticated]`, so all REST API endpoints require authentication by default. Individual API views do not need explicit `permission_classes` — their absence is intentional, not a security gap.

## Key Patterns

### Object-level permissions
All major models use Django Guardian for object-level permissions. The mixin `GuardianPermissionsMixin` (in `library/django_utils/guardian_permissions_mixin.py`) provides `can_view()`, `can_write()`, `filter_for_user()`. Standard groups are `all_users` and `public`. Use `assign_permission_to_user_and_groups()` from `library/guardian_utils.py`.

### Frontend
The project uses **Bootstrap 4**. Use `data-toggle` (not `data-bs-toggle`) and `data-target` (not `data-bs-target`) for collapse, modal, and other Bootstrap JS components.

### Grid/table views
Two systems coexist:
- **jQGrid** (legacy): `JqGridUserRowConfig` base class in `library/jqgrid/`. Still used in many places.
- **DataTables** (current): `DatatableConfig` + `RichColumn` in `snpdb/views/datatable_view.py`. Preferred for new tables.

### Celery task queues
Four worker queues: `analysis_workers`, `annotation_workers`, `db_workers` (default), `web_workers`, plus `scheduling_single_worker`. Assign tasks to appropriate queues via `@app.task(queue='...')`.

### Preview system
Models implement `PreviewModelMixin` to support hover-card previews. Apps connect to `preview_request_signal` and `preview_extra_signal` (in `library/preview_request.py`) to register their handlers. The `PreviewKeyValue` dataclass carries key/value pairs for the preview.

### Model readmes
Several apps have `__<app>_readme.md` files documenting architecture (e.g., `snpdb/__snpdb_readme.md`, `classification/__classification_readme.md`).

## Git Commits

Do NOT add "Co-Authored-By: Claude" or any similar co-author trailer to commit messages.

Reference GitHub issues in commit messages (e.g., `#1400`) but do NOT use keywords that auto-close issues (e.g., "fix", "close", "resolve"). Issues must go through a testing pipeline before being closed manually.

Before committing, check `git status` for already-staged changes unrelated to the current task. If any exist, stop and confirm with the user before proceeding — do not include them in the commit.

## GitHub Comments

When writing any comment on a GitHub issue or pull request, always preface it with 🤖 Written by Claude.

Do NOT close GitHub issues. Issues must go through a testing lifecycle before being closed by the user.

## Testing

Tests extend `django.test.TestCase`. URL tests use `URLTestCase` from `library/django_utils/unittest_utils.py`, which:
- Overrides settings for static files, disables Celery async, disables annotation caching
- Provides `_test_urls()` helper for batch URL status code testing

Fake/fixture data helpers are in `annotation/tests/test_data_fake_genes.py`, `snpdb/tests/utils/`, etc.

`UNIT_TEST = sys.argv[1:2] == ['test']` is set in default_settings and used to conditionally skip expensive setup.

## Database

PostgreSQL via `psqlextra` backend (`psqlextra.backend`), which adds PostgreSQL-specific features (partitioning, upserts). Redis is used for caching. `CACHE_VERSION` in settings must be incremented to flush caches after breaking changes.

## Classification App Notes

`Classification` records store evidence as JSON keyed by `EvidenceKey` slugs. Each edit creates a `ClassificationModification`. Only "published" modifications are visible outside the owning lab. The `ImportedAlleleInfo` model resolves HGVS → `Allele` during import. Discordance is auto-detected when classifications for the same allele span multiple clinical significance buckets (B/LB vs VUS vs LP/P).
