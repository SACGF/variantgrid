# Circular Import Reduction Plan

How the import graph looks today, which cycles are worth breaking, and how to stop new ones
appearing. Ordered so that each stage is independently shippable.

## How the measurements were taken

Three tools, each answering a different question:

```bash
# 1. pylint's own cyclic-import checker (already part of our lint stack)
DJANGO_SETTINGS_MODULE=variantgrid.settings pylint \
  --load-plugins pylint_django --rcfile=config/pylint3.rc \
  --disable=all --enable=cyclic-import --ignore-paths='.*/migrations/.*' \
  analysis annotation classification eventlog flags genes library manual mme ontology \
  pathtests patients pedigree seqauto snpdb sync uicore upload variantgrid variantopedia \
  review email_manager user_messages

# 2. grimp — full import graph, strongly connected components
uv pip install import-linter          # brings in grimp
```

`grimp` builds the graph (1886 modules) and `networkx.strongly_connected_components` gives the
cycles. A third pass with `ast` distinguishes **top-level** imports (the ones that actually run at
`import` time and can raise `ImportError`) from **function-level** imports (the workarounds already
in the tree). Migrations and tests are excluded throughout — tests are leaves and never imported by
production code, so a test importing another app does not create a runtime cycle.

pylint reports 357 `R0401` chains, but most are the same handful of cycles re-reported through
different paths, so the SCC view below is the one to work from.

## Current state

**App level:** all 23 first-party apps sit in a single strongly connected component. 47 app pairs
import each other. This is largely inherent — Django apps with cross-app ForeignKeys are mutually
dependent by construction, and no amount of refactoring changes that `snpdb.Variant` and
`genes.Transcript` reference each other. The app-level graph is *not* the thing to fix.

**Module level (top-level imports only) — 6 genuine cycles, 28 modules:**

| Size | Cycle | Cause |
|---|---|---|
| 9 | `classification.models` ↔ 8 of its own submodules | submodules import the package `__init__` aggregator |
| 8 | `sync.models` ↔ `sync.alissa` ↔ `sync.shariant` ↔ `sync.sync_runner` | registration-by-import via `__init__` |
| 4 | `classification.views.classification_export_view` ↔ `classification.views.exports.*` | aggregator |
| 3 | `library.utils` ↔ `diff_utils`, `os_utils` | submodules import the star-export package |
| 2 | `analysis.models` ↔ `analysis.models.models_candidate_search` | aggregator |
| 2 | `classification.signals` ↔ `classification_hooks_import_notifications` | aggregator |

**Module level (counting the inline-import workarounds too) — 9 cycles, 48 modules.** The three
extra ones are the genuinely hard, model-level knots:

- 12 modules: `genes.hgvs` ↔ `genes.models` ↔ `annotation.models.models` ↔ `snpdb.models.models_variant` ↔ `snpdb.clingen_allele`
- 7 modules: `analysis.models.nodes.sources` ↔ `_stats_cache` ↔ the source node classes
- 4 modules: `patients.models` ↔ `annotation.models.models_phenotype_match` / `phenotype_matching`
- 3 modules: `classification.models.classification` ↔ `discordance_models` ↔ `clinical_context_models`
- 2 modules: `snpdb.models.models` ↔ `snpdb.models.models_user_settings`

**Workaround debt:** 183 function-level first-party imports. Excluding the legitimate ones
(`apps.py` `ready()` signal registration, ~22), that is roughly 160 imports that exist to dodge a
cycle. Concentrated in `snpdb` (49), `classification` (37), `analysis` (19), `annotation` (16),
`ontology` (16), `genes` (13). CLAUDE.md's rule is that inline imports are only for genuine
circular deps — so this number is a fair proxy for "cycle debt", and it should go down.

---

## Stage 1 — Move framework code out of `snpdb` into `library`

Highest value for the effort: three modules in `snpdb` are generic Django infrastructure with **no
snpdb dependencies at all**, yet 16 apps import them, which is what drags `snpdb` into the middle
of everything.

| Module | Depends on | Imported by |
|---|---|---|
| `snpdb/admin_utils.py` | `library.log_utils`, `library.utils`, guardian, django only | 24 modules across 17 apps |
| `snpdb/views/datatable_view.py` | `library.*` + `datatable_mixins` only | 34 modules across 14 apps |
| `snpdb/views/datatable_mixins.py` | django only | via `datatable_view` |

Move them to:

- `library/django_utils/admin_utils.py`
- `library/django_utils/datatables/datatable_view.py`
- `library/django_utils/datatables/datatable_mixins.py`

and update the ~58 import sites (mechanical `sed`, then run the test suite).

This on its own fully breaks `email_manager ↔ snpdb` and `snpdb ↔ user_messages`, and removes the
only `snpdb` dependency from `manual`, `mme`, `review`, and `sync`'s admin layers. It also stops
`library` looking like a leaf while every app has to reach into `snpdb` for a `ModelAdmin` base
class.

CLAUDE.md's "Grid/table views" section names `snpdb/views/datatable_view.py` as the home of
`DatatableConfig` / `RichColumn` — update that reference as part of the move.

## Stage 2 — Make `library` a true base layer

`library` is documented as "shared utilities (not a Django app)", but six top-level imports reach
back up into apps:

| Import | Fix |
|---|---|
| `library.keycloak` → `snpdb.models.models`, `snpdb.models.models_user_settings` | move `library/keycloak.py` → `snpdb/keycloak.py` (3 importers) |
| `library.jqgrid.jqgrid_user_row_config` → `snpdb.models.UserGridConfig` | move → `snpdb/jqgrid_user_row_config.py` (10 importers) |
| `library.genomics.vcf_utils` → `snpdb.models` | see below |
| `library.log_utils` → `eventlog.models.Event` | see below |
| `library.preview_request` → `variantgrid.perm_path` | leave; see Stage 5 |

Plus four inline workarounds in `library` that disappear once the same rule is applied:
`library/utils/export_utils.py` → `snpdb.user_settings_manager`, `library/guardian_utils.py` →
`snpdb.models.UserSettings`, `library/django_utils/django_partition.py` →
`snpdb.models.models_partition_archive`, `library/preview_request.py` → `genes.models_enums`.

**`vcf_utils`** (20 importers, so moving it is the expensive option). It uses `VariantCoordinate`,
`Variant.REFERENCE_ALT`, `Sequence.allele_is_symbolic`, `GenomeFasta`, `SequenceRole`.
`VariantCoordinate` is a `pydantic.BaseModel` value type living in `snpdb/models/models_variant.py`
— it is not a Django model and is the natural thing to hoist into `library/genomics/`. Doing that
plus taking `genome_fasta` as a parameter rather than looking it up removes most of the
dependency; what remains (`Variant.REFERENCE_ALT`, `SequenceRole`) is small enough to pass in or
duplicate as a genomics-level constant. If that turns out to be more surgery than it is worth,
moving the whole module to `snpdb/vcf_utils.py` is the fallback.

**`log_utils` → `eventlog.Event`**: `report_event()` is the only user. Either give `library` a
tiny sink registry that `eventlog.apps.ready()` populates, or move `report_event` itself into
`eventlog` (it is an event-log function living in the wrong app). The second is simpler and is the
recommendation.

`library/utils/export_utils.py`'s inline `UserSettingsManager` import is for the user's timezone —
pass the tzinfo in from the caller instead.

## Stage 3 — Break the package-aggregator self-cycles

Every one of the six top-level cycles has the same shape: `app/models/__init__.py` re-exports the
submodules, and a submodule then imports `from app.models import X` instead of naming the module
that defines `X`. Python survives it only because of `__init__` execution order, and it is the
single most common reason a new import here explodes.

The fix is purely mechanical — import the defining module:

| Cycle | Change |
|---|---|
| `classification.models` (9 modules) | 8 submodules: `from classification.models import X` → `from classification.models.<defining_module> import X` |
| `library.utils` (3) | `diff_utils`: `first` from `library.utils.collection_utils`; `os_utils`: `FormerTuple` from `library.utils.collection_utils` |
| `classification.signals` (2) | `classification_hooks_import_notifications`: import `send_prepared_discordance_notifications` from its defining module |
| `analysis.models` (2) | `models_candidate_search` imports the defining modules |
| `classification.views.exports` (4) | `classification_export_formatter_mvl` / `_lab_compare` import from the module that defines what they need, rather than `classification_export_view` |
| `snpdb.models.models_user_settings` (2) | line 26 `from snpdb.models import AlleleOriginFilterDefault, UserAwards` → name the defining modules |

`sync` (8 modules) is the one aggregator cycle with real intent behind it:
`sync/models/models.py` inline-imports `run_sync`, and `sync_run.py` does
`from sync.alissa import *  # to get decorators to register`. Replace the import-for-side-effect
with explicit registration in `sync/apps.py::ready()`, and the whole component collapses.

Expected result: **top-level module cycles go from 6 to 0.**

## Stage 4 — The model-level knots

These are real design coupling, not import hygiene. Each is worth doing on its own merits, but
none is urgent and each is a bigger change than Stages 1–3.

**`analysis.models.nodes.sources._stats_cache` ↔ the source nodes (7 modules)** — the cheapest of
the four. `_stats_cache` inline-imports `CohortNode`, `PedigreeNode`, `SampleNode`, `TrioNode` at
line 120 to dispatch on node type. Invert it: let each node class hand `_stats_cache` what it
needs (or register itself), so the dependency runs one way only.

**`patients.models` ↔ `annotation` phenotype matching (4 modules)** — `patients.models` imports
`HasPhenotypeDescriptionMixin` from `annotation`, while `annotation.models.models_phenotype_match`
and `annotation.phenotype_matching` both import `patients.models.Patient`. The mixin is the odd
one out; moving `HasPhenotypeDescriptionMixin` into `patients` (it is a patients-side mixin that
happens to be defined in annotation) leaves a clean `annotation → patients` edge.

**`classification.models.classification` ↔ `discordance_models` ↔ `clinical_context_models` (3)** —
already handled by one inline import at `classification/models/classification.py:203`. Genuine
bidirectional domain coupling (a Classification knows its DiscordanceReports, a DiscordanceReport
knows its Classifications). Leave as is; the inline import is the correct tool here.

**`genes.hgvs` ↔ `genes.models` ↔ `annotation.models.models` ↔ `snpdb.models.models_variant` ↔
`snpdb.clingen_allele` (12 modules)** — the HGVS/Variant/Allele knot, the hardest in the codebase
and the one that generates the most inline imports (`genes.hgvs` is inline-imported 10 times).
Untangling it means separating HGVS *resolution* (needs transcripts, ClinGen, a genome build) from
the HGVS *value types*, which is a substantial piece of work. Recommendation: leave it, and freeze
it with a contract (Stage 5) so it does not grow.

## Stage 5 — Enforce, so this does not regress

`import-linter` (installed in Stage 0) reads a config file and fails a CI run when a contract is
violated. Add `.importlinter` at the repo root:

```ini
[importlinter]
root_packages =
    analysis
    annotation
    classification
    ...
exclude_type_checking_imports = True

[importlinter:contract:library-is-a-base-layer]
name = library must not import any app
type = forbidden
source_modules = library
forbidden_modules =
    analysis
    annotation
    classification
    eventlog
    genes
    snpdb
    ...
ignore_imports =
    library.preview_request -> variantgrid.perm_path

[importlinter:contract:no-new-module-cycles]
name = No import cycles within an app
type = independence
modules = ...
```

Notes on scoping the contracts:

- `variantgrid` is the project package, not an app. Every `urls.py` importing
  `variantgrid.perm_path` and every `tasks.py` importing `variantgrid.celery` is normal Django
  layout, and `variantgrid/views.py` and `variantgrid/deployment_validation/` legitimately sit
  above the apps. Exempt `variantgrid` from the app-level contracts entirely.
- Tests are leaves. Keep them out of scope so a test fixture importing another app never fails a
  contract.
- Start with contracts that pass *today* (after Stages 1–3) rather than aspirational ones, so the
  build stays green and the file is a ratchet.

Add it to the lint script alongside pylint and ruff:

```bash
lint-imports --config .importlinter
```

Optionally also flip `cyclic-import` on in `config/pylint3.rc` once Stage 3 lands, though
import-linter gives better signal for this specific job.

## Recommended order and rough sizing

| Stage | Effort | Payoff |
|---|---|---|
| 1. `admin_utils` + datatables → `library` | ~58 mechanical import updates | breaks 2 app cycles, removes snpdb from 17 apps' admin layer |
| 3. Aggregator cycles | ~25 import-line edits across 6 packages + `sync/apps.py` | top-level module cycles 6 → 0 |
| 5. `.importlinter` + CI | small | stops regression — do this as soon as 1 and 3 land |
| 2. `library` base layer | 3 module moves + ~35 import updates + a `vcf_utils` decision | `library` becomes a genuine leaf |
| 4a. `_stats_cache` inversion | contained to `analysis/models/nodes/sources/` | 7-module cycle gone |
| 4b. `HasPhenotypeDescriptionMixin` move | contained to `patients` + `annotation` | 4-module cycle gone |

Explicitly out of scope: the 23-app SCC (inherent to cross-app FKs), the HGVS/Variant/Allele
12-module knot, and the classification discordance triangle. Those are documented above so the
next person does not re-derive them.
