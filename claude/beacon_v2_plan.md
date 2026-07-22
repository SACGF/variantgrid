# Beacon v2 Genomic Data-Sharing Endpoint — Implementation Plan

Implementation plan for GitHub issue **#1661** — expose VariantGrid variants via a
GA4GH **Beacon v2** (v2.2.0, July 2025) REST API. Split from #40. Sibling of the
MatchMaker Exchange plan (`claude/matchermaker_exchange_plan.md`, #1662).

## 1. Background

Beacon is a **variant-centric** discovery protocol: "does this database contain an
allele at *chrom:pos ref>alt* on *assembly*?" Responses are tiered by requester
authorisation:

- **boolean** — yes/no presence (anonymous, lowest disclosure)
- **count** — number of observations / datasets (anonymous or semi-trusted)
- **record** — full per-variant / per-dataset detail (authorised users)

Beacon v2 replaced the dead v1. It uses proper `/api/...` paths (so the old #40
"must live at `/`" blocker is gone — we can host at e.g. `beacon.shariant.org`),
and defines a **framework** (info/config/map endpoints) plus **model** endpoints
(`g_variants`, and optionally `biosamples`, `individuals`, `cohorts`, `datasets`).

- Spec / docs: https://docs.genomebeacons.org/
- OpenAPI: https://github.com/ga4gh-beacon/beacon-v2
- Reference impls: [`EGA-archive/beacon2-ri-api`](https://github.com/EGA-archive/beacon2-ri-api),
  [`ga4gh/ga4gh-starter-kit-beacon`](https://github.com/ga4gh/ga4gh-starter-kit-beacon)

**Federation** (ELIXIR Beacon Network / EGA registry) is a *later, optional* step —
it just registers our public URL; there is no central server to join. Out of scope
for the initial build.

## 2. Design decisions (from planning)

- **All three granularity tiers** (boolean → count → record) are in scope.
- **Permission-driven, anonymous allowed.** Anonymous external requests are served
  the **public** records only. This maps *directly* onto the existing
  `GuardianPermissionsMixin.filter_for_user()` behaviour: for `AnonymousUser` it
  falls back to the `PUBLIC_GROUP_NAME` ("public") group's permissions
  (`library/django_utils/guardian_permissions_mixin.py:93-94`,
  `PUBLIC_GROUP_NAME` at `default_settings.py:318`). No new permission machinery.
- **Counts come from `VariantZygosityCount` collections** for fast lookups
  (`snpdb/models/models_zygosity_counts.py`), *scoped by permission* — see §5.
- **Feature-flagged, off by default, off in Shariant** (issue comment requirement).

## 3. Feature flag & configuration

Add to `variantgrid/settings/components/default_settings.py` (near other `*_ENABLED`
flags, e.g. `SEQAUTO_ENABLED` at `:332`):

```python
BEACON_ENABLED = False  # GA4GH Beacon v2 genomic data-sharing endpoint (#1661)
BEACON_CONFIG = {
    # Static metadata served by the framework endpoints (/info, /service-info, /map)
    "beacon_id": "org.variantgrid.beacon",       # reverse-DNS id; override per deployment
    "name": "VariantGrid Beacon",
    "api_version": "v2.0.0",
    "environment": "prod",                        # prod | test | dev | staging
    "organization": {
        "id": "variantgrid",
        "name": "VariantGrid",
        "welcome_url": "https://variantgrid.com/",
        "contact_url": "mailto:...",
    },
    "default_granularity": "boolean",             # tier for anonymous requests
    "max_granularity": "record",                  # ceiling; per-request clamped by auth
}
```

Turn it **off** in the Shariant env files (`variantgrid/settings/env/shariant.py`,
`shariantcommon.py`, `sharianttest.py`, `shariantdemo.py`) — explicit
`BEACON_ENABLED = False`, even though it is the default, so intent is on the record.
`variantgrid.com` env can set `BEACON_ENABLED = True`.

## 4. App structure

Create a **new Django app `beacon`** (self-contained; mirrors how `seqauto` etc. are
gated). This keeps the spec-conformance surface out of `snpdb`.

- `beacon/apps.py`, `beacon/__init__.py`
- `beacon/urls.py` — framework + model endpoints (see §6)
- `beacon/views_rest.py` — DRF views (project uses plain `APIView` / generics, **no
  routers** — follow `snpdb/views/views_rest.py`)
- `beacon/serializers.py` — DRF serializers producing the Beacon response envelope
- `beacon/schema.py` — static framework metadata (entry types, filtering terms,
  service-info) built from `BEACON_CONFIG`
- `beacon/models.py` — see §5 (dataset registry; may be thin / config-only initially)
- `beacon/response.py` — helpers to build the standard `meta` + `responseSummary` +
  `response` envelope and clamp granularity by auth
- `beacon/tests/` — see §9

Wire-up:
- Add `"beacon"` to `INSTALLED_APPS` (`default_settings.py`) and to `APPS_WITH_URLS`
  (`variantgrid/urls.py:15`).
- Gate the URL include behind the flag. Two options; **prefer the explicit `if`**
  pattern already used for `INBOX_ENABLED`/`CONTACT_US_ENABLED`
  (`variantgrid/urls.py:53,63`):

  ```python
  if settings.BEACON_ENABLED:
      urlpatterns += [path('beacon/', include('beacon.urls'))]
  ```
  (Alternatively toggle via `URLS_APP_REGISTER["beacon"]`, but the `if` is clearer
  for a whole-app on/off.)

- **`PUBLIC_PATHS`**: Beacon must answer anonymous requests, but
  `GlobalLoginRequiredMiddleware` rejects them *before DRF runs*
  (`default_settings.py:579`, note at `:619-620`). Add the Beacon prefix to
  `PUBLIC_PATHS` (`default_settings.py:832`):
  ```python
  r'^/beacon/.*',
  ```
  The view still enforces tiering via `filter_for_user` (anonymous → public group),
  so "public path" ≠ "public data".

## 5. Data model mapping

### 5.1 Variant identity → Beacon `g_variant`

VariantGrid stores coordinates normalised/interned (no chrom/ref/alt strings on
`Variant`):

| Beacon field | VariantGrid source |
|---|---|
| `assemblyId` / `referenceName` | `GenomeBuild` (PK is name, e.g. "GRCh37"); `Contig.name` / `ucsc_name` (`models_genome.py`) |
| `start` (0-based) | `Variant.locus.position` (VG is 1-based → convert) |
| `referenceBases` | `Variant.locus.ref.seq` |
| `alternateBases` | `Variant.alt.seq` |
| `variantType` / `svlen` | `Variant.svlen` (symbolic/SV variants) |
| `variantInternalId` | `Variant.pk` (or `Allele.pk`) |
| stable/canonical id (`identifiers`) | `Allele.clingen_allele` (ClinGen Allele Registry ID) — natural cross-build id (`models_variant.py:53`) |

Query helpers already exist: `VariantCoordinate` / `Variant.get_from_string`
(`models_variant.py:224+`), `Variant.get_chrom_q()` (`:562`),
`Variant.get_contigs_q(genome_build)` (`:566`). Build-independence via `Allele`
(`.grch37()`, `.grch38()`, `variant_for_build()`).

### 5.2 Two-stage count: fast VZC gate → permission-scoped cohort walk

This is the crux of correct + fast counting. Two data sources, used in sequence:

**Stage 1 — quick lookup: `VariantZygosityCount`** (`snpdb/models/models_zygosity_counts.py`).
The **global** collection (`GLOBAL_ALIAS`, `get_global_germline_counts()` `:63`,
`annotate_global_germline_counts()` `:72`) stores one precomputed indexed
`ref/het/hom` row per variant. Beacon queries are dominated by "not in the DB at
all", so this gives a cheap early exit:

- If the variant has **no** VZC row / global count `== 0` → `exists: false`, return
  immediately, no `CohortGenotype` touched.

The global count includes **private** samples, so it is used only as a **negative
gate** — it can prove *absence* cheaply, but a positive global count must **not** be
reported to a requester as presence, or it leaks that a variant exists in private-only
data. Positive answers always come from stage 2.

**Stage 2 — permission-scoped exact count: `CohortGenotype`.** Only runs when stage 1
says "present". This is the existing logic in
`snpdb/variant_sample_information.py` (`VariantSampleInformation`), which already:

- gates by permission: `Sample.filter_for_user(user)` → `user_sample_ids`
  (`:28-30`), anonymous → `public` group;
- unpacks the packed per-VCF `CohortGenotype` row into per-sample zygosities
  (`_cohort_genotype_to_sample_genotypes`, `:101`);
- separates visible from hidden: `num_observations` vs `num_visible_observations` vs
  `num_invisible_observations` (`:60-67`) and per-zygosity `visible_zygosity_counts`.

`exists` = visible count `> 0`; the Beacon `count` granularity reports the visible
het+hom count. This is safe for every tier: presence/counts are only ever asserted
from samples the requester may read.

### 5.3 Reuse: extract the permission-aware core from `VariantSampleInformation`

`VariantSampleInformation` mixes a **reusable core** (the permission-gated per-variant
sample/zygosity lookup, §5.2 stage 2) with **presentation** (pandas `locus_counts_df`
`:186`, checkbox-formatted `visible_zygosity_counts` `:55-58`, phenotype-match graphs
`:73-76`). Refactor to pull the core out so Beacon and the existing Variantopedia view
(`variantopedia/views.py:743`) share it:

- Extract a function/class returning, for `(user, variant, genome_build)`: visible
  per-sample genotype rows + `{num_observations, num_visible, num_invisible}` +
  per-zygosity visible counts — **without** pandas/template concerns.
- `VariantSampleInformation` keeps only the presentation layer, calling the core.
- Beacon's `g_variants` view calls the same core for stage 2.

This is a small, low-risk refactor (single caller today) and avoids a second copy of
the fiddly `CohortGenotype`-array unpacking + permission logic.

### 5.4 Beacon datasets (optional finer model / scaling)

The §5.2 model gives **per-sample** permission granularity from one global VZC + a
cohort walk. If per-dataset advertising or faster positive-path counts are needed
later, add a thin **`BeaconDataset`** registry in `beacon/models.py` mapping a Beacon
dataset id → a permissioned VG `Cohort` + its own `VariantZygosityCountCollection`
(collections are named & individually partitioned, `:16`). Then per-dataset counts are
fully precomputed (no cohort walk) and permission is per-dataset. It also makes Beacon
exposure an **explicit opt-in** separate from general Guardian read perms ("readable in
the app" ≠ "advertised on Beacon"). Defer unless needed — start with §5.2.

## 6. Endpoint map (`beacon/urls.py`)

**Framework (mostly static, from `BEACON_CONFIG`):**

| Path | Purpose |
|---|---|
| `GET /` and `/info` | Beacon identity/metadata |
| `GET /service-info` | GA4GH service-info profile |
| `GET /configuration` | Beacon configuration object |
| `GET /entry_types` | Supported entry types (phase 1: `genomicVariant`) |
| `GET /filtering_terms` | Supported ontology filters (phase 1: minimal) |
| `GET /map` | Endpoint map |

**Model — genomic variants:**

| Path | Purpose |
|---|---|
| `GET/POST /g_variants` | Query by `referenceName`, `start`, `referenceBases`, `alternateBases`, `assemblyId`, range/bracket queries; returns tiered `responseSummary` + `response.resultSets` |
| `GET /g_variants/{id}` | Single variant detail |
| `GET /g_variants/{id}/biosamples` etc. | Deferred (later phase) |

Follow existing REST wiring: per-app `path('...')` in `beacon/urls.py`, `APIView`/
generics in `beacon/views_rest.py`, `@extend_schema` for drf-spectacular docs (already
configured; docs at `/api/schema`, `/api/docs`).

### Request → query flow (`g_variants`)

1. Parse Beacon params → `VariantCoordinate` / build a `Variant` queryset filtered by
   `get_contigs_q(genome_build)` + position/ref/alt (support exact and range queries).
2. **Stage 1 (fast gate)** — read global `VariantZygosityCount`. If absent/zero →
   `exists: false`, done (no cohort walk). See §5.2.
3. **Stage 2 (permission count)** — for present variants, call the extracted
   permission-aware core (§5.3): `Sample.filter_for_user(request.user)` (anonymous →
   public group) over `CohortGenotype` → visible het/hom count. `exists` = visible
   count > 0.
4. Clamp granularity to `min(requested, allowed_by_auth, BEACON_CONFIG.max)`.
5. Build the Beacon envelope (`beacon/response.py`).

## 7. Response envelope

Every response wraps a standard `meta` (beaconId, apiVersion, returnedGranularity,
receivedRequestSummary). Query responses add `responseSummary` (`exists`,
`numTotalResults`) and, above boolean, a `response.resultSets[]` with counts (a single
result set in the §5.2 model; per-dataset once §5.4 `BeaconDataset`s exist).
`beacon/response.py` centralises building this and the granularity clamp so no view
leaks more than the tier allows.

## 8. Phasing

**Phase 1 — MVP (this issue):**
- New `beacon` app, flag, `PUBLIC_PATHS`, Shariant off.
- Framework endpoints (static from `BEACON_CONFIG`).
- `g_variants` query supporting **boolean + count + record** granularity for the
  `genomicVariant` entry type, via the §5.2 two-stage lookup (global
  `VariantZygosityCount` gate → permission-scoped `CohortGenotype` count) using the
  core extracted from `VariantSampleInformation` (§5.3), tiered via `filter_for_user`.
- drf-spectacular docs; tests (§9).

**Phase 2 — richer model (later):**
- `datasets` / `cohorts` entry-type endpoints with metadata.
- `filtering_terms` backed by real ontology (`ontology` app).
- gnomAD/population frequency in record responses (from `annotation`
  `VariantAnnotation`).

**Phase 3 — individuals/phenotypes + federation (later):**
- `individuals` / `biosamples` entry types (needs `patients` + HPO mapping).
- Register with ELIXIR Beacon Network / EGA registry.

## 9. Testing

- Extend `django.test.TestCase`; URL smoke tests via `URLTestCase`
  (`library/django_utils/unittest_utils.py`).
- Framework endpoints return spec-shaped JSON.
- `g_variants`: known variant present/absent; boolean vs count vs record output.
- **Permission tiering** (most important): anonymous sees only public-group datasets;
  a variant that exists *only* in a private dataset returns `exists:false` /
  count 0 to anonymous but is visible to an authorised user. Guards the leak in §5.2.
- Fake data via `snpdb/tests/utils/`, `annotation/tests/test_data_fake_genes.py`.
- Validate a sample response against the Beacon v2 OpenAPI schema.

## 10. Open questions

1. **Dataset model** — start with the §5.2 global-VZC-gate + per-sample cohort walk
   (per-sample permission, no dataset registry), or go straight to the §5.4
   `BeaconDataset` per-dataset model with explicit opt-in? (Plan starts with §5.2.)
2. **Public count floor** — do we suppress/round small counts (k-anonymity) for the
   anonymous tier, as some Beacons do, or report exact counts?
3. **Hostname** — dedicated `beacon.<deployment>` vhost vs `/beacon/` on the main app.
   Plan assumes `/beacon/` path (spec-compatible either way).
4. **Auth for the `record` tier** — token vs session (both configured in
   `REST_FRAMEWORK`). Any OIDC-scoped access needed, or is Guardian-per-dataset
   sufficient?
