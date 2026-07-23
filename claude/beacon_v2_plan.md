# Beacon v2 Genomic Data-Sharing Endpoint — Implementation Plan

Implementation plan for GitHub issue **#1661** — expose VariantGrid variants via a
GA4GH **Beacon v2** (v2.2.0, July 2025) REST API. Split from #40. Sibling of the
MatchMaker Exchange plan (`claude/matchmaker_exchange_plan.md`, #1662).

## 1. Background

Beacon is a **variant-centric** discovery protocol: "does this database contain an
allele at *chrom:pos ref>alt* on *assembly*?" Responses are tiered by requester
authorisation:

- **boolean** — yes/no presence (anonymous, lowest disclosure)
- **count** — number of observations / datasets (anonymous or semi-trusted)
- **record** — full per-variant / per-dataset detail (authorised users)

Beacon v2 replaced the dead v1. It uses proper `/api/...` paths (so the old #40
"must live at `/`" blocker is gone — hosting is settled in §3.2: `/beacon/` on the
existing host, no subdomain), and defines a **framework** (info/config/map endpoints)
plus **model** endpoints (`g_variants`, and optionally `biosamples`, `individuals`,
`cohorts`, `datasets`).

- Spec / docs: https://docs.genomebeacons.org/
- OpenAPI: https://github.com/ga4gh-beacon/beacon-v2
- Reference impls: [`EGA-archive/beacon2-ri-api`](https://github.com/EGA-archive/beacon2-ri-api),
  [`ga4gh/ga4gh-starter-kit-beacon`](https://github.com/ga4gh/ga4gh-starter-kit-beacon)

Beacon is **bidirectional** here: we both **serve** a Beacon (inbound — §4–§8) and
**query** other Beacons (outbound client — §9). The outbound half lets a user viewing
a variant see which external Beacons report it — a live read-enrichment, not a push of
our data (only a public coordinate leaves). It mirrors the outbound half of the
MatchMaker Exchange app (`mme/client.py`), which is why §9 reuses that shape directly.

**Federation** (registering *our* public URL in the ELIXIR Beacon Network / EGA
registry so others discover us) is a separate, *later, optional* step — there is no
central server to join. Out of scope for the initial build. Note this is distinct from
outbound querying (§9): federation advertises us; the outbound client consumes others.

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

Follow the split the MME app established (commit `e6397dc3a`): the **example / off**
defaults live in `default_settings.py`; the **real, enabled** values live only in
`vgaws.py` (variantgrid.com). Beacon is **off everywhere by default** and turned on by a
single env file.

### 3.1 `default_settings.py` — example, off

Add near other `*_ENABLED` flags (e.g. `SEQAUTO_ENABLED` at `:332`). Placeholder
`organization`/`contact_url` here are overridden per-deployment (same as MME's blank
`MME_CONTACT` template):

```python
BEACON_ENABLED = False  # GA4GH Beacon v2 genomic data-sharing endpoint (#1661). On only in vgaws.py.
BEACON_CONFIG = {
    # Static metadata served by the framework endpoints (/info, /service-info, /map).
    # Placeholders — a deployment that enables Beacon overrides these in its env file.
    "beacon_id": "org.variantgrid.beacon",       # reverse-DNS id; override per deployment
    "name": "VariantGrid Beacon",
    "api_version": "v2.0.0",
    "environment": "prod",                        # prod | test | dev | staging
    "organization": {
        "id": "",
        "name": "",
        "welcome_url": "",
        "contact_url": "",                        # e.g. "mailto:..."
    },
    "default_granularity": "boolean",             # tier for anonymous requests
    "max_granularity": "record",                  # ceiling; per-request clamped by auth
}

# Small-count anonymity floor for the observations dataset (§5.2): exact counts below
# this are suppressed (the resultSet drops to boolean presence) to limit membership-
# inference risk. Behind a setting so the second security pass (open question #2) can
# tune it — best-practice Beacons avoid exposing exact small counts. The classification
# dataset (§5.5) is exempt (each record is already a deliberate public share).
BEACON_MIN_REPORTABLE_COUNT = 5
```

No explicit `BEACON_ENABLED = False` is needed in the Shariant (or any other) env files —
`False` is the default and they simply never turn it on. Only the enabled deployments
(§3.2 — `vgaws.py` and `vgtest2.py`) turn it on.

### 3.2 Enabled deployments — `vgaws.py` (prod) and `vgtest2.py` (test)

Two deployments turn Beacon on: **`vgaws.py`** (variantgrid.com, production) and
**`vgtest2.py`** (test.variantgrid.com — for manual testing). Everywhere else it stays
off by default. Each sets `BEACON_ENABLED = True` and overrides identity (mirrors how
`MME_CONTACT` gets real values only where MME is enabled):

```python
# vgaws.py — production
BEACON_ENABLED = True
BEACON_CONFIG = {
    **BEACON_CONFIG,                              # inherit defaults, override identity
    "beacon_id": "com.variantgrid.beacon",
    "environment": "prod",
    "organization": {
        "id": "variantgrid",
        "name": "VariantGrid",
        "welcome_url": "https://variantgrid.com/",
        "contact_url": "mailto:admin@variantgrid.com",
    },
}
```

```python
# vgtest2.py — test box (test.variantgrid.com), for manual testing
BEACON_ENABLED = True
BEACON_CONFIG = {
    **BEACON_CONFIG,
    "beacon_id": "com.variantgrid.test.beacon",
    "environment": "test",                        # spec signal: non-production instance
    "organization": {
        "id": "variantgrid",
        "name": "VariantGrid (test)",
        "welcome_url": "https://test.variantgrid.com/",
        "contact_url": "mailto:admin@variantgrid.com",
    },
}
```

`environment: "test"` is the spec's own flag that the instance is non-production, so
anything crawling or registering it knows not to treat its data as real. (Note
`vgtest2.py` sets `HEALTH_CHECK_ENABLED = False`, so the §8 Slack line will not post
there — inspect `BeaconInboundQuery` rows via admin when testing on the test box.)

**No subdomain.** Beacon v2's `/api/...` paths remove the v1 "must live at `/`" blocker
(the original reason #40 floated `beacon.shariant.org`), so Beacon is served at
`/beacon/` on the existing `variantgrid.com` host — no new vhost, TLS cert, or
`ALLOWED_HOSTS` entry. A dedicated `beacon.variantgrid.com` is purely cosmetic for
federation registration and can be added later without app changes. (Resolves open
question #3.)

## 4. App structure

Create a **new Django app `beacon`** (self-contained; mirrors how `seqauto` etc. are
gated). This keeps the spec-conformance surface out of `snpdb`.

**Inbound serving:**
- `beacon/apps.py`, `beacon/__init__.py`
- `beacon/urls.py` — framework + model endpoints (see §6)
- `beacon/views_rest.py` — DRF views (project uses plain `APIView` / generics, **no
  routers** — follow `snpdb/views/views_rest.py`)
- `beacon/serializers.py` — DRF serializers producing the Beacon response envelope
- `beacon/schema.py` — static framework metadata (entry types, filtering terms,
  service-info) built from `BEACON_CONFIG`
- `beacon/response.py` — helpers to build the standard `meta` + `responseSummary` +
  `response` envelope and clamp granularity by auth
- `beacon/signals/beacon_health_check.py` — nightly Slack health-check receiver (§8),
  imported from `beacon/apps.py` `ready()` (mirrors `sync/apps.py`)

**Outbound querying (§9):**
- `beacon/client.py` — `BeaconClient` (`requests` wrapper, mirrors `mme/client.py`)
- `beacon/tasks.py` — outbound fan-out task (network-bound, `web_workers`)
- `beacon/views.py` — `external_beacons` HTML-fragment view for the variant page (§9.4)

**Shared / cross-cutting:**
- `beacon/variant_mapping.py` — `Variant`/`Allele` ↔ Beacon `g_variant` translation,
  used by both the inbound view and the outbound client (§5.1, §9.2)
- `beacon/models.py` — `BeaconInboundQuery` audit row (§8, mirrors `MMEInboundQuery`);
  optional `BeaconQueryCache` for outbound results (§9.5). No `BeaconDataset` registry —
  rejected (§5.4).
- `beacon/tests/` — see §11

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

This translation lives in `beacon/variant_mapping.py` and is shared by the inbound
`g_variants` view and the outbound `BeaconClient` (§9.2).

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

**Small-count anonymity floor (observations only).** To limit membership-inference /
re-identification risk (the Shringarpure–Bustamante concern), exact counts below
`settings.BEACON_MIN_REPORTABLE_COUNT` are not reported for this dataset — the resultSet
drops to boolean presence. It sits behind a setting so the exact policy (including
whether to also gate *presence*, not just the count) can be tuned in the planned second
security pass (open question #2). The classification dataset (§5.5) is exempt: each
record is already a deliberate record-level public share.

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

### 5.4 Rejected: per-dataset `BeaconDataset` registry

A precomputed per-dataset registry (a `BeaconDataset` mapping a Beacon dataset id → a
permissioned VG `Cohort` + its own `VariantZygosityCountCollection`) was considered for
faster positive-path counts, but **we are not building it**. Its only real benefit is
optimisation, and it would have to be invalidated whenever *any* sample's permission
changed — a stale-cache liability for little gain, since the §5.2 path is dominated by
the cheap "not in the DB at all" fast-fail (global VZC gate) and only does the cohort
walk on the rare positive hit. Revisit only if real traffic shows the positive-path
count is a bottleneck. (Resolves open question #1.)

### 5.5 Two datasets served: observations + public classifications

Beacon serves **two datasets**, returned as independent `response.resultSets[]`:

1. `variantgrid_observations` — sample genotype presence/count (§5.2), gated by `Sample`
   Guardian perms (anonymous → `public` group).
2. `variantgrid_classifications` — variants that carry a **published `ShareLevel.PUBLIC`
   classification**, i.e. the *same* consented set MME submits/serves. This gives Beacon
   a curated, explicitly-shared dataset alongside the raw observation counts, and
   resolves the consent-surface question (open question #5) by modelling both rather than
   picking one.

**Reuse (from `mme/`):** the classification dataset is a thin Beacon wrapper over code
MME already has:
- eligible set: `mme_eligible_classifications`
  (`mme/serializers/patient_profile.py`) — the published + `ShareLevel.PUBLIC` filter.
- coordinate → allele/variant: `classification_genomic_feature`, plus the `Allele` build
  helpers (`.grch37()`/`.grch38()`/`variant_for_build()`) — classifications attach to
  `Allele` (build-independent via `ImportedAlleleInfo`), so a per-assembly Beacon query
  maps *coordinate → Allele → its public classifications*.
- `record`-tier attribution: `mme/contact.py` (`lab_mme_contact` /
  `mme_contact_for_classification`) — extract to a generic `lab_contact()` shared by both.

**Tiers for the classification dataset:**
- *boolean* — a public classification exists for the allele.
- *count* — number of public classification **records** for the allele. This is not the
  number of labs: a single lab may classify the same allele multiple times (separate
  samples/patients, or re-classification over time), each a distinct `Classification`
  record, so the count can exceed the distinct-lab count. No k-anonymity floor needed:
  each record is already a deliberate record-level public share (unlike genotype counts,
  §5.2 / open question #2).
- *record* — ACMG clinical significance + condition/phenotype + owning-lab contact.
  Phenotype (condition-under-curation + any linked HPO terms) reuses MME's ontology
  routing (`classification_ontology_slots`, `mme/serializers/patient_profile.py`), the
  same features/disorders mapping MME sends outward. This is the same data already
  flowing to ClinVar/MME under that share level, so exposing it is consistent with
  existing consent, not new disclosure.
  **Patient phenotype (extra gate):** when the classification's linked `Patient` is
  itself shared publicly (patient-level Guardian permission → `public` group, the same
  gate as everywhere else), also include that patient's HPO ontology phenotype terms in
  the record-tier phenotype. MME already sources HPO features from a linked patient, so
  this reuses the same path — just gated on the patient being public. With no public
  linked patient, phenotype stays limited to the classification's condition.

Advertise both datasets via the `/datasets` entry type (§6, phase 2); phase 1 can return
both resultSets from `g_variants` with static dataset metadata.

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
4. Classification dataset (§5.5): resolve coordinate → `Allele` → public classifications.
5. Clamp granularity to `min(requested, allowed_by_auth, BEACON_CONFIG.max)`.
6. Build the Beacon envelope (`beacon/response.py`); record the query for audit (§8).

## 7. Response envelope

Every response wraps a standard `meta` (beaconId, apiVersion, returnedGranularity,
receivedRequestSummary). Query responses add `responseSummary` (`exists`,
`numTotalResults`) and, above boolean, a `response.resultSets[]` with **one entry per
dataset** (§5.5: `variantgrid_observations` + `variantgrid_classifications`).
`beacon/response.py` centralises building this and the granularity clamp so no view
leaks more than the tier allows.

## 8. Inbound query audit & health check

Inbound Beacon traffic is expected to be **low volume**, so — like MME's
`MMEInboundQuery` — we store the **full `request_json`** for each query. That is cheap at
this volume and genuinely useful: you can inspect *what people are actually looking for*,
not just how many queries arrived. The nightly Slack health check aggregates these rows
over its window.

### 8.1 Audit row (`BeaconInboundQuery`) — mirrors `MMEInboundQuery`

One row per inbound `g_variants` query, storing the full request plus a few denormalised
columns so the health check can break down by type without digging into JSON:

```python
class BeaconInboundQuery(models.Model):
    """ Audit row for one inbound Beacon /g_variants query: the exact request we were
        sent (so we can see what people search for) plus a small response summary for
        health-check aggregation. Mirrors mme.models.MMEInboundQuery. """
    created = models.DateTimeField(default=timezone.now)
    request_json = models.JSONField()                   # exact Beacon query we received
    granularity = models.CharField(max_length=8)        # boolean | count | record (returned)
    authenticated = models.BooleanField(default=False)  # anon (public tier) vs authed requester
    observations_exists = models.BooleanField(default=False)
    observations_count = models.IntegerField(default=0)
    classifications_exists = models.BooleanField(default=False)
    classifications_count = models.IntegerField(default=0)
```

Written by the `g_variants` view. A matching `logging.info` line is still worth emitting
for live tailing, but the DB row is now the primary record.

### 8.2 Health-check receiver (Slack)

Connect a receiver to `health_check_overall_stats_signal` (framework in
`library/health_check.py`; pattern from `sync/signals/sync_health_check.py`), registered
via `beacon/apps.py` `ready()`. **Gated on `settings.BEACON_ENABLED`** so it only appears
in the Slack post where Beacon is actually on:

```python
@receiver(signal=health_check_overall_stats_signal)
def beacon_health_check(sender, health_request: HealthCheckRequest, **kwargs):
    if not settings.BEACON_ENABLED:
        return []
    qs = BeaconInboundQuery.objects.filter(created__gte=health_request.since,
                                           created__lt=health_request.now)
    total = qs.count()
    by_gran = dict(qs.values_list("granularity").annotate(n=Count("pk")))  # boolean/count/record split
    obs_hits = qs.filter(observations_exists=True).count()
    cls_hits = qs.filter(classifications_exists=True).count()
    return [HealthCheckRecentActivity(
        emoji="📡", name="Beacon queries", amount=total, sub_type="received",
        extra=(f"granularity {by_gran}; hits — observations {obs_hits}, "
               f"classifications {cls_hits}"),
        stand_alone=True,
    )]
```

`HealthCheckRecentActivity` is the right stat type — it is the "N things happened in the
window" shape (`is_recent_activity=True`), so the number of queries and the
type-breakdown ("results of what type went out over a period") land in the nightly post.
The `health_request` window supplies the period; no extra bookkeeping needed.

### 8.3 Optional live view (later)

The health-check framework notes it "might be extended to show the data live on a
webpage". Since `BeaconInboundQuery` is a normal model, a small admin/`DatatableConfig`
listing (counts over selectable windows) is a cheap later add — out of scope for MVP.

## 9. Outbound Beacon querying (client)

Let a user viewing a variant see which **external** Beacons report it. We send only a
public coordinate (`referenceName`/`start`/`referenceBases`/`alternateBases`/
`assemblyId`) and read back presence/count; nothing sensitive leaves. This is a live
read-enrichment, **not** a submission of our data — so unlike MME outbound there is no
draft/confirm workflow, no per-node submission rows, and no curator notification.

### 9.1 Reuse from `mme/` (outbound half)

The transport mirrors `mme/client.py` directly:

| MME piece | Beacon outbound analog |
|---|---|
| `MMEClient` (`requests` wrapper: per-node `base_url`/`api_version`/headers, `raise_for_status`) | `BeaconClient` — same shape, issues `POST /g_variants` (or `GET`) |
| `settings.MME_NODES` (per-node dict, token via `get_secret(mandatory=False)`) | `settings.BEACON_QUERY_NODES` — same structure; **token optional** (boolean/count tier is anonymous at most Beacons; only some `record` tiers need a registered token) |
| `@app.task(queue='web_workers')` network task | Same queue for the fan-out task (network-bound) |
| Admin-notify-on-failure (`AdminNotificationBuilder`) | Same, on per-node error |

**Deliberately not reused:** `MMESubmission` draft/confirm/submit, persisted
`MMEMatchResult` candidate rows, curator `Message` notification. Those exist because
MME *pushes our consented data*; a Beacon query pushes only a coordinate.

### 9.2 Shared coordinate mapping (symmetry with inbound)

Just as MME shares `mme/serializers/patient_profile.py` between its inbound matcher and
outbound client, the `Variant`/`Allele` ↔ Beacon `g_variant` translation lives in one
**`beacon/variant_mapping.py`** (§5.1), used by **both**:
- the inbound `g_variants` view (Variant → response), and
- the outbound `BeaconClient` (Variant → query params; response → parsed presence/count).

One coordinate-mapping module, both directions.

### 9.3 Settings

Same split as §3 — example/off in `default_settings.py`, real values (and enable) in
`vgaws.py`:

```python
# default_settings.py — example, off. Outbound: query external Beacons from the variant
# page. Off by default; vgaws.py turns it on, and can turn it back off if remote Beacons
# are slow/unreliable (§9.5).
BEACON_OUTBOUND_ENABLED = False
BEACON_QUERY_TIMEOUT = 5          # per-node seconds; keep small — remote Beacons vary wildly
BEACON_QUERY_CACHE_DAYS = 7       # cache each (variant, node) result; live-refresh on expiry
BEACON_QUERY_NODES = {}           # populated in vgaws.py
```

```python
# vgaws.py — real values, enabled. base_url + api_version are plain public config; token
# (if a node requires one for count/record) is the only secret, via get_secret(mandatory=False).
BEACON_OUTBOUND_ENABLED = True
BEACON_QUERY_NODES = {
    # Configure whichever beacons/aggregators we confirm are LIVE and reachable — the
    # current genomebeacons.org guide documents no single official central query endpoint,
    # so there is no canonical aggregator URL to hard-code. Prefer an aggregator (one
    # gateway that fans out to member beacons) as a single entry when a reachable one is
    # confirmed; otherwise list specific beacons directly. Open (boolean/count) tiers need
    # no token; protected tiers use OIDC. Verify reachability before enabling (the older
    # ELIXIR Beacon Network host has been unreliable).
    # "some_beacon": {"base_url": "https://<confirmed-live-host>", "api_version": "v2.0.0",
    #                 "token": get_secret("BEACON.some_beacon_token", mandatory=False)},
}
```

Because the ecosystem is federated with **no documented central query endpoint**, the
node list is deployment-curated: we point at whatever aggregator/beacons we have verified
are live. An aggregator, when available, collapses to a single entry that fans out; absent
one, we list specific beacons. Confirm the exact API base + version from each target's own
`/service-info` before wiring it. (Registering *our* beacon so others' queries reach us is
the separate manual step in §10, phase 3.)

For manual testing, `vgtest2.py` (test.variantgrid.com) can also set
`BEACON_OUTBOUND_ENABLED = True` with a `BEACON_QUERY_NODES` pointing at any
confirmed-live beacon — including our own prod beacon — to exercise the variant-page
section end to end.

### 9.4 Variant-page section (async `.load()`)

Follow the existing lazy-load pattern in
`variantopedia/templates/variantopedia/variant_details.html` — e.g. line 399
`$("#variant_sample_information").load('{% url ... %}')`, and the `gene_coverage`
`{% settings_value %}` + `{% if %}` gate (`:246`):

1. Add `beacon/views.py:external_beacons(request, variant_id, genome_build_id)` — a
   normal HTML-fragment view (not the DRF `views_rest.py` inbound endpoint). It reads
   the cache (§9.5), and for uncached/expired nodes runs the fan-out, then renders a
   small table (node → exists / count / link-out).
2. Add a URL `path('external_beacons/<int:variant_id>/<int:genome_build_id>', ...)`.
3. In `variant_details.html`, gate a `#external_beacons` div behind the flag and
   `.load()` it **after** page render so remote latency never blocks the page:
   ```django
   {% settings_value 'BEACON_OUTBOUND_ENABLED' as beacon_outbound_enabled %}
   {% if beacon_outbound_enabled %}
     <div id="external_beacons"></div>
     <script>$(function(){ $("#external_beacons").load('{% url 'external_beacons' variant.id genome_build.pk %}'); });</script>
   {% endif %}
   ```
   Flipping `BEACON_OUTBOUND_ENABLED` off removes the section entirely — the "turn it
   off if it's too slow" lever.

### 9.5 Caching & timing

Remote Beacon latency is the whole reason this is async + cached:
- **One direct node**, exact lookup: ~sub-second–2s (indexed presence lookup; network
  round-trip dominates).
- **Fan-out across a few configured nodes in parallel:** bounded by slowest node +
  `BEACON_QUERY_TIMEOUT` → a few seconds.
- **Beacon Network aggregator:** commonly 5–20s, tail cases 30s+ (it waits on many
  heterogeneous member Beacons, some slow/down).

Therefore:
- Query nodes **concurrently** with a hard per-node `BEACON_QUERY_TIMEOUT`; a slow/dead
  node yields "unavailable" for that row, never blocks the others.
- **Cache** each `(variant, node)` result for `BEACON_QUERY_CACHE_DAYS` (thin
  `BeaconQueryCache` model, or Redis via the existing cache) so repeat views are
  instant; refresh live on expiry. This is the "cached with live refresh" default of
  open question #6.
- Because the section is `.load()`ed after render, even a 20s network query only delays
  that one panel, not the variant page.

### 9.6 Testing

- Extend `django.test.TestCase`; **mock `requests`** (mirror `mme/tests/fakes.py` /
  `test_client.py`) — never hit live Beacons in tests.
- `BeaconClient` builds correct `g_variants` params from a `Variant` and parses
  boolean/count responses (shared `variant_mapping.py`, exercised both directions).
- Fan-out: one node timing out / erroring still returns the others; results cached and
  the second call is served from cache (no second `requests` call).
- `external_beacons` view renders the fragment and is gated by `BEACON_OUTBOUND_ENABLED`.

## 10. Phasing

**Phase 1 — MVP (this issue):**
- New `beacon` app, flag, `PUBLIC_PATHS`, off by default; enabled on variantgrid.com
  (prod) and test.variantgrid.com (manual testing) per §3.2. Shariant never enables it.
- Framework endpoints (static from `BEACON_CONFIG`).
- `g_variants` query supporting **boolean + count + record** granularity for the
  `genomicVariant` entry type, via the §5.2 two-stage lookup (global
  `VariantZygosityCount` gate → permission-scoped `CohortGenotype` count) using the
  core extracted from `VariantSampleInformation` (§5.3), tiered via `filter_for_user`.
- Second dataset `variantgrid_classifications` (§5.5): boolean/count/record over the
  published `ShareLevel.PUBLIC` set, reusing `mme_eligible_classifications` +
  `classification_genomic_feature` + the shared `lab_contact()`. Low-cost given the MME
  reuse, and it delivers the richer (significance-bearing) record tier from day one.
- Inbound audit + health check (§8).
- drf-spectacular docs; tests (§11), including a per-dataset resultSet assertion and the
  classification `record`-tier significance/contact payload.

**Phase 2 — richer model (later):**
- `datasets` / `cohorts` entry-type endpoints with metadata.
- `filtering_terms` backed by real ontology (`ontology` app).
- gnomAD/population frequency in record responses (from `annotation`
  `VariantAnnotation`).

**Phase 3 — individuals/phenotypes + federation (later):**
- `individuals` / `biosamples` entry types (needs `patients` + HPO mapping).
- **Register our beacon so others discover it** — the manual federation step (an operator
  process, not code), and lightweight in current practice. Per the
  [genomebeacons.org Deployment Guide](https://genomebeacons.org/Deployment-Guide/),
  discovery is a **self-report**: submit our public `/beacon/` endpoint via the registration
  form linked there. To be usable once listed: keep our framework endpoints
  (`/info`, `/service-info`, `/map`, `/filtering_terms`) reachable and spec-conformant, and
  wire OIDC only if we ever expose protected (non-public) tiers — the public boolean/count
  tier needs none. There is no single global server to "join"; if a live network aggregator
  (e.g. an ELIXIR/EGA Beacon Network, when reachable) is confirmed, register with it too so
  federated queries reach us. (EGA's own beacons use the B2RI reference implementation; ours
  only needs to be conformant and reachable.)

**Outbound querying (§9)** is an independent track that can land alongside phase 1 — it
only needs the shared `beacon/variant_mapping.py` coordinate translation, not the inbound
serving stack. Ship it whenever the client + variant-page section are ready.

## 11. Testing

- Extend `django.test.TestCase`; URL smoke tests via `URLTestCase`
  (`library/django_utils/unittest_utils.py`).
- Framework endpoints return spec-shaped JSON.
- `g_variants`: known variant present/absent; boolean vs count vs record output.
- **Permission tiering** (most important): anonymous sees only public-group datasets;
  a variant that exists *only* in a private dataset returns `exists:false` /
  count 0 to anonymous but is visible to an authorised user. Guards the leak in §5.2.
- Classification dataset: per-dataset resultSet present; `record`-tier significance +
  phenotype + lab contact payload (§5.5).
- Outbound client + variant-page section: see §9.6.
- Fake data via `snpdb/tests/utils/`, `annotation/tests/test_data_fake_genes.py`.
- Validate a sample response against the Beacon v2 OpenAPI schema.
- **Post-deploy conformance check:** once the beacon is deployed at a reachable URL, run
  the EGA **beacon-verifier** against it (demo: https://beacon-verifier-demo.ega-archive.org/)
  — it checks our framework + `g_variants` endpoints against the Beacon v2 spec. This is a
  deployment/CD step, not a unit test (it needs a live URL), so it lands after the
  variantgrid.com rollout and is worth re-running whenever the response shape changes or
  before registering with any network (§10, phase 3).

## 12. Open questions

1. **Dataset model** — RESOLVED (§5.4): use the §5.2 global-VZC-gate + per-sample cohort
   walk; do **not** build the `BeaconDataset` per-dataset registry. It is only an
   optimisation and would go stale on any sample permission change, for little gain given
   the fast-fail-dominated workload. Revisit only if the positive-path count proves a
   bottleneck.
2. **Public count floor** — RESOLVED for now (§5.2, §3.1): suppress exact small counts on
   the observations dataset via `BEACON_MIN_REPORTABLE_COUNT` (best practice), behind a
   setting; the classification dataset is exempt. The exact policy (and whether to gate
   presence too) is revisited in the planned second security pass.
3. **Hostname** — RESOLVED (§3.2): serve at `/beacon/` on the existing host, **no
   subdomain**. Beacon v2's `/api/...` paths make a dedicated vhost unnecessary; a
   cosmetic `beacon.variantgrid.com` for federation can come later without app changes.
4. **Auth for the `record` tier** — RESOLVED: no OIDC scoping needed yet. The existing
   DRF token/session auth plus Guardian permission tiering (anonymous → public group) is
   sufficient; revisit if a future consumer needs scoped access.
5. **Consent surface (inbound)** — RESOLVED (§5.5): expose *both* the sample-permission
   surface (§5.2) and the `ShareLevel.PUBLIC` classifications (the MME set) as **two
   separate Beacon datasets** rather than choosing one. Beacon v2's per-dataset
   `resultSets[]` model carries them independently. The classification dataset's `record`
   tier exposes clinical significance **and** phenotype (§5.5) — publicly sharing a
   classification consents to sharing both; it is the same set already flowing to
   ClinVar/MME. Additionally, when the linked `Patient` is itself shared publicly, its
   HPO phenotype terms are included too (§5.5, gated on patient-level public permission).
6. **Outbound: live vs cached** (§9) — plan defaults to **cached with live refresh**.
   Confirm cache TTL and whether a manual "refresh" control is wanted.
