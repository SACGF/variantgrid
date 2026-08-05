# Plan: Variant details — samples grid expansion (#98)

**Issue**: [#98](https://github.com/SACGF/variantgrid/issues/98) — *Variant details - samples grid - classifications, VCF filter, and other builds / inheritance*

Requested features:

1. **Other genome builds** — the variant details page resolves an `Allele`, which links variants in other builds. Load samples from those builds into the same grid, with a genome build column.
2. **REST API + progressive loading** — replace the single server-rendered blob with per-build API calls so the current build renders immediately and other builds fill in as they arrive. Show a per-build loader/status with counts. Cacheable.
3. **Classifications** — show classification clinical significance pills (B/LB/VUS/LP/P) per sample row.
4. **VCF FILTER** — show FILTER, which is already on the `CohortGenotype` records being retrieved.

Inheritance (raised in a follow-up comment on #98) is tracked separately in #1693 — *sample membership (cohort/trio/pedigree) and trio inheritance*. This plan covers 1–4.

---

## 1. Current implementation

### Server

`snpdb/variant_sample_information.py`

- `VariantZygosityCounts` (line 18) — the permission-scoped core. Given `(user, variant, genome_build)` it:
  - Queries every `CohortGenotype` at the variant's **locus** (`_get_sample_values_for_variant_via_cohort_genotype`, line 87), restricted to VCF-backed cohorts.
  - Unpacks each packed `CohortGenotype` row into per-sample dicts (`_cohort_genotype_to_sample_genotypes`, line 106), annotating patient HPO/OMIM/MONDO via `StringAgg`.
  - Builds `locus_counter` (per-variant zygosity counts across the locus), `visible_rows` (rows for *this* variant, for samples in `Sample.filter_for_user`), and observation counts.
  - Also consumed by the Beacon `g_variants` endpoint — `beacon/datasets.py:82`. **This class's behaviour must not regress.**
- `VariantSampleInformation` (line 191) — presentation subclass: builds a pandas `locus_counts_df`, the zygosity checkbox counts, `hidden_samples_details`, and `patient_ids`.

`variantopedia/views.py:764` `variant_sample_information(request, variant_id, genome_build_name)` — renders the whole fragment, including `vsi.visible_rows` inlined as JSON.

### Client

`variantopedia/templates/variantopedia/variant_sample_information.html`

- Locus counts table (server rendered from the pandas DataFrame).
- Multi-allelic "other loci variants" list.
- A **free-jqGrid** (`#genotype-grid`) with `datatype: "local"` and `data: visible_rows` inlined — sample name, zygosity, enrichment kit, patient HPO/OMIM/MONDO, VCF, AF, AD, DP, PL.
- Zygosity checkbox filter (`filterGenotypes`), custom CSV export (`exportGenotypeGridCSV`), Plotly allele-frequency scatter + histogram, and server-rendered patient phenotype/OMIM graphs.

`variantopedia/templates/variantopedia/variant_details.html:399` loads the fragment:

```js
$("#variant_sample_information").load('{% url 'variant_sample_information' variant.id genome_build.pk %}');
```

The allele (and hence the other builds) is fetched separately and asynchronously in the same `$(document).ready` block via `variant_allele_for_variant` → `populateVariantAllele(data)`.

### Existing pieces we can reuse

| Need | Existing code |
|---|---|
| REST view patterns, per-user caching | `snpdb/views/views_rest.py:19` / `:62`; `@method_decorator([cache_page(MINUTE_SECS), vary_on_cookie], name='dispatch')` |
| Allele → variants in other builds | `Allele.variant_alleles()` (`snpdb/models/models_variant.py:97`); `AlleleSerializer` (`snpdb/serializers.py:69`) needs a `variants_by_build` field adding (§3.1) |
| Builds a variant belongs to | `Variant.genome_builds` (`snpdb/models/models_variant.py:675`), `Contig.get_genome_builds()` (`snpdb/models/models_genome.py:398`) — handles contigs shared between builds |
| Decode FILTER codes → filter ids | `VCFFilter.get_formatter(vcf)` (`snpdb/models/models_vcf.py:294`); per-sample FT via `SampleGenotype.get_vcf_filters` (`snpdb/models/models_cohort.py:703`) |
| Classifications for an allele, permission scoped | `ClassificationModification.latest_for_user(user, allele=..., published=True, **kwargs)` (`classification/models/classification.py:2370`) |
| Clin sig pill label/CSS | `clinical_significance_values` tag (`classification/templatetags/classification_tags.py:210`), CSS `.c-pill.cs-*` / `.c-pill.scs-*` in `global.scss:3066+`. A JS version already exists at `analysis/templates/.../view_candidate_search_run.html:133` |
| Classification ↔ sample link | `Classification.sample` FK (`classification/models/classification.py:469`) |
| DataTables + Buttons + HTML5 CSV export | Bundled globally — `uicore/templates/uicore/page/base.html:111` loads DataTables 1.13.1 with Buttons 2.3.2 and HTML5 export |

### Problems found in the current code

- **`vsi.has_hidden_samples` does not exist.** The template branches on it three times (lines 86, 383, and the `{% if %}` at 100 uses `hidden_samples_details`), but `VariantSampleInformation` never defines it — only the local variable inside `__init__`. So the "you can see N" message never renders. Fix as part of this work.
- The help panel at the top of the fragment says *"locus counts and sample genotypes are from VCFs in the current genome build only"* — must be removed once multi-build lands.
- `visible_rows` is inlined into the HTML with no cap. A variant present in thousands of samples produces a very large page. Worth a bounded response (see §6).
- The locus-counts pandas DataFrame is used only by this template; it can be dropped in favour of plain dicts once the table is rendered client-side.
- ~~**Shared contigs are mis-reported.**~~ **FIXED (done ahead of this plan).** `VariantZygosityCounts.__init__` scoped visible samples to the passed build but counted `num_observations` across *every* sample at the variant. For a variant on a contig shared by two builds (MT — 65 shared contigs in this deployment: MT plus unplaced scaffolds), samples from the other build's VCFs inflated `num_observations` while never entering `visible_rows`, so the template rendered **"You cannot see N samples"** — a permissions claim about samples the user could read.

  Fixed by dropping the `genome_build` constructor parameter and scoping to `variant.genome_builds`. `has_hidden_samples` (which the template branched on but never existed) added as a property. Call sites updated in `beacon/datasets.py` and `variantopedia/views.py`; template now names the builds actually searched. Regression tests in `snpdb/tests/test_variant_sample_information.py`. **Do not re-introduce the build parameter when implementing the rest of this plan.**

---

## 2. Proposed architecture

Split the fragment into **a shell** (containers + JS, rendered once) and **one REST call per distinct variant** on the allele, streaming rows into a single client-side grid.

```
variant_details.html
  └── $("#variant_sample_information").load(variant_sample_information)     <- shell, no data
         │
         ├── fetch api/variant_sample_genotypes/<this_variant_id>              (fires immediately)
         │
         └── on populateVariantAllele(data)  →  data.allele.variants_by_build
                └── fetch api/variant_sample_genotypes/<variant_id>            (one per *distinct* variant id)
```

The variant the page is on never waits on the ClinGen/allele round trip. The allele's other variants are requested as soon as the allele response lands, and each one's rows are appended to the grid on arrival.

**The unit of work is a variant, not a genome build.** A `Variant` is already tied to a build through `locus.contig`, so the build is derivable and does not belong in the request. This matters because *the mapping is not one variant per build*:

- Normal case — GRCh37 and GRCh38 have separate `Contig` rows (different RefSeq accessions), so the allele has a distinct variant per build. One call each.
- **Shared contigs** — builds can share contigs (chrM/MT, and some unplaced scaffolds); `Contig` exists precisely so this is possible (`snpdb/models/models_genome.py:378`), and `VariantAllele` documents it (`snpdb/models/models_variant.py:865`): *"Some builds share contigs (eg GRCh37/38 share MT and some unplaced scaffolds) - in those cases we'll have the same variant linked through different VariantAlleles"*. For an MT allele, `variant_for_build_optional(GRCh37)` and `variant_for_build_optional(GRCh38)` return **the same variant**.

Keying calls on build would fetch that variant twice and double every MT row in the grid. Keying on distinct variant id makes one call, and the response covers all builds sharing the contig.

### 2.1 API

New endpoint in `snpdb/views/views_rest.py`, registered in `snpdb/urls.py` alongside the existing `api/variant_allele_for_variant/...`:

```
GET /snpdb/api/variant_sample_genotypes/<int:variant_id>
```

No build in the URL — the server derives it from `variant.genome_builds` (`snpdb/models/models_variant.py:675`, which walks `GenomeBuildContig` from the variant's contig). That returns a *set*: one build normally, two for a shared contig like MT.

Sample scoping therefore uses `vcf__genome_build__in=variant.genome_builds`, and **each row's build comes from that row's sample's VCF**, not from the request. The shared-contig case then falls out for free: one call, rows from both builds' VCFs, each correctly labelled.

Responses:

- `200` — payload below.
- `404` — no such variant.

The client never asks for a build the allele has no variant in; it derives its call list from the allele's variants (§3.1), so "build not lifted over" is a client-side absence, not a server round trip.

```jsonc
{
  "genome_builds": ["GRCh37", "GRCh38"],
  "variant_id": 4567,
  "variant": "1:12345 A>G",
  "variant_url": "/variantopedia/view_variant/4567",
  "g_hgvs": "NC_000001.11:g.12345A>G",
  "samples": {"total": 1200, "visible": 900},
  "observations": {"total": 12, "visible": 10, "invisible": 2},
  "zygosity_counts": {"E": 8, "O": 2},
  "locus_counts": [
    {"variant_id": 4567, "variant": "1:12345 A>G", "url": "...", "description": "This variant",
     "total": 10, "REF": 0, "HET": 8, "HOM_ALT": 2, "Unknown": 0}
  ],
  "multiallelic": [
    {"multiallelic": "1:12345 A>G,T", "variants": [{"id": 4568, "label": "1:12345 A>T", "url": "..."}]}
  ],
  "patient_ids": [12, 15],
  "truncated": true,
  "rows": [
    {
      "genome_build": "GRCh38",          // from this sample's VCF, not from the request
      "sample": 88, "sample_name": "S1", "sample_url": "/snpdb/view_sample/88",
      "vcf": "run1.vcf",
      "zygosity": "E",
      "allele_frequency": 0.48, "allele_depth": 30, "read_depth": 62, "phred_likelihood": 99,
      "filters": "PASS",
      "sample_filters": null,
      "enrichment_kit": "idt_exome",
      "patient_hpo": "...", "patient_omim": "...", "patient_mondo": "...",
      "classifications": [
        {"id": 55, "url": "/classification/view/55", "label": "VUS A",
         "css_class": "cs cs-vus_a", "title": "Clinical Significance", "lab": "Example Lab"}
      ]
    }
  ]
}
```

Row keys are deliberately flattened/renamed (`sample_name` rather than `sample__name`) — the ORM `values()` names are an implementation detail and should not be the API contract. Rows are a flat list of independent objects, so #1693 can add keys without reshaping anything.

The endpoint takes an optional `?limit=` query parameter (§6). `truncated` reports whether `rows` was cut; every count in the payload is always exact.

### 2.2 Server-side structure

Keep everything in `snpdb/variant_sample_information.py`:

- **`VariantZygosityCounts`** — two changes:
  - Add `filters` / `samples_filters` and the VCF id to the `values()` in `_get_sample_values_for_variant_via_cohort_genotype` (line 96), and decode them in `_cohort_genotype_to_sample_genotypes` (line 106) using `VCFFilter.get_formatter(vcf)` built once per `CohortGenotypeCollection`.
  - ~~Drop the `genome_build` constructor parameter~~ — **already done** (§1), scoping now derives from `variant.genome_builds`. Nothing further needed here.
- **`VariantSampleGenotypes(VariantZygosityCounts)`** — replaces `VariantSampleInformation`. Produces the JSON payload:
  - locus counts as a list of dicts (drop pandas — `_get_locus_counts` at line 226 becomes a plain list builder; `has_unknown_zygosity` becomes a check on the counter rather than DataFrame columns),
  - the enrichment passes below, applied to `visible_rows` only,
  - the counts/observations summary, fixing the missing `has_hidden_samples`.
- **`_add_classifications(user, allele, rows)`** in the same module — the one enrichment pass, applied after row selection.

`snpdb/serializers.py` gains a `VariantSampleGenotypesSerializer` (or the view returns the dict directly with `@extend_schema(responses=OpenApiTypes.OBJECT)`, matching `VariantZygosityForSampleView`).

**Beacon and cross-build results.** The scoping fix does **not** make Beacon aggregate across builds via the allele. Shared contigs are not liftover: GRCh37 and GRCh38 share `NC_012920.1` (rCRS), so an MT variant is *one* variant record at *identical* coordinates, and samples from either build's VCFs are native observations of it. Beacon still counts exactly one variant.

Aggregating genuinely lifted-over observations (different variant records linked by an allele) is a separate question, and lifted-over counts should not be merged into a native count. If that is wanted, the natural shape is a **third resultSet** rather than a wider number — and there is precedent in this codebase: `classifications_dataset` is already allele-scoped (`beacon/datasets.py:141`, `classification__allele=allele`) while `observations_dataset` is variant-scoped, so Beacon already returns a cross-build set alongside a build-native one. A `variantgrid_observations_other_builds` dataset would follow the same pattern, keeping native and lifted-over counts separately attributable. Out of scope here — noted so the option isn't lost.

**Import direction check**: `snpdb/variant_sample_information.py` will import from `classification.models`. Several `snpdb` non-model modules already do this at top level (`snpdb/common_variants.py:7`, `snpdb/grid_columns/custom_columns.py:9`), and nothing in `snpdb.models` imports `variant_sample_information`, so a top-level import is correct here. Verify with a `manage.py check` after wiring it up.

### 2.3 Client-side structure

Replace the jqGrid with a **client-side DataTables** instance (`data: []`, rows appended with `rows.add(...).draw()`):

- DataTables is the current preferred grid (per CLAUDE.md) and the bundle already includes Buttons + HTML5 CSV export, which replaces the hand-written `exportGenotypeGridCSV`.
- `render` callbacks give us clean HTML cells for classification pills, sample links and filters.
- Incremental `rows.add()` is the natural fit for streaming one response at a time; the jqGrid local-data path would need re-initialising on each.
- The zygosity checkbox filter becomes a `$.fn.dataTable.ext.search` predicate over the accumulated data (currently it uses jqGrid `postData.filters`).

Columns:

| Column | Notes |
|---|---|
| Genome Build | Hidden until a second build returns rows |
| Sample Name | link, as now |
| Zygosity | as now |
| **FILTER** | new — §3.2 |
| **Classification** | new — §3.3 |
| Enrichment Kit, Patient HPO/OMIM/MONDO, VCF, AF, AD, DP, PL | as now |

Move the fragment's JavaScript into a static file (e.g. `variantgrid/static_files/default_static/js/variant_sample_information.js`) with the template supplying only the URLs and settings — the current fragment mixes ~200 lines of JS with Django template variables.

### 2.4 Per-build loader

Above the grid, a row of build badges:

```
Genome builds:  [ GRCh38  ✔ 10 samples ]   [ GRCh37  ⏳ loading… ]   [ T2T-CHM13  — no variant in this build, not searched ]
```

Badges are per **build** (that is what users think in), while requests are per **variant**. The client keeps a `build → variant_id` map from the allele and a `variant_id → request` map, so:

- one build per variant — badge tracks its own request;
- **shared contig** — GRCh37 and GRCh38 point at the same variant id, so both badges track the *same* in-flight request and resolve together. Render these merged, e.g. `[ GRCh37 + GRCh38  ✔ 10 samples — same variant (MT) ]`, so it is obvious the counts are not being double-counted;
- no variant for a build — no request is made at all; the badge shows "not searched" from the outset.

States: `pending` (spinner) → `loaded` (visible observation count) / `empty` (0 samples) / `missing` (no variant in that build) / `error`. Each badge updates as its request resolves.

The badge row is the answer to "what was actually searched" — a build with no variant is visibly *not searched* rather than silently absent.

The locus-counts table and the multi-allelic section are also rendered client-side, one block per response, labelled with that response's build(s), from `locus_counts` / `multiallelic`. A shared-contig response yields a single block covering both builds — the locus is shared too, so there is genuinely only one set of counts.

---

## 3. Feature detail

### 3.1 Other genome builds

- Build discovery needs per-build **variant ids**, not just build names, so the client can collapse duplicates. `AlleleSerializer` (`snpdb/serializers.py:69`) has an explicit `fields` tuple, so add a `variants_by_build` `SerializerMethodField` built from `Allele.variant_alleles()` — additive, `build_names` stays as-is for existing consumers:

  ```json
  "variants_by_build": [{"genome_build": "GRCh37", "variant_id": 4567},
                        {"genome_build": "GRCh38", "variant_id": 4567}]
  ```

  The client requests `new Set(variants_by_build.map(v => v.variant_id))` minus the variant already loaded. For an MT allele that set has **one** element, so the grid gets one copy of each row.
- `variant_details.html` currently calls `populateVariantAllele(data)` in two places (inline data, and the ajax success path). Add a single hook there that forwards the allele data to the sample-information module if it has loaded (`window.variantSampleInformation?.addBuilds(data)`), and have the module handle being called before/after its own initialisation.
- Builds with no variant for the allele **return 404 and the client stops there** — no liftover is offered or triggered. A variant that has just been created by liftover cannot have samples attached to it, so a liftover button would do nothing useful in a samples panel. Liftover belongs on the allele page; the badge simply records that the build was not searched.
- The *locus* scan is per variant: `VariantZygosityCounts` walks `variant.locus`, so locus counts are scoped to whatever builds that locus belongs to — one build normally, both for a shared contig.

### 3.2 VCF FILTER column

`CohortGenotype` has two relevant fields (`snpdb/models/models_cohort.py:738`, `:747`):

- `filters` — record-level FILTER, stored as a string of per-VCF single-char codes (`VCFFilter.filter_code`); `null`/empty means PASS.
- `samples_filters` — per-sample FT array, `MISSING_FT_VALUE = "NULL"` when absent.

Both must be decoded with the **owning VCF's** `VCFFilter` map — codes are assigned per VCF (`upload/vcf/vcf_import.py:100`), so a shared lookup across VCFs would be wrong. `_cohort_genotype_to_sample_genotypes` already loads the `CohortGenotypeCollection` per row; build the formatter there once per collection (`VCFFilter.get_formatter(cgc.cohort.vcf)`).

Emit `filters` (record) always and `sample_filters` (FT) only when the VCF stored them. Render PASS in muted text and non-PASS values highlighted, so filtered calls stand out.

### 3.3 Classifications per sample

For the allele of the variant:

```python
qs = ClassificationModification.latest_for_user(
    user, allele=allele, published=True,
    classification__sample__in=visible_sample_ids,
)
```

`latest_for_user(allele=...)` filters on `classification__variant__in=allele.variants`, so it naturally covers classifications made against *any* build's variant — which is exactly right for a multi-build grid.

Per classification, take `classification.summary_typed["pathogenicity"]["classification"]` (values like `B`, `LB`, `VUS_A`, `LP`, `P`) and derive label + CSS class the same way `clinical_significance_values` does (`classification/templatetags/classification_tags.py:210`), including the pending-change case. Include the somatic tier pill when `allele_origin_bucket != "G"` or a somatic value is set. Do **not** use `Classification.clinical_significance` — the model comments flag it as an optimisation field that is "relatively out of date".

Group by `classification.sample_id` into a `{sample_id: [pill, ...]}` map and attach to rows. Render client-side as `<span class="c-pill {css_class}">{label}</span>` linked to the classification, reusing the shape of `classification_clinical_significance_renderer`.

Classifications on this allele that are **not** linked to a sample are already shown by the `{% classification_groupings %}` block higher up the page (`variant_details.html:559`) — no duplication needed here.

### 3.4 Phenotype/OMIM graphs

These are server-rendered from `vsi.patient_ids` via `{% patient_phenotypes_graph %}` / `{% patient_omim_graph %}` (`patients/templatetags/patient_graph_tags.py:47`), which cannot work from the async payload as-is.

Each build response returns `patient_ids`. The client accumulates them and, once all build requests have settled, loads a small fragment endpoint that renders both graphs for those patients. `match_graph` already intersects with `Patient.filter_for_user`, so accepting patient ids from the client is safe.

---

## 4. Caching

Apply the established pattern from `snpdb/views/views_rest.py:42`:

```python
@method_decorator([cache_page(MINUTE_SECS), vary_on_cookie], name='dispatch')
```

`vary_on_cookie` keys on the session, so permission-scoped payloads are not shared between users. Keep the TTL short (a minute, as with `TrioView`/`QuadView`) — a longer TTL would hide newly imported VCFs and newly published classifications from a page users refresh precisely to see those.

---

## 5. Permissions

Nothing changes about how visibility is decided, and no new decorators are needed (global login + CSRF middleware, DRF `IsAuthenticated` by default — see CLAUDE.md *Security*). The invariants to preserve and test:

- Rows are limited to `Sample.filter_for_user`; hidden observations are only ever reported as a count.
- Classifications come through `ClassificationModification.latest_for_user` (read permission on the modification, not the classification).
- Resolving the allele's other variants goes through the allele; each variant's response is permission-scoped exactly as the page's own variant is.
- The row cap is applied after permission filtering, so it can never widen visibility.

---

## 6. Row limit

This view is for rare disease work — there is no value in spending database and browser resources materialising every row for a common variant. So **rows are capped by default here, while counts stay exact everywhere.**

- New setting `VARIANT_DETAILS_SAMPLES_MAX_ROWS` (default 1000) in `variantgrid/settings/components/default_settings.py`.
- The endpoint accepts `?limit=<n>`, defaulting to the setting. `limit=0` means no cap.
- **Counts are never capped.** `num_observations`, `num_visible_observations`, the zygosity counter and the locus counts all come from the full locus walk in `VariantZygosityCounts`, which runs regardless. Only the `rows` list is sliced, and `"truncated": true` says so.
- **Beacon is unaffected.** `beacon/datasets.py` uses `VariantZygosityCounts` directly and never touches rows, so its counts stay exact with no change and no opt-out needed.
- **"Load all N rows" button** per truncated build, which re-requests that build with `limit=0` and replaces its rows in the grid. The badge for a truncated build shows e.g. `GRCh38 ✔ 1,000 of 4,312 shown`.
- **CSV export gets everything.** DataTables' HTML5 CSV button exports what is loaded, so when any build is truncated the export handler first loads the full set for those builds (with a spinner), then exports. A download must not silently be a truncated view of the data.

Applying the cap *after* the permission filter (i.e. to `visible_rows`) keeps it a pure presentation limit — it can never change which samples a user is allowed to see, only how many are drawn at once.

## 7. Settings

- `VARIANT_DETAILS_SHOW_SAMPLES` (existing) still gates the whole section.
- `VARIANT_DETAILS_SAMPLES_MAX_ROWS` (new) — see above.

---

## 8. Testing

- `variantopedia/tests/test_urls.py:81` — keep the `variant_sample_information` shell URL test (still a 200), add the new API URL.
- New `snpdb/tests/test_variant_sample_genotypes.py`:
  - rows restricted to visible samples; invisible observation counts correct,
  - FILTER decoding: PASS, a single filter, multiple filters, and two VCFs with **different** code→id maps in one response,
  - **shared contig (MT)**: a variant whose contig is in both builds returns rows from both builds' VCFs in one response, each row labelled with its sample's VCF build; `num_invisible_observations` is **0** when the user can read all of them (the bug in §1);
  - the allele payload's `variants_by_build` collapses to a single variant id for an MT allele, so the client issues one request;
  - classifications attached to the right sample, and a classification the user cannot read is absent,
  - row cap: `truncated` set, `rows` sliced, **and every count still exact**; `limit=0` returns everything.
- Shared-contig scoping is already covered by `snpdb/tests/test_variant_sample_information.py` (shared contig counts both builds; build-specific contig does not; `has_hidden_samples`). Keep these passing.
- Counts must not depend on the row cap — assert this for Beacon and for the endpoint.
- Fixtures: `snpdb/tests/utils/` has VCF/cohort helpers; classifications-linked-to-samples may need a small addition.

---

## 9. Work breakdown

Each phase leaves the page working.

0. ~~**Shared-contig scoping fix.**~~ **Done** — see §1.
1. **API + client-side DataTables port, single variant only.** New REST endpoint returning the payload for one variant; shell fragment; DataTables replacing jqGrid (zygosity filter, CSV via Buttons, Plotly plots fed from accumulated rows); locus counts and multi-allelic rendered client-side; drop pandas from the module; fix `has_hidden_samples`. Adds the FILTER column (§3.2) since it is pure payload work.
2. **Row cap.** Setting, `?limit=`, `truncated`, "load all" button, full-data CSV export. Early, so the multi-build phase never multiplies an unbounded payload.
3. **Multi-build.** `variants_by_build` on `AlleleSerializer`, distinct-variant-id request set, build badges/loader (including the merged shared-contig badge), Genome Build column, per-build locus tables, remove the "current genome build only" help text.
4. **Classifications** (§3.3).
5. **Phenotype graphs fragment, caching, tests.**

Sample membership and trio inheritance are tracked separately in #1693.

---

## 10. Decisions

| # | Decision |
|---|---|
| 1 | **Port to DataTables.** Confirmed — desirable generally and here in particular. |
| 2 | **No variant in a build → stop.** No liftover offered from the samples panel; a freshly lifted-over variant cannot have samples. Liftover stays on the allele page. The build badge shows it was not searched. |
| 3 | **Cap rows in this view by default**, exact counts always, Beacon unaffected *by the cap*, plus a "load all" button and full-data CSV export. |
| 4 | **Trio inheritance and sample membership split into their own issue** — raised as #1693. Different feature, own design questions (confirmed vs putative de novo, multi-trio roles, the `Zygosity` refactor). Row objects take the new keys additively when it lands. |
| 5 | **Endpoint is keyed on variant id, not (variant, build).** A variant already knows its build via `locus.contig`, and the mapping is not 1:1 — builds sharing a contig (MT, unplaced scaffolds) resolve to the *same* variant. Requests are per distinct variant id so shared-contig rows are fetched once; build labels come from each row's sample's VCF. The underlying scoping fix is already implemented (§1). Beacon still counts one variant — shared contigs are identical coordinates, not liftover; genuine cross-build aggregation would be a separate resultSet (§2.2). |
