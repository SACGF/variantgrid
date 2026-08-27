# #1554 — Pfam protein domains from the InterPro API, fetched lazily per gene

[#1554](https://github.com/SACGF/variantgrid/issues/1554). The bulk ingest of protein domains died
when EBI stopped publishing per-organism Pfam files; this replaces it with a per-gene fetch off the
InterPro REST API, done on demand and cached in the existing tables.

## The problem

`store_pfam_from_web()` (`genes/cached_web_resource/pfam.py:29`) raises `NotImplementedError` and Pfam
sits in `DISABLED_CACHED_WEB_RESOURCES` (`variantgrid/settings/components/annotation_settings.py:438`,
surfaced by `annotation/views.py:283`). It was disabled because
`proteomes/9606.tsv.gz` no longer exists — confirmed, `current_release/` is now Pfam38.2 with no
`proteomes/` subdirectory at all, and per-organism breakdowns are gone for good.

Two of the three steps still work. `store_pfam()` (`:41`) reads `Pfam-A.clans.tsv.gz`, still published.
The UniProt half of `store_pfam_sequences_and_domains()` (`:63`) reads
`HUMAN_9606_idmapping.dat.gz`, still published, and is what fills `PfamSequence` and
`PfamSequenceIdentifier` — the transcript → UniProt mapping everything else hangs off. Only the last
step, `insert_domains()` (`:131`), lost its source.

Everything downstream reads through one method, `TranscriptVersion.protein_domains_and_accession`
(`genes/models.py:1257`), and only as the fallback — a `ProteinDomainTranscriptVersion` overrides it.
Three call sites: the hotspot graph's transcript search (`genes/views/views_hotspot_graphs.py:129`),
the same view's context build (`:174`), and the transcript autocomplete's `has_protein_domains` filter
(`genes/views/views_autocomplete.py:106`). The classification `pfam_protein_domain` evidence key comes
from VEP's `DOMAINS` field (`classification/autopopulate_evidence_keys/evidence_from_variant.py:157`)
and never touches these tables.

## What the numbers say

Measured against the live API and this deployment's data:

| | |
|---|---|
| `GET /interpro/api/entry/pfam/protein/uniprot/<acc>/` | mean **1.60s**, median 1.62s, min 1.13s, max 2.40s (50 serial calls) |
| proteins with no Pfam match | **11 of 50** — HTTP **204**, empty body |
| domains per protein that has any | ~2.7 |
| `PfamSequence` rows | 90,149 |
| distinct `PfamSequence` per gene symbol | mean 4.2, **median 3**, p90 9, p99 22, max 84 (GBF1) |

So a full prefetch of all 90,149 is 40 hours serial and not worth doing, while a single gene is a
median of 3 calls. That is the whole argument for going lazy.

Two alternatives were measured and rejected: `Pfam-A.regions.tsv.gz` is 5.05 GB / ~142M rows and takes
16–85 min to stream depending on EBI's mood (HTTPS at `https://ftp.ebi.ac.uk/...` runs 1.0–5.2 MB/s,
FTP only 0.8 MB/s); the bulk-by-Pfam-entry API endpoint is 28.6s per entry cold × 9,155 human entries
= 73 hours. Both pay for 90k proteins to serve the few thousand anyone looks at.

## The design

Fetch **all of a gene's sequences in one batch**, not one per transcript version. The hotspot graph's
`SEARCH_ORDER` loop walks every candidate transcript version until one has domains, and the transcript
autocomplete walks every transcript for the symbol — batching by gene means both are satisfied by a
single warm-up, and neither needs to know a fetch ever happened.

### 1. `PfamSequence.domains_imported`

```python
class PfamSequence(models.Model):
    seq_id = models.TextField(primary_key=True)
    domains_imported = models.DateTimeField(null=True)
```

11 of 50 proteins legitimately have zero domains, so "no `PfamDomains` rows" cannot mean "not fetched
yet" — without this stamp every view of a domain-less protein re-hits InterPro forever. It is also the
staleness handle if a refresh policy is ever wanted.

### 2. `genes/interpro.py` — the fetcher

Named for the service, alongside `genes/panel_app.py` and `genes/cdot_data_release.py`.

```python
def store_domains_for_sequences(seq_ids: Iterable[str]) -> int
def store_domains_for_transcripts(transcript_ids: Iterable[str]) -> int
```

`store_domains_for_transcripts` resolves transcripts → `PfamSequence` via `PfamSequenceIdentifier`,
drops any already stamped, and hands the rest to `store_domains_for_sequences`. That signature suits
every caller — the hotspot view has a transcript version list, the autocomplete has a transcript
queryset.

Per sequence: `GET {PFAM_INTERPRO_API_URL}/entry/pfam/protein/uniprot/{seq_id}/?page_size=200`,
following `next` (`count` is the number of matching Pfam entries and the default page size is 20, so a
protein matching many families would otherwise be truncated).

Response shape:

```json
{"count": 3, "next": null, "results": [
  {"metadata": {"accession": "PF00097", "name": "Zinc finger, C3HC4 type (RING finger)"},
   "proteins": [{"accession": "p38398", "protein_length": 1863,
                 "entry_protein_locations": [{"fragments": [{"start": 24, "end": 64,
                                                             "dc-status": "CONTINUOUS"}]}]}]}]}
```

One `PfamDomains` row **per fragment** — discontinuous domains then draw as separate boxes, which the
model already allows ("No constraints as there can be multiple of same domains within a sequence",
`genes/models.py:2561`). Map `metadata.accession` through `Pfam.get_pk_from_accession()`; skip and
count anything not in `Pfam`, as `insert_domains()` did.

Writes per sequence are delete-then-`bulk_create`-then-stamp, in a transaction, so a sequence that
already holds stale rows converges instead of doubling up.

Outcomes:

- **200 with results** → rows created, `domains_imported` stamped.
- **204** → no Pfam matches. Stamped with zero rows. This is the common case that makes the stamp
  necessary.
- **anything else, timeout, connection error** → logged, **left unstamped**, retried on the next view.
  The caller never raises; a graph renders without a domain track rather than 500ing.

Concurrency: a `ThreadPoolExecutor` doing the HTTP, results applied to the DB on the main thread — no
Django connections in workers. At `PFAM_INTERPRO_MAX_WORKERS = 8` the median gene is under 2s, p90
~2s, p99 ~4.5s. The tail is the problem: GBF1's 84 sequences would be ~17s. So the batch also carries
`PFAM_INTERPRO_DEADLINE_SECONDS = 20` — collect via `as_completed`, stop draining at the deadline,
keep everything already returned, and leave the rest unstamped for the next view to pick up. Page
latency is bounded and no work is thrown away. (If the tail proves annoying in practice, the same
function moves behind a Celery task on `web_workers` and the first view renders domain-less — worth
knowing, not worth building yet.)

HTTP goes through a module-level `requests.Session` mounted with the same `Retry` config as
`library/utils/misc_utils.py:137`, shared across the pool (GET-only, no per-request session state).

Guarded by `settings.PFAM_INTERPRO_LAZY_DOMAINS` and a `settings.UNIT_TEST` short-circuit, so tests
never reach the network by accident.

### 3. Where it is called from

**Hotspot graph** — `genes/views/views_hotspot_graphs.py`, immediately after
`tv_list = sorted(tv_qs, ...)` (`:108`) and before the `SEARCH_ORDER` walk (`:113`). `tv_list` is the
complete candidate set for all three lookup routes (transcript / gene / gene symbol), so one call
there covers `_pick_transcripts` and the `get_context_data` re-read at `:174`.

**Transcript autocomplete** — `genes/views/views_autocomplete.py`, at the top of the
`has_protein_domains` block (`:101`), over the same `qs` it is about to iterate. This one **does**
fetch rather than read-only. It is always scoped to a chosen gene symbol (`GeneAndTranscriptForm`,
`genes/forms.py:43`, forwards `gene_symbol`; the only user is `cohort_hotspot`,
`snpdb/views/views.py:1642`), so it is the same bounded per-gene batch — and the view is wrapped in
`cache_page(WEEK_SECS)` (`:79`), so a partial first answer would be served for a week. Better to pay
once and cache the complete list.

### 4. The bulk task keeps its working half

- `store_pfam_from_web()` (`:29`) — restore the original body minus the domains step; description
  becomes `f"{num_pfam} Pfam. {num_sequences} sequences."`.
- `store_pfam()` (`:41`) — unchanged.
- `store_pfam_sequences_and_domains()` (`:63`) → `store_pfam_sequences()`, returning the sequence
  count; drop the FTP fetch of `9606.tsv.gz` and the `insert_domains()` call.
- `insert_domains()` (`:131`) — deleted, superseded by `genes/interpro.py`.
- Module docstring loses the `proteomes/9606.tsv.gz` line and gains the InterPro API.
- Remove `CACHED_WEB_RESOURCE_PFAM` from `DISABLED_CACHED_WEB_RESOURCES`.

Note that `store_pfam_from_web()` opens with `Pfam.objects.all().delete()` and
`PfamSequence.objects.all().delete()`, so "Update from web" cascades away every lazily-fetched domain.
That is the right behaviour — a new UniProt/Pfam release invalidates them anyway and they refill on
demand — but it should be deliberate, and worth a comment at the delete so nobody reads it as a bug.

### 5. Existing deployments

A data migration stamping `domains_imported` on every `PfamSequence` that already has `PfamDomains`
rows. Boxes that ran the ingest successfully before the upstream file vanished keep their data and
skip a pointless 90k-call re-fetch; boxes with nothing (this dev DB is `PfamDomains` 0) start clean and
fill in as pages are viewed. No `ManualOperation` — nothing needs running by hand.

### 6. Settings

In `variantgrid/settings/components/default_settings.py`, near `PANEL_APP_CACHE_DAYS` (`:477`):

```python
PFAM_INTERPRO_LAZY_DOMAINS = True
PFAM_INTERPRO_API_URL = "https://www.ebi.ac.uk/interpro/api"
PFAM_INTERPRO_TIMEOUT_SECONDS = 30
PFAM_INTERPRO_MAX_WORKERS = 8
PFAM_INTERPRO_DEADLINE_SECONDS = 20
```

## Suggested order

1. `PfamSequence.domains_imported` + migration, and the data migration stamping existing rows (§1, §5).
2. `genes/interpro.py` with its tests against mocked responses (§2).
3. Wire up the hotspot graph, then the autocomplete (§3).
4. Re-enable the cached web resource and drop it from `DISABLED_CACHED_WEB_RESOURCES` (§4).
5. Run "Update from web" on a dev box, then load a gene page and confirm the domain track draws.

## Tests

Worth keeping, all against mocked InterPro responses:

- a 204 stamps `domains_imported` and creates no rows — the case the whole design turns on
- an already-stamped sequence makes no HTTP call
- multi-fragment `entry_protein_locations` produce one row per fragment
- a `metadata.accession` with no `Pfam` row is skipped and counted, not an error
- a failed request leaves `domains_imported` null and returns cleanly
- re-fetching a sequence that already has rows replaces rather than duplicates them

Not worth keeping: anything asserting the field is nullable, or that the thread pool runs things in
parallel.

## Done when

A gene symbol page draws its protein domain track on first view without a bulk import, the transcript
autocomplete's "has protein domains" filter returns the same transcripts it did before the upstream
file disappeared, and Pfam no longer shows as disabled on the annotation admin page.
