# #1690 Replace search_signal with an explicit search registry

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-10-09
Status: in progress

[#1690](https://github.com/SACGF/variantgrid/issues/1690). Search receivers are today a Django `Signal` used as a
plugin registry: `snpdb/search.py:search_receiver` wraps the function and calls `search_signal.connect()` at import
time, each app imports its `*_search` modules from `AppConfig.ready()` purely for that side effect, and
`snpdb/search.py:SearchInput.search` fires `send_robust()` and sifts `(receiver, result)` tuples. The thread on the
issue settled the shape: the decorator builds a receiver object, and each app's `ready()` registers its receivers
**explicitly by name**, so the imports are real uses that no unused-import autofix can strip (the failure mode the
issue is worried about: a dropped import silently removes a search, nothing fails, nothing logs).

What stays the same: every receiver module's body, the `@search_receiver(...)` arguments, `SearchResponse` and
everything downstream of it (`SearchResponsesCombined`, the search page template, the tests that call
`search_data`). What changes: the decorator's return value, how receivers are registered, how `SearchInput.search`
dispatches, and the ready() of nine apps.

Inventory at the time of writing (52 receivers in 28 modules; the issue said 45, more have landed since):

| App | Module | Receivers |
|---|---|---|
| analysis | `analysis/signals/analysis_search.py` | search_analysis |
| annotation | `annotation/signals/citation_search.py` | search_citations |
| classification | `classification/signals/classification_search.py` | classification_search |
| classification | `classification/signals/discordance_report_search.py` | discordance_report_search |
| genes | `genes/signals/gene_search.py` | gene_search, gene_version_search |
| genes | `genes/signals/gene_symbol_search.py` | gene_symbol_alias_search |
| genes | `genes/signals/transcript_search.py` | search_transcript |
| ontology | `ontology/models/ontology_search.py` | ontology_search_id, ontology_search_hgnc, omim_name_search, mondo_name_search, hpo_name_search |
| patients | `patients/signals/external_pk_search.py` | search_external_pk |
| patients | `patients/signals/patient_search.py` | patient_search |
| patients | `patients/signals/specimen_search.py` | specimen_search, extraction_search |
| pedigree | `pedigree/signals/pedigree_search.py` | search_pedigree |
| seqauto | `seqauto/signals/enrichment_kit_search.py` | enrichment_kit_search |
| seqauto | `seqauto/signals/experiment_search.py` | experiment_search |
| seqauto | `seqauto/signals/sequencing_run_search.py` | sequencing_run_search |
| seqauto | `seqauto/signals/tso500_pair_search.py` | tso500_pair_search |
| snpdb | `snpdb/signals/clinvar_export_search.py` | clinvar_id_search, clinvar_export_batch_search |
| snpdb | `snpdb/signals/cohort_search.py` | search_cohort |
| snpdb | `snpdb/signals/duo_search.py` | search_duo |
| snpdb | `snpdb/signals/genomics_search.py` | genome_build_search, contig_search |
| snpdb | `snpdb/signals/lab_search.py` | lab_search |
| snpdb | `snpdb/signals/organization_search.py` | organization_search |
| snpdb | `snpdb/signals/quad_search.py` | search_quad |
| snpdb | `snpdb/signals/sample_search.py` | sample_search |
| snpdb | `snpdb/signals/scv_search.py` | scv_search |
| snpdb | `snpdb/signals/trio_search.py` | search_trio |
| snpdb | `snpdb/signals/user_search.py` | user_search |
| snpdb | `snpdb/signals/variant_search.py` | variant_cosmic_search, search_variant_locus_no_ref, search_variant_locus_with_ref, allele_search, variant_search_vcf, search_variant_gnomad, search_variant_variant, search_variant_symbolic, search_variant_db_snp, search_hgvs, search_variant_id, search_allele_id, search_variant_gene_fusion, search_variant_gene_copy_number, search_variant_splice_event |
| snpdb | `snpdb/signals/vcf_search.py` | vcf_search |

The implementer regenerates this list rather than trusting it (the AST scan in §4 is the same scan the new test
runs). `sapath` (the sibling repo symlinked into this checkout) has no search receivers; the Shariant repo is not on
this box and is checked after merge (§6).

## Data

Nothing in the database changes. The in-memory holders, in `snpdb/search.py`:

```python
@dataclass(frozen=True)
class SearchReceiver:
    """ One search, as declared by @search_receiver. Registered from the owning app's AppConfig.ready() """
    func: Callable[[SearchInputInstance], Iterable]   # the decorated generator, untouched
    search_type: PreviewCoordinator                   # required - today's Optional is never None in practice
    pattern: Pattern
    admin_only: bool
    sub_name: Optional[str]
    example: Optional[SearchExample]
    match_strength: Optional[SearchResultMatchStrength]
    enabled: bool                                     # settings-driven, evaluated once at import as today


class SearchRegistry:
    _receivers: list[SearchReceiver]    # registration order; responses are sorted later, so order carries no meaning


search_registry = SearchRegistry()
```

`SearchResponse` gains one field, filled in by the dispatcher:

```python
    duration_seconds: float = 0.0
    """ Wall time the receiver took, 0 when the pattern did not match. The timing the concurrency follow-on needs """
```

## Behaviour

### 1. `snpdb/search.py`

- `search_signal` and the `django.dispatch.Signal` import go. `send_robust` and the `(caller, response)` sifting in
  `SearchInput.search` go with them.
- `@search_receiver(...)` keeps its signature and its docstring (the "how to write a search" guide), and returns a
  `SearchReceiver` instead of a signal-connected closure. It registers nothing. `enabled=False` now yields a receiver
  with `enabled=False` rather than `None`, so the decorated name is always a `SearchReceiver`.
- `SearchReceiver.visible_to(search_input) -> bool`: the four gates that today sit at the top of the wrapper, in
  the same order - `enabled`, `admin_only` against `user.is_superuser`, `search_type.preview_enabled()`, and the
  classify gate (`search_input.classify` only admits receivers whose `preview_category()` is `"Variant"`).
- `SearchReceiver.search(search_input) -> SearchResponse`: the rest of today's wrapper body verbatim - pattern match,
  `SearchInputInstance`, the result cap (`MAX_VARIANT_RESULTS` for Variant, else `MAX_RESULTS_PER_SEARCH`),
  `INVALID_INPUT`, the `_SearchResultFactory` loop, the default `match_strength`, exception capture into a
  `SearchMessageOverall` with `report_exc_info()`, the `SearchResponse` construction and the
  `settings.PREFER_ALLELE_LINKS` allele conversion. Two small fixes while it moves: the `"returned None"` error
  names `self.func.__name__` (today it reads `sender.__name__`, which is always `SearchInput`), and the receiver
  times the call with `time.perf_counter()` into `duration_seconds`. When the pattern does not match the response
  is built as today with `matched_pattern=False` - that response is what the search page's "Accepted Inputs" card
  lists, so it stays.
- `SearchRegistry.register(*receivers)`: appends each; a receiver already registered (same object) raises
  `ValueError` naming the function, since Django runs each `ready()` once and a second registration means two apps
  claim one search. `SearchRegistry.receivers -> tuple[SearchReceiver, ...]`.
- `SearchInput.search`:
  ```python
  return [receiver.search(self) for receiver in search_registry.receivers if receiver.visible_to(self)]
  ```
  No type-checking of results and no "doesn't happen" branch: `SearchReceiver.search` owns error capture. The
  responses are sorted by `SearchResponse.__lt__` in `SearchResponsesCombined` as before.
- `SearchResponsesCombined.summary` (the string the search view records in the EventLog, parsed by nothing) appends
  the slowest receiver when any took over a second, e.g. `slowest: Variant/HGVS 1.3s`. This is the data the
  issue's concurrency follow-on wants before deciding anything.
- Module docstring: rewritten to say receivers are `SearchReceiver` objects registered from each app's
  `AppConfig.ready()` via `search_registry.register`, and that `SearchInput.search` runs the visible ones.

### 2. Each app's `ready()`

Nine apps register their receivers by name. The existing ready() convention (imports inside the method under
`# pylint: disable=import-outside-toplevel`, because `apps.py` cannot import models at module level) stays; the
`# noqa: F401` / "Registers receivers on import" comments on the `*_search` imports go, because the names are now
used. Non-search modules in the same import lists (health checks, previews, hooks) keep their `# noqa: F401`.

```python
# snpdb/apps.py (the others follow the same shape)
from snpdb.search import search_registry
from snpdb.signals.variant_search import (
    allele_search, search_allele_id, search_hgvs, ...
)
...
search_registry.register(
    allele_search, search_allele_id, search_hgvs, ...,
    search_cohort, search_duo, ...,
)
```

- `classification/apps.py` currently does `import classification.signals` (whose `__init__` star-imports every
  hooks module). That import stays for the hooks; the two search receivers are imported by name from
  `classification/signals/classification_search.py` and `classification/signals/discordance_report_search.py`
  and registered.
- `ontology/apps.py` registers the five receivers from `ontology/models/ontology_search.py` (the module lives
  under `models/` and is already imported through `ontology/models/__init__.py`; that stays as is).
- `genes`, `annotation`, `analysis`, `patients`, `pedigree`, `seqauto`, `snpdb`: register the receivers of the
  modules their ready() already imports.

Registration order is INSTALLED_APPS order, and it does not matter: results are sorted.

### 3. Call sites

`variantopedia/views.py:search` keeps calling `search_data` with an empty string to list the accepted inputs; the
comment above it is rewritten to say what now happens (every visible receiver returns an unmatched response; no
pattern matches the empty string). `analysis/variant_text.py`, `classification/utils/clinvar_matcher.py` and the
tests that call `search_data` or build a `SearchInput` need no change.

### 4. Tests

- New test module snpdb/tests/test_search_registry.py:
  - **Every declared receiver is registered.** AST-scan every installed app (reuse the walk in
    `library/tests/test_signal_receiver_registration.py`, skipping the same directories) for module-level functions
    decorated `@search_receiver`, and assert each `(module, name)` pair is in `search_registry.receivers` by
    `func.__module__` / `func.__name__`. The failure message says to add the name to that app's
    `AppConfig.ready()` `search_registry.register(...)` call. Guard the scan as the existing test does (more than
    40 found, `snpdb.signals.trio_search.search_trio` among them).
  - **Dispatch gates**, with a stub receiver (a `PreviewProxyModel` as `search_type`, a generator that yields one
    `PreviewData`) run through `SearchReceiver.visible_to` and `SearchReceiver.search` directly, no registry and no
    database: `enabled=False` is not visible; `admin_only` is visible only to a superuser; `classify=True` hides a
    non-Variant receiver; a non-matching pattern gives `matched_pattern=False` with no results; an exception in
    the generator becomes an `ERROR` `SearchMessageOverall` and the response still comes back. `register` twice
    with the same receiver raises.
- `library/tests/test_signal_receiver_registration.py`: drop `search_receiver` from `REGISTRATION_DECORATORS` and
  its mention in the docstring, since the decorator no longer registers on import. Modules that also carry
  `@receiver` handlers stay covered by it; the pure search modules are covered by the new test.
- Existing search tests (`variantopedia/tests/test_search.py`, `variantopedia/tests/test_gene_level_search.py`,
  `snpdb/tests/test_variant_search.py`, `patients/tests/test_specimen_extraction.py`,
  `patients/tests/test_sample_preview.py`, `analysis/tests/test_intersection_node_variant_text.py`) pass unchanged
  and are the behavioural regression check; `scripts/vg tests --explain` lists the rest.

### 5. Docs

- `snpdb/AGENTS.md`: the "Search handlers register with search.py:search_receiver" line (the one that says every
  other receiver is connected in ready, not at import) becomes: every receiver, search included, is connected in
  `snpdb/apps.py:SnpdbConfig.ready`; a search is declared with `snpdb/search.py:search_receiver` and registered
  there with `search_registry.register`, and the registry test fails when one is declared and not registered.
- `claude/research/snpdb.md`: the two paragraphs that describe `ready()` importing signals modules "for its
  `@search_receiver` side effect" and `SearchInput.search` sending `search_signal` are updated to the registry.
- `scripts/vg docs check` passes afterwards.

### 6. After merge

The Shariant repo (`SACGF/variantgrid_shariant`) is grepped for `@search_receiver`; any it declares needs a
`search_registry.register(...)` in its `ready()` or those searches silently stop. The `sapath` sibling has none.

## Out of scope

Running receivers concurrently (the issue's follow-on). The `duration_seconds` field and the EventLog summary are
the measurement that decision is waiting on.
