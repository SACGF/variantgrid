# MME Gene Identifiers Plan

Follow-on to `claude/matchmaker_exchange_plan.md`. Covers how VariantGrid identifies
genes on the MatchMaker Exchange wire, in both directions.

## 1. Why

Two reasons, one urgent and one scheduled.

**Urgent — we are silently losing matches today.** `mme/matching.py:_genes_from_profile()`
reduces each `genomicFeatures[].gene.id` to an upper-cased string and intersects the two
sets. We always emit an HGNC **symbol** (`mme/serializers/patient_profile.py:96`,
`feature["gene"] = {"id": str(gene_symbol)}`). A peer that sends
`{"gene": {"id": "ENSG00000133392"}}` for the same gene we hold as `MYH11` produces an
empty intersection, so the gene term contributes **0.0** to the score. Gene carries
`_GENE_WEIGHT = 0.5` — the largest single component — so the strongest evidence type is
the one most likely to silently evaluate to zero. No error, no log line; the match just
scores low and drops below `MIN_SCORE`.

**Scheduled — v1.0 → v1.1 changed the guidance and v2.0 changes the rule.** v1.1 reordered
the accepted `gene.id` forms to put Ensembl first and states that Ensembl gene IDs are
"***strongly*** encouraged and **will become mandatory in 2.0**". Symbols are the form we
currently emit and the form that stops being legal.

Related: MME metrics `numberOfUniqueGenes` / `numberOfUniqueGenesMatched`
(`mme/metrics.py`) count upper-cased symbols, so those published figures inherit whatever
identity model we land on here.

## 2. Design: one canonical gene identity, three published forms

The mistake to avoid is swapping one identifier for another — that just moves which peers
we fail to match. Instead:

- **Internally**, resolve every gene reference (ours *and* an inbound peer's) to a
  canonical identity, and compare on that.
- **On the wire**, publish Ensembl as `id` plus the alternatives as extension fields.

The spec explicitly blesses the second half. `search-api.md` §Non-Standard Fields gives
this as its worked example — near-verbatim what we want:

```json
"gene": {
    "id": "ENSG00000133392",
    "_ensemblGeneID": "ENSG00000133392",
    "_entrezGeneID": "4629",
    "_geneName": "MYH11"
}
```

So a v1.0 peer matching on symbols can still find us via `_geneName`, and a peer doing
proper resolution gets an unambiguous ID. We already use this `_`-prefixed convention for
`_derivedFrom` in `classification_ontology_slots()`, so it is an established pattern in
this codebase rather than a new liberty.

### What already exists — reuse it

The app already answers "which gene does this classification relate to", in two places,
and MME should build on them rather than beside them.

**`classification_gene_symbol_filter()`** (`classification/models/classification_utils.py:325`)
— the query side, used by the classification datatable (`classification_datatables.py:382`)
and candidate search. Given a symbol it builds a `Q` matching on three routes at once:
VEP annotation gene (`variantannotation__gene__in`), the `gene_symbol` evidence key across
**every alias**, and transcript (`transcriptversion__gene_version__gene__in` against
`allele_info__grch37/grch38__transcript_version__transcript`).

**`_gene_symbol_search_for()`**
(`classification/signals/classification_hooks_grouping_search_terms.py:16`) — the
denormalising side. Per `ClassificationGrouping` it extracts gene symbols from **four**
routes, each tagged in `extra`:

| tag | source |
|---|---|
| `from_normalized` | `ImportedAlleleInfo.resolved_builds[].gene_symbol` |
| `from_imported` | the imported c.HGVS symbol, falling back through `GeneSymbolAlias` |
| `from_transcript` | `GeneVersion.objects.filter(transcriptversion__in=…).gene_symbol` |
| `from_allele` | VEP `VariantAnnotation.gene.get_symbols()` |

and writes them upper-cased into `ClassificationGroupingSearchTerm`, indexed on
`("term_type", "term")`.

Two things follow:

- **Alias handling is solved** — `GeneSymbol.alias_meta` (`genes/models.py:235` →
  `GeneSymbolAliasesMeta`) gives `alias_symbol_strs` and carries a `different_genes` flag
  marking aliases that are *not* safe for automatic expansion. Use it rather than querying
  `GeneSymbolAlias` directly, which would miss that distinction.
- **The gap is real and narrower than it looked.** All of the above is keyed on
  **symbols**. Nothing in the app currently resolves a classification to an Ensembl gene
  id, so §1's problem stands — but the work is *adding an identifier dimension* to
  established extraction, not building extraction from scratch.

### Resolution sources for the Ensembl id (in `genes/`)

**`GeneVersion` is the workhorse** — it already maps identifiers to symbols in both
directions (`genes/models.py:510`, FKs to `Gene`, `GeneSymbol` and `GenomeBuild`). Filter
on the consortium and the Ensembl id falls out:

```python
GeneVersion.objects.filter(gene_symbol__symbol=symbol,
                           gene__annotation_consortium=AnnotationConsortium.ENSEMBL)
# -> .gene.identifier is the ENSG
```

and the reverse, `filter(gene__identifier="ENSG00000133392") -> .gene_symbol`, covers
inbound. Note `AnnotationConsortium.ENSEMBL` is stored as `"E"` (`genes/models_enums.py:4`
— `TextChoices`, `max_length=1`), with `"Ensembl"` as the display label.

Two properties worth having on purpose:

- It is populated by the ordinary gene annotation import that every deployment runs, so it
  needs no extra data load.
- Rows are per annotation version, so **superseded symbols still resolve**. A peer querying
  with a symbol we have since renamed hits an older `GeneVersion` row and lands on the
  right gene — which is precisely the inbound case we are trying to stop losing.

Supporting sources:

1. **Preferred entry point — the classification's resolved transcript.**
   `ResolvedVariantInfo` (`classification/models/classification_variant_info_models.py:93`)
   carries `transcript_version` and `gene_symbol` FKs; `TranscriptVersion.gene_version`
   reaches the `GeneVersion` directly. Unambiguous, because it starts from the curated
   transcript rather than a symbol string. Clinical labs curate against RefSeq (`NM_`), so
   this lands on the RefSeq `GeneVersion` — take its `gene_symbol` and re-query for the
   Ensembl row.
2. **Fallback — the gene symbol** from `SpecialEKeys.GENE_SYMBOL`, via `GeneSymbol.cast()`
   (cached, `genes/models.py:205`). Ambiguous when a symbol maps to several Ensembl genes;
   see §5.
3. **Entrez** — a `Gene` with `annotation_consortium == REFSEQ` has the Entrez gene id as
   its `identifier` (`Gene.prefixed_identifier` renders it `GeneID:####`). MME's
   `_entrezGeneID` example is bare (`"4629"`), so strip the prefix.
4. **Aliases** — via `GeneSymbol.alias_meta` (`genes/models.py:235`) for inbound symbols
   outside the `GeneVersion` history, honouring its `different_genes` flag.
5. **HGNC** — `HGNC.ensembl_gene_id` (`genes/models.py:88`), reachable via
   `GeneVersion.hgnc`. A useful cross-check and a second opinion on ambiguity, but not the
   primary route: it is a separate data load, so a deployment without it should still
   resolve normally.

## 3. Work

### 3a. `mme/genes.py` (new) — the resolver

One module owning gene identity for MME, so outbound, inbound and metrics cannot drift.

- `resolve_gene(classification) -> GeneIdentity | None` — start from the symbol the
  existing extraction already produces (`ResolvedVariantInfo.gene_symbol`, the same
  `from_normalized` route as `_gene_symbol_search_for`), then add the identifier dimension
  via `GeneVersion`.
- `gene_identity_from_id(gene_id: str) -> GeneIdentity | None` — for inbound: sniff the
  form (`ENSG…` / all-digits / `GeneID:…` / symbol), resolve to the same canonical
  identity via `GeneVersion`, then `GeneSymbol.alias_meta`.
- `GeneIdentity` dataclass: `ensembl_gene_id`, `entrez_gene_id`, `symbol`, plus
  `match_keys() -> set[str]` (every identifier form this gene is known by, normalised —
  including superseded symbols from older `GeneVersion` rows and
  `alias_meta.alias_symbol_strs`) and `as_mme_gene() -> dict` (the `id` + `_`-prefixed
  block above).

`match_keys()` is the crux: it makes comparison identifier-agnostic without needing both
sides to agree on a preferred form. It is the same expansion
`classification_gene_symbol_filter` performs for the datatable, widened from symbols to
identifiers — so honour `alias_meta`'s `different_genes` flag the way that code does, and
keep the two consistent.

### 3b. Outbound — `mme/serializers/patient_profile.py`

`classification_genomic_feature()` currently builds `{"id": str(gene_symbol)}`. Replace
with `resolve_gene(classification).as_mme_gene()`, keeping the bare-symbol form when
resolution yields nothing (§5).

### 3c. Inbound — `mme/matching.py`

`_genes_from_profile()` becomes: for each `genomicFeatures[].gene.id`, resolve via
`gene_identity_from_id()` and union its `match_keys()`; unresolvable ids fall back to the
current upper-cased-string behaviour so an unknown gene still self-compares. Both sides of
the Jaccard go through the same function, so `ENSG00000133392` vs `MYH11` intersects.

Also read the `_ensemblGeneID` / `_entrezGeneID` / `_geneName` extension fields when a
peer supplies them — free precision, no resolution needed.

Watch the N+1: `find_matches()` scans up to `MAX_CLASSIFICATIONS_SCANNED = 2000`
classifications per inbound query, resolving each one's profile in Python. The query's own
genes resolve once; ours want a bulk prefetch or a `timed_cache` keyed on symbol, in the
spirit of `GeneSymbol._cast`.

There may be a better answer than caching. `ClassificationGroupingSearchTerm` is an indexed
`(term_type, term)` lookup that already answers "which records relate to these gene
symbols" — `ClassificationGroupingSearchTerm.filter_q(GENE_SYMBOL, match_keys())` would
turn the gene half of the scan into one indexed query, and the term rows are maintained by
the signal rather than by us. Worth evaluating first, with one thing to confirm: those rows
key on `ClassificationGrouping`, while MME eligibility is per `ClassificationModification`,
so the two sets need to line up (a grouping spans an allele × lab, which may be coarser
than the record-level consent gate MME requires). If they do not line up cleanly, prefer
the cache and leave the scan record-level — the consent gate outranks the speed-up.

### 3d. Metrics — `mme/metrics.py`

`unique_genes` should key on the canonical identity rather than the upper-cased `id`, so
one gene counts once regardless of the form it was published in. This changes published
numbers; note it in the changelog since the figures are public.

### 3e. Settings

`MME_GENE_ID_PREFERENCE = "ensembl"` (`"symbol"` to retain current behaviour). Deployments
with no Ensembl `GeneVersion` rows set `"symbol"` and are unaffected. Default `"ensembl"`,
per the spec's direction of travel.

## 4. Decision to confirm before building

**Changing `id` from symbol to ENSG is peer-visible.** A v1.0 peer that matches on
symbols and ignores extension fields will stop matching us the day we switch. The spec
wants Ensembl; six of the eight deployed nodes are v1.0.

Recommendation: switch `id` to Ensembl and carry `_geneName` — the spec is unambiguous
about direction, and any peer sophisticated enough to be running v1.0 in 2026 handles
symbols in `_geneName` or does its own resolution. But it is a judgement call about other
people's implementations, and worth raising at the September consortium meeting rather
than deciding unilaterally. `MME_GENE_ID_PREFERENCE` exists so the answer is one setting.

## 5. Degradation rule

Never fail a submission or a match over gene resolution. If a gene cannot be resolved:
emit the symbol as `id` exactly as today, omit the extension fields, and log at debug.
A resolver returning `None` must be an ordinary path, not an exception — annotation
completeness varies per deployment and MME participation cannot depend on it.

## 6. Tests

- `resolve_gene` from a RefSeq-transcript classification reaches the right ENSG by
  re-querying `GeneVersion` on the symbol with `annotation_consortium = ENSEMBL`.
- `gene_identity_from_id` maps ENSG / bare Entrez / `GeneID:` / current symbol /
  superseded symbol / alias to one identity, and an alias flagged `different_genes`
  stays out of the automatic expansion.
- **The regression that motivates this**: an inbound query with `gene.id = ENSG…` scores
  > 0 against a classification we hold by symbol. Assert against `find_matches()`, since
  that is where the loss occurs today.
- Unresolvable gene → symbol emitted, no extension fields, submission still builds.
- `MME_GENE_ID_PREFERENCE = "symbol"` reproduces current output byte for byte.

Fakes live in `mme/tests/fakes.py`; gene/transcript fixtures in
`annotation/tests/test_data_fake_genes.py`.

## 7. Out of scope

- `genomicFeatures[].variant` coordinates — unchanged.
- Transcript-level identifiers on the wire; MME has no slot for them.
- v2.0 conformance generally. This removes the gene-id blocker; the rest of v2.0 is not
  specified yet.
