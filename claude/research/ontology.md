# ontology — research notes

Verified against 7c4408c62 on 2026-09-06

The ontology app is the local mirror of the condition and phenotype vocabularies VariantGrid classifies against:
MONDO (the preferred disease ontology, because it is hierarchical and cross-references the others), OMIM, HPO, and
HGNC gene symbols stored as terms purely so gene-disease edges have a node to point at. Terms are a single unversioned
table; the *relations* between them - is-a, synonym and cross-reference edges from the MONDO and HPO files, gene-disease
assertions from GenCC and PanelApp Australia, phenotype-gene links from HPO's annotation file - are partitioned per
import so an old AnnotationVersion can still read the graph it was built against. Downstream, gene annotation
walks the graph to list a gene's conditions, patient phenotype text is matched to terms, and the classification app's
condition text matching resolves a lab's free-text condition to MONDO/OMIM ids. Fields, URLs, commands and settings are in
the generated maps ([models](../maps/models.md#ontology), [urls](../maps/urls.md#ontology),
[commands](../maps/commands.md), [settings](../maps/settings.md));
`ontology/__ontology_readme.md` is the short introduction.

## Flows

### One builder, many sources

Every source goes through `ontology/ontology_builder.py:OntologyBuilder`, which is the whole import contract: one
builder per (`import_source`, `context`, `processor_version`); `OntologyBuilder.ensure_hash_changed` raises
`OntologyBuilderDataUpToDateException` when the file's md5 matches the last completed import of that context;
`OntologyBuilder.cache_everything` pulls every existing term into memory; `add_term` / `add_ontology_relation` mutate
cached objects that `OntologyBuilder.complete` writes back with one bulk create and one bulk update per model. A previous import made with a different `processor_version` is ignored, so bumping the number
in a loader is how a parsing change forces a reload of an unchanged file. `add_term` distinguishes the
`primary_source` for a term (MONDO for MONDO ids, mimTitles.txt for OMIM, hp.owl for HPO) from a mere reference:
a non-primary caller never overwrites a term another file owns, and an untrusted reference that creates a term leaves
it with status `STUB`. Each builder creates its own `ontology/models/models_ontology.py:OntologyImport`; when that
import is one of the files an `OntologyVersion` tracks, `OntologyImport.save` adds a list partition of
`OntologyTermRelation` for it, so relations are always inserted into a fresh partition and old ones are never rewritten.

`ontology/management/commands/ontology_import.py:Command` runs the file loaders in one go. `load_mondo` reads the
MONDO JSON graph: nodes become terms with synonyms as aliases, `basicPropertyValues` and `xrefs` become one relation
per (term, other term) in priority order (exact, related, xref - one edge per pair, since several get in the way of
traversal), and axioms and edges become `is_a` edges between MONDO terms and gene edges to HGNC nodes named by the RO
relations in `GENE_RELATIONS`. `load_hpo` reads hp.owl through pronto for terms and
`is_a` edges; `load_phenotype_to_genes` reads HPO's phenotype_to_genes.txt for HPO->OMIM associations and
OMIM->gene `mim2gene` edges, deleting the relations of the older `OMIM_ALL_FREQUENCIES` file it replaced. `load_omim`
reads mimTitles.txt as the primary source for OMIM names and aliases and marks gene and removed entries
`NON_CONDITION`; deployments without an OMIM licence use `load_biomart` instead, which fills the same terms from
Ensembl BioMart's morbid descriptions. `sync_hgnc` mirrors the genes app's HGNC table into HGNC terms (previous
symbols as aliases) and is unversioned. The command ends by calling `OntologyVersion.latest`, which is what creates
the new version (below).

GenCC and PanelApp Australia are the two gene-disease sources and arrive differently. `ontology/gencc.py:load_gencc`
groups the GenCC submissions CSV by (gene, disease) and writes one `RELATED` edge MONDO->HGNC whose `extra` holds
every submitter's classification, MOI and PubMed ids plus `strongest_classification`; it runs from a file via the
command or from the web via `ontology/tasks/cached_web_resource_tasks.py:ClinGenCCWebResourceTask`. PanelApp is
live rather than a file: `ontology/panel_app_ontology.py:_update_gene_relations` fetches one gene's panels, deletes
that gene's existing `panelappau` edges and rewrites them (one per condition id found in the panel's phenotypes, with
the strongest panel confidence), under an *unversioned* builder whose context is the gene symbol, throttled by
`PANEL_APP_CACHE_DAYS` and the response hash. The per-symbol hook `ontology/panel_app_ontology.py:update_gene_relations`
runs on read paths only when `GENE_RELATION_PANEL_APP_LIVE_UPDATE` is set; batch work uses
`ontology/panel_app_ontology.py:bulk_update_gene_relations`, one paginated crawl of `/api/v1/genes/` that identifies
genes by the HGNC id PanelApp reports rather than its symbol (#1667), skipped when `panel_app_bulk_data_age` says the
crawl is fresh.

### Versions

`ontology/models/models_ontology.py:OntologyVersion` is one FK per tracked file (GenCC, MONDO, hp.owl,
phenotype_to_genes.txt, and optionally OMIM), unique together, and it is created lazily: `OntologyVersion.latest`
picks the newest completed import per slot, `get_or_create`s the version, and on creation calls
`AnnotationVersion.new_sub_version` so every build's AnnotationVersion now points at it. That is the link that makes
an ontology import a deploy event: `annotation/models/models.py:AnnotationVersion.validate` refuses a version whose
GeneAnnotationVersion was built for a different OntologyVersion, so `manage.py gene_annotation` must follow every
import (the command prints the `--ontology-version` to pass). Readers pick their graph through
`OntologyVersion.get_ontology_term_relations` (relations whose `from_import` is one of the version's five) or, for
live pages, `OntologyVersion.get_latest_and_live_ontology_qs`, which is the latest version's relations plus every
`panelappau` edge - PanelApp edges are outside versioning on purpose, since they are refreshed per gene in place.
Only relations are versioned; a term's name, aliases and status are whatever the last import wrote.

### Traversal

`ontology/models/models_ontology.py:OntologySnake` is an immutable path (source term, leaf term, list of relations
walked), and `OntologySnake.snake_from` is a breadth-first search from a term toward a target ontology within
`max_depth` steps, over `ontology/ontology_traversal.py:bfs_to_ontology`. Edges are filtered by an
`OntologyRelationshipQualityFilter`, which applies a minimum GenCC strength and a minimum PanelApp confidence only to
edges from those two sources and passes everything else; the presets run from `ONTOLOGY_RELATIONSHIP_NO_QUALITY_FILTER`
to `ONTOLOGY_RELATIONSHIP_STANDARD_QUALITY_FILTER` (GenCC Strong, PanelApp Green), which is the default everywhere a
gene-disease claim is shown. `OntologySnake.terms_for_gene_symbol` (gene to conditions) and
`OntologyVersion.gene_symbols_for_terms` (conditions to genes, used by the analysis phenotype node) are the two
directions; `OntologySnake.check_if_ancestor` and `OntologySnake.all_descendants_of` walk `is_a` only, and
`OntologySnake.has_gene_relationship` also tries the MONDO/OMIM `exact` counterpart of the term, because GenCC asserts
against MONDO while labs often pick OMIM. `ontology/ontology_traversal.py:get_ontology_traverser` gives batch callers
the same API over a `MemoryOntologyTraverser`, which loads one version's whole edge list into adjacency dicts once;
`gene_annotation` uses it across every HGNC term, where a query per step was the bottleneck. The HGNC bridge is
`ontology/models/models_ontology.py:OntologyTerm.get_gene_symbol`: name match, then alias match, then
`genes/gene_matching.py:HGNCMatcher` - and if HGNC knows the gene but no term exists it creates one on the spot under
an ad-hoc import. Its docstring is the rule: pass GeneSymbols around and only cross to the HGNC term at the edge of a
traversal.

### Condition text matching

A classification's condition arrives as free text. `classification/models/condition_text_matching.py:ConditionText`
is one row per (normalised text, lab), and under it a `ConditionTextMatch` hierarchy - root, gene symbol, mode of
inheritance, classification - where a level with terms set resolves every classification below it that has none of
its own (`condition_text_match_for_classification_id` walks up to the first valid level).
`ConditionTextMatch.sync_condition_text_classification` builds or repairs that tree for one published modification
(gene symbol from the evidence key, else an unambiguous alias, else the c.HGVS), and runs from the publish signal
(`classification/models/condition_text_matching.py:published`) and the withdrawn-flag signal;
`classification/management/commands/sync_condition_text_matches.py:Command` rebuilds everything and, with
`--obsolete`, reports matches whose terms an ontology import has since deprecated.

Suggestions come from `classification/models/condition_text_matching.py:top_level_suggestion`: first
`embedded_ids_check` (an `OMIM:123456` or `MONDO:` id written in the text, with a validation that the surrounding
words resemble the term, and a trailing "uncertain" / "co-occurring" setting the multi-condition operation), then
`search_suggestion` over the local terms using `ontology/ontology_matching.py:SearchText` (tokenised, de-pluralised,
roman numerals normalised, split at "type") and the Monarch MONDO search. `ConditionTextMatch.attempt_automatch`
assigns a suggestion without a human only when `ConditionMatchingSuggestion.is_auto_assignable` holds: a single term,
not via an alias, no warnings, and either an embedded id at the root or, at gene level, a leaf term with a
standard-quality relationship to that gene. `apply_condition_resolution` then copies the resolved terms onto each
affected classification's `condition_resolution`, adding a flag comment recording the change, so exports and
discordance read a cached column rather than the tree. `ontology/ontology_matching.py:OntologyMatching` is the same
machinery for the interactive picker, ranking terms by gene relationship, search score and direct reference.

## Why it is shaped this way

Terms and relations are split because they age differently: a MONDO id is stable and a rename should be visible to
every old record, while whether gene X causes disease Y is an assertion that changes with each GenCC or PanelApp
refresh and must be reproducible for an old AnnotationVersion. Partitioning relations by import (#206, 2022) made
"which edges did this version see" a partition pruning rather than a filter over one huge table, and made dropping a
superseded import cheap. Storing gene symbols as HGNC terms was the price of a single relation table - every edge has
two `OntologyTerm` ends - and `get_gene_symbol` exists to keep that leakage contained. `OntologyVersion` is created
lazily from the latest imports rather than by a command because the files are loaded one at a time over several
deploys, and the version should reflect whatever combination is complete. PanelApp is unversioned and live because it
is refreshed per gene on demand; GenCC carries PanelApp Australia's submissions too, so a deployment with live updates
off still sees them, and `load_gencc` drops the PanelApp rows when live updates are on to avoid double counting.
Condition matching lives in classification, not here, because it is about a lab's text and permissions (each
`ConditionText` is Guardian-scoped to its lab); the ontology app only supplies terms and search.

## History

The app was split out of annotation in January 2021 (`ontology/migrations/0001_initial.py`). Relation partitioning
and `OntologyVersion` came together in June 2022 (#206, `ontology/migrations/0013_one_off_ontology_term_relationship_partitions_and_copy.py`,
`ontology/migrations/0015_ontologyimport_version_ontologyversion.py`), followed by `status` on terms so every OMIM
entry could be stored but only conditions offered (`ontology/migrations/0020_ontologyterm_status.py`, 2022), a
de-duplication of PanelApp edges (`ontology/migrations/0022_one_off_remove_panelappau_duplicates.py`) and gene-disease
relations per version (#648). External ontologies were added as non-local services - DOID, Orphanet, MedGen, MeSH
(2024) - and internal term links became a setting (#700). Strength filtering was tightened in 2024
(variantgrid_private#3646), the in-memory traverser landed in May 2026, unsupported Monarch prefixes (MPATH) stopped
raising (#1603), and PanelApp gene identity moved to HGNC ids (#1667). On the classification side, condition text accepted multiple terms with auto-assignment in
2023 (#897) and the resolved-condition history began to be tracked (#881).

## Traps

`OntologyTerm.is_stub` means the object is unsaved (`_state.adding`), not that `status == STUB`:
`ontology/models/models_ontology.py:OntologyTerm.get_or_stub` returns a transient term for any well-formed id it does
not hold, and callers must check `is_stub` before joining on it. `OntologyTerm.is_obsolete` is "not valid as a
condition": true for deprecated terms, stubs and OMIM gene or removed entries (`NON_CONDITION`), while HGNC terms
are stored as `CONDITION` and pass. `OntologyVersion.latest` writes (it may create
a version and bump every AnnotationVersion) and raises `OntologyVersion.DoesNotExist` until GenCC, MONDO and both
HPO files have all been imported once; the unique constraint includes the nullable `omim_import`, so two versions
differing only in OMIM-null are both allowed. Aliases on a term are replaced wholesale by its primary source, but
`get_gene_symbol` appends to an HGNC term's aliases in place - run `ontology_import --hgnc_sync` after an HGNC
refresh or the two drift. The `OntologyTermRelation` default manager `select_related`s both ends and the import.
A relation's direction is conventional (MONDO->OMIM, MONDO->HGNC,
OMIM->HGNC) but every reader uses `OntologyTermRelation.other_end` and checks both directions, so never assume
`source_term` is the disease. `sync_condition_text_matches --obsolete` is the only thing that notices when an
ontology import deprecates a term a lab has already matched.
