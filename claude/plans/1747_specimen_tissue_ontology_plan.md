# #1747 — specimen tissue as a UBERON ontology term

[#1747](https://github.com/SACGF/variantgrid/issues/1747), split out of
[#1706](https://github.com/SACGF/variantgrid/issues/1706) once the "there's currently no way to create
tissue" comment turned out to be bigger than that issue's scope.
[`1706_remaining_work_plan.md`](1706_remaining_work_plan.md) is the rest of #1706 — the specimen and
extraction grids — and the two are independent.

## The problem

`Tissue` (`patients/models.py:315`) is inert. `Specimen.tissue` is a nullable FK to it, and both
`SpecimenForm` (`patients/forms.py:137`) and `PatientSpecimenFormSet` (`patients/forms.py:122`) render
a tissue `<select>` that is empty on every deployment. There is no creation path outside
`patients/admin.py:13`, no seed data (`Tissue` appears only in `patients/migrations/0001_initial`), no
API field (`SpecimenSerializer`, `patients/serializers.py:184`), and the patient CSV importer has it
commented out at `patients/import_records.py:360`.

#1706 suggests seeding it from `HumanProteinAtlasTissueSample`. That does not hold up: those rows are
populated by `human_protein_atlas_import`
(`annotation/management/commands/human_protein_atlas_import.py:23`) off a downloaded HPA TSV, so they
only exist on deployments that ran that annotation import, and they are HPA's expression panel — full
of things nobody biopsies for a diagnostic workup, thin on the ones labs actually accession.

**The answer is not a better seed list, it is the ontology app.** A local `Tissue` table with a
hand-written list is a vocabulary we would own, version and reconcile ourselves, which is exactly what
`ontology` exists to avoid. DOID is queued to be added there in the standard manner, so the machinery
is being touched anyway and this rides along.

`Tissue` keeps no page either way — settled in TSO 500 Phase 3, and an `OntologyTerm` brings its own.

## 1. Which ontology — UBERON

`Specimen.tissue` is an **anatomy** question: which tissue the material is, with `tissue_status`
already carrying the role it plays and `Extraction.nucleic_acid_source` the arm taken off it. UBERON
is the anatomy ontology for exactly this slot:

- **GA4GH Phenopackets recommends UBERON for `Biosample.sampled_tissue`** — the same field, in the
  schema VariantGrid would emit if it ever exchanges biosamples.
- **CC BY 3.0**, so it can be referenced and shipped without a licence negotiation.
- **OBO/OWL, which `ontology_import` already parses with `pronto`** (`load_hpo`,
  `ontology/management/commands/ontology_import.py:294`) — the loader is the same shape as the ones
  already there.
- Covers what accessioning actually sends: blood `UBERON:0000178`, bone marrow `UBERON:0002371`,
  buccal mucosa `UBERON:0006956`, saliva `UBERON:0001836`, skin of body `UBERON:0002097`, plus every organ a
  solid-tumour block comes from.
- It is multi-species and ~21k classes, but that is the same order as MONDO and HPO, which the app
  already holds and serves through autocomplete rather than a dropdown. UBERON also ships subsets
  (`uberon_slim`, `organ_slim`, `major_organ`) if the candidate list ever needs narrowing; `pronto`
  reads the `subset:` tags, so that is a filter at import time rather than a different source.

Considered and not chosen:

- **SNOMED CT specimen hierarchy** (`123038009 |Specimen|`) is the LIMS-native vocabulary and the best
  fit for what a lab actually receives, since it carries the material and its topography together.
  It needs a SNOMED Affiliate licence plus registration to download — free for Australia through the
  NCTS, but not redistributable, and variantgrid.com deployments outside Australia would each need
  their own. It cannot be a default. Worth keeping as a mapping target for a deployment that already
  has it, which the design below allows.
- **BRENDA Tissue Ontology (BTO)** — tissues and cell lines, freely available, but less adopted for
  clinical specimens and with no phenopacket alignment.
- **NCIt** — its value is the tumour axes (progression, grade, histology), which is what phenopackets
  assigns it; for plain anatomy it is a larger and less precise UBERON.
- **EFO** — phenopackets' recommendation for `sample_type`, i.e. the *material* axis below, not this
  one.

**Preservation and processing are a different axis.** FFPE, fresh frozen, an EDTA tube — none of that
is anatomy, and none of it belongs on a UBERON FK. Phenopackets splits it out as
`sample_processing` / `sample_type`. Out of scope here; §3 keeps the raw text, which is where an
`FFPE block` value with no site in it will sit until someone needs that axis modelled.

## 2. The import, in the standard manner

Alongside the DOID work, and following MONDO/HPO:

- `OntologyService.UBERON = "UBERON", "UBERON"` in `ontology/models/models_ontology.py:43`, with a
  7-digit `EXPECTED_LENGTHS` entry, an `IMPORTANCE` entry and a `URLS` entry pointing at OLS.
- Add it to `LOCAL_ONTOLOGY_PREFIXES` (`ontology/models/models_ontology.py:88`) and leave it **out of**
  `CONDITION_ONTOLOGIES` (line 95). That one line is what keeps anatomy out of condition text matching,
  the classification autocomplete and the condition search receivers — all of which already gate on
  that set.
- `OntologyImportSource.UBERON`, and a `load_uberon()` beside `load_hpo` taking `--uberon`.
- Leave it out of `OntologyVersion.ONTOLOGY_IMPORTS` (`ontology/models/models_ontology.py:799`).
  Those imports are pinned because `AnnotationVersion` replays historical `OntologyTermRelation`
  traversals; a specimen stores a term id and traverses nothing, so pinning it would add a required
  import to every deployment's ontology version for no gain.
- A search receiver via `_ontology_search_name` (`ontology/models/ontology_search.py:85`) and an
  autocomplete subclass in `ontology/views_autocomplete.py:60` onwards.

## 3. Take what the LIMS gives

A hard FK to an `OntologyTerm` would mean rejecting or dropping whatever a LIMS sends that does not
resolve, and LIMS send `FFPE`, `EDTA blood`, `BM aspirate`, `Buccal` and local codes. A specimen with
an unrecognised tissue string is still a perfectly good specimen, so **the text is the record and the
term is the interpretation**:

- `Specimen.tissue_text` — what arrived, unmodified.
- `Specimen.tissue_term` — FK to `OntologyTerm`, `limit_choices_to` UBERON, null until resolved.
- Match status, so an unresolved value is visibly unresolved rather than silently absent.

Two precedents to follow rather than invent against:

- **`ExtractionMatchMixin`** (`patients/models.py:478`) — already in this app, already the shape: a
  claim that may not be resolvable yet, carrying status, error and the date the claim was parked, with
  a settled link left alone so it never flaps back.
- **`ConditionText` / `ConditionTextMatch`** (`classification/models/condition_text_matching.py:61`) —
  normalized text keyed per lab, resolved to standard terms, human-confirmable. The same LIMS sends
  the same twenty strings forever, so a per-lab mapping resolves every future specimen once someone
  has mapped it once. `normalize_condition_text` and `ontology/ontology_matching.py` are the existing
  machinery.

The specimen page shows the text with its matched term, and lets a user with write permission pick the
term — an `OntologyTerm` autocomplete forwarding `ontology_service=UBERON`, which replaces the empty
`<select>` that started this. A deployment holding SNOMED can map its own codes into the same field
without VariantGrid shipping SNOMED.

Nothing needs migrating: `Tissue` has no rows to carry across on any deployment we know of — no seed,
no importer path, no API field ever wrote one. Confirm that against the real deployments, then drop
the model with the FK.

## 4. Where tissue gets set

- **Specimen page** — the autocomplete above. This is the path that makes the field usable at all.
- **Patient CSV** — no tissue column exists. `PatientRecord` has `specimen_tissue_status` but no
  `specimen_tissue` (`patients/models.py:835`), and `PatientColumns.COLUMN_DETAILS`
  (`patients/models.py:661`) has no entry, which is why `patients/import_records.py:360` is commented
  out. Adding one is a column constant, a `PatientRecord` field, a migration and a match attempt —
  cheap once §3 exists, because an unmatched name parks as text instead of failing the row.
- **API** — `SpecimenSerializer` (`patients/serializers.py:184`) omits tissue. Same story: accept the
  text, return the matched term. #1707's client work has not been written, so this is cheapest to
  settle before that ships.

## Suggested order

1. **§2** — the UBERON import, riding with DOID while `ontology_import` is open anyway. Independent of
   the rest and useful on its own.
2. **§3 and §4's specimen page** — the model change and the autocomplete. This is what turns a visibly
   broken dropdown on a page that already shipped into a working one.
3. **§4's CSV column and API field** — whenever a consumer needs them. Both are cheap once §3 exists,
   because an unmatched name parks as text rather than failing a row.

While this is outstanding, #1706 hides the field: `'tissue'` joins the excludes in `SpecimenForm` and
`PatientSpecimenFormSet`, so the page stops offering a select that can never be filled.

## Done when

- `ontology_import --uberon` loads UBERON terms, they are searchable and autocompletable, and they
  reach neither the condition text matching UI nor the classification condition autocomplete.
- The specimen page offers a UBERON autocomplete, saves the choice, and shows the LIMS text beside the
  matched term — including where nothing matched, which reads as unresolved rather than empty.
- A LIMS string that has been mapped once resolves by itself on the next specimen that carries it, and
  one that has never been seen parks the specimen's tissue as text without failing its import.
