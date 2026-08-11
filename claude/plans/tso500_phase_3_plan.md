# TSO 500 Phase 3 — tissue status, and pages for Specimen / Extraction

Design for the Phase 3 row of [`tso500_overall_plan.md`](tso500_overall_plan.md): two Phase 1 leftovers,
[SACGF/variantgrid_private#2447](https://github.com/SACGF/variantgrid_private/issues/2447) (replace
`Specimen.mutation_type` with a tissue status) and [#1706](https://github.com/SACGF/variantgrid/issues/1706)
(specimen and extraction pages, search, preview).

Both are markedly cheaper before Phase 4 puts `Specimen` behind a public API, and both are the last
model-shaped work in `patients` before #1707 serializes it.

---

## Scope

In:

- `Specimen.mutation_type` → `Specimen.tissue_status`, everywhere it reaches: autopopulate, the patient
  CSV, `PatientRecord`, admin, forms.
- Detail pages for `Specimen` and `Extraction`, with permissions, `PreviewModelMixin`, a search handler,
  and links each way along `Patient → Specimen → Extraction → Sample`.

Out, riding with Phase 6 where the grid column pass already is:

- Top-level specimen / extraction grids under the patients menu.
- An extraction link on the sample grid.
- private#2837's specimen/patient IDs in the variant page's sample table.

The split is on "does Phase 4 need it": Phase 4's done-when is *the measures show on the specimen page*,
so the page is a prerequisite; the grids are not.

---

## Part A — `Specimen.tissue_status` (private#2447)

### What is there now

`Specimen.mutation_type` (`patients/models.py:321`) is a two-value Germline/Somatic `CharField` with
`default=Mutation.GERMLINE`. Its one substantive consumer is classification autopopulate
(`classification/autopopulate_evidence_keys/evidence_from_sample_and_patient.py:91`):

```python
data[SpecialEKeys.ALLELE_ORIGIN] = specimen.get_mutation_type_display()
```

`allele_origin` feeds `AlleleOriginBucket.bucket_for_allele_origin`
(`classification/enums/classification_enums.py:24`), which partitions `AlleleOriginGrouping` and the
cross-lab overlaps. So every TSO 500 tumour block, having never had the field set, stamps its
classifications "Germline" and files them alongside other labs' real germline records — corrupting
exactly the somatic classifications Phase 8 exists to produce.

The rest of the reach is mechanical: `patients/admin.py:20`, the specimen formset
(`patients/forms.py:121`, which excludes only `external_pk` so it picks the field up automatically),
`PatientColumns.SPECIMEN_MUTATION_TYPE` and its queryset path (`patients/models.py:533`),
`PatientRecord.specimen_mutation_type` (`patients/models.py:662`), `patients/import_records.py:241,
363, 373, 419`, and `patients/templates/patients/view_patient_record.html:53`.

### The field

In `patients/models_enums.py`, replacing `Mutation` (which has no other consumer):

```python
class TissueStatus(models.TextChoices):
    REFERENCE = 'R', 'Reference (unaffected)'   # hair, buccal, blood in a solid-tumour workup
    AFFECTED = 'A', 'Affected / lesional'       # tumour block, affected skin in mosaic disease
    UNKNOWN = 'U', 'Unknown'
```

On `Specimen`, `default=TissueStatus.UNKNOWN`, not null. The current `null=True, blank=True` plus a
Germline default is where the silent damage comes from: an unset field is indistinguishable from a
deliberate one. `UNKNOWN` is an explicit value that says so.

This lives on `Specimen` rather than on the `Tissue` lookup (`patients/models.py:303`) because the same
tissue plays different roles — blood is the reference in a solid-tumour workup and *is* the tumour in
leukaemia; skin is reference normally and is the affected tissue in Proteus. That is the mosaic case the
issue was opened about, and it is why a per-specimen field earns its keep instead of tagging tissue types.

### Three levels, and one correction to the issue

| Level | Question | Where | Values |
|---|---|---|---|
| Specimen | What is this material, and what role does it play in the test? | `Specimen.tissue_status` (new) | Reference / Affected / Unknown |
| Call set | What's in this VCF? | `Sample.variants_type` (exists) | Unknown / Germline / Mixed / Somatic only |
| Variant | What is *this* variant's origin? | `allele_origin` → `AlleleOriginBucket` (exists) | germline / somatic / other |

The issue says `VCF.variants_type`. It is actually on **`Sample`** (`snpdb/models/models_vcf.py:343`),
alongside `Sample.extraction` and `Sample.is_somatic`. That matters for where the derivation goes.

### Autopopulate becomes a derivation

`get_evidence_fields_for_extraction` (`:87`) takes only an `Extraction`, so it cannot see
`variants_type`. Move the origin decision up into `get_evidence_fields_for_sample_and_patient` (`:34`),
which already holds both the sample and `sample.extraction` (`:45`), and leave the extraction helper
doing the things that are genuinely per-extraction (specimen ID, sample type, age, nucleic acid source).

The rule — assert an origin only where both levels agree:

| `tissue_status` | `variants_type` | `allele_origin` |
|---|---|---|
| Reference | Germline | `germline` |
| any | Somatic only (tumour minus normal) | `somatic` |
| anything else | anything else | unset |

Unset falls back to `settings.ALLELE_ORIGIN_NOT_PROVIDED_BUCKET` so the curator sets it per variant.
That last row is the heart of the issue: for a mixed tumour sample the origin is genuinely per-variant
and cannot be known at accession.

Emit the lowercase option strings (`"germline"`, `"somatic"`) rather than a `get_..._display()` label —
`bucket_for_allele_origin` looks the value up through `EvidenceKeyMap.cached_key(...).matched_options()`,
and the evidence key's own option keys are what that matches.

### Migration

Three migrations' worth of work, in one file per Django convention (next is `patients/0016_`):

1. Add `tissue_status` with `default=TissueStatus.UNKNOWN`.
2. `RunPython` mapping `mutation_type` `S` → `AFFECTED`, `G` and null → `UNKNOWN`.
3. Remove `mutation_type`.

`S` → `AFFECTED` is clean. `G` is ambiguous because it was the default — it may mean "a curator said
germline" or "nobody touched it", and there is no way to tell the two apart from the column alone. On
this dev box every `PatientRecord.specimen_mutation_type` is null, but that says nothing about SA
Pathology or variantgrid.com.

So rather than guess, map `G` → `UNKNOWN` and register a `ManualOperation` whose `test=` callable only
fires on deployments that actually have evidence of a deliberate germline declaration — a
`PatientRecord` row with `specimen_mutation_type='G'`, meaning the CSV column was populated rather than
defaulted. Following `snpdb/migrations/0188_one_off_migrate_common_filter_gnomad_versions.py`. On a
deployment that never populated the column the task is never registered and nothing is lost; on one that
did, a human reviews those specimens.

`PatientRecord.specimen_mutation_type` → `specimen_tissue_status` with the same value mapping, in the
same migration. `PatientRecord` is an audit of what an uploaded CSV row said, so it follows the
vocabulary the CSV now speaks.

### The patient CSV

`PatientColumns.SPECIMEN_MUTATION_TYPE = 'Specimen Mutation type (Germline/Somatic)'` carries its
vocabulary in the header text, so it becomes:

```python
SPECIMEN_TISSUE_STATUS = 'Specimen Tissue status (Reference/Affected/Unknown)'
```

`parse_choice` (`patients/import_records.py:119`) already accepts either the key or the display value in
any case, so `R`, `Reference (unaffected)`, `AFFECTED` all parse without further work.

Import validates that *every* column in `PatientColumns.COLUMNS` is present
(`patients/import_records.py:455`) and raises otherwise, so a saved spreadsheet with the old header
fails the import naming the column it is missing. That is the intended behaviour — Germline → Unknown
is a loss of meaning rather than a translation, so silently reinterpreting an old column would be
worse than making people re-download the template from `example_upload_csv_empty`
(`patients/views.py:221`). The existing error message already names the missing column, so the failure
says what to do.

### Everything else in Part A

- `patients/admin.py:20` — `mutation_type` → `tissue_status` in `list_display`.
- `patients/forms.py:121` — nothing, the formset excludes only `external_pk`.
- `patients/templates/patients/view_patient_record.html:53` — label and accessor.
- `patients/models.py:533` — `SAMPLE_QUERYSET_PATH` entry, which drives the "download all samples as CSV"
  round trip (`patients/views.py:230`).
- `patients/tests/test_specimen_extraction.py:50` — uses `mutation_type=Mutation.SOMATIC`.

---

## Part B — pages, search, preview (#1706)

### Permissions first

`Specimen` and `Extraction` are plain `ExternallyManagedModel` (`patients/models.py:311`, `:363`): no
guardian permissions of their own, which is right — a specimen's confidentiality is its patient's. The
codebase already has the pattern for delegating, in `Trio` (`snpdb/models/models_cohort.py:780`, `:795`,
`:800`):

```python
class Specimen(GuardianPermissionsMixin, ExternallyManagedModel, PreviewModelMixin):
    @classmethod
    def get_permission_class(cls):
        return Patient

    def get_permission_object(self):
        return self.patient

    @classmethod
    def _filter_from_permission_object_qs(cls, queryset):
        return cls.objects.filter(patient__in=queryset)
```

and `Extraction` the same through `specimen__patient__in`. That buys `can_view`, `can_write`,
`check_can_view`, `filter_for_user` and `get_for_user` for free, and `get_for_user` loads the full
instance where permissions delegate (`library/django_utils/guardian_permissions_mixin.py:99`), so the
delegation resolves correctly.

One wrinkle: `ExternallyManagedModel.can_write` already exists on both models and returns whether the
external manager permits modification. `Patient` resolves the clash by ANDing the two
(`patients/models.py:180`); do the same on both.

This also lets `patients/views_autocomplete.py:26,40` drop their hand-rolled
`filter(patient__in=Patient.filter_for_user(user))` for `Specimen.filter_for_user(user)`.

### Pages

Two views and two templates, plus URLs:

```python
path('view_specimen/<int:specimen_id>', views.view_specimen, name='view_specimen'),
path('view_extraction/<int:extraction_id>', views.view_extraction, name='view_extraction'),
```

Each view is `Specimen.get_for_user(request.user, pk)`, a `ModelForm` with `set_form_read_only` where
`can_write` is false, saved on POST — mirroring `view_patient` (`patients/views.py:33`). The fields are
the same ones the patient tabs' formsets already edit, so the forms are a `modelform_factory` over the
same field lists; the tabs stay as the bulk editor and the page is where one specimen lives.

`external_pk` is excluded from both forms and shown read-only beside the editable fields, the way the
specimen formset already excludes it (`patients/forms.py:121`). It is the tracking system's identity for
the record and Phase 4's resolver keys on it, so it is displayed rather than typed — and
`ExternallyManagedModel.can_write` separately locks the whole form where the external manager forbids
modification, which `ExternalModelManager.explaination` (`patients/models.py:48`) already has text for.

The specimen page shows: patient (linked), tissue, collection/received dates, age at collection,
`tissue_status`, external PK where set, then a table of its extractions (linked) with their samples
(linked). The extraction page shows: specimen (linked), patient (linked), reference, nucleic acid source,
extraction date, and its samples — plus its `SequencingSample` rows, which Phase 1 already pointed here
and Phase 4 will start populating.

`get_absolute_url` on both, which is what `PreviewData.for_object` and every future link need.

Add both URL names to the `patients` app's block in the settings URL register, so the shariant
deployment that turns `patients` off (`variantgrid/settings/env/shariantcommon.py:193`) turns these off
with it.

### Preview

`PreviewModelMixin` on both, following `Patient` (`patients/models.py:157`):

- `preview_if_url_visible()` → `'patients'`, matching `Patient`, `Sample` and `ExternalPK`, so the
  hover-cards vanish wherever the patients app is off.
- `preview_icon()` — a vial for `Specimen`, a flask/droplet for `Extraction`. Distinct from `Patient`'s
  `fa-solid fa-user-injured` and `Sample`'s `fa-solid fa-microscope`.
- `preview` — identifier is `reference_id or external_pk or f"({pk})"`, title is the patient, and
  `summary_extra` carries tissue status and collection date for a specimen, nucleic acid source and
  extraction date for an extraction.

### Search

A new `patients/signals/specimen_search.py`, registered the way `patient_search.py` is. `HAS_3_ANY`
rather than `HAS_3_ALPHA_MIN`, because a TSO 500 reference like `2600000001` is all digits.

```python
@search_receiver(search_type=Specimen, pattern=HAS_3_ANY,
                 example=SearchExample(note="Specimen reference", examples=["2600000001"]))
def specimen_search(search_input: SearchInputInstance):
    yield Specimen.filter_for_user(search_input.user).filter(search_input.q_words('reference_id'))
```

and an extraction receiver over `reference_id` — matching what `ExtractionAutocompleteView` already
searches (`patients/views_autocomplete.py:43`), which also searches `specimen__reference_id` so a
container suffix and its parent both find it.

Both go through `filter_for_user`, which Part B's first section just gave them. `ExternalPK` search
(`patients/signals/external_pk_search.py`) walks a hardcoded `RELATED_OBJECT_FIELDS` list of
`["case", "pathologytestorder", "patient"]` — add `specimen` and `extraction` so an external identifier
finds them too.

### Links each way

- Patient page: the specimens and extractions tabs are formsets, so add a link out of each row to its
  own page. (`patients/templates/patients/view_patient_specimens.html`,
  `view_patient_extractions.html`.)
- Specimen page → patient, extractions, samples. Extraction page → specimen, patient, samples.
- Sample page (`snpdb/templates/snpdb/data/view_sample.html:73`) already has an extraction autocomplete;
  give it a "View Extraction" cross-link beside the existing "View Patient" one (`:49`), the same
  `setCrossLink` pattern.

---

## Order of work

1. Part A's model change, migration and `ManualOperation` — it is the one ordering constraint in the
   phase, since #1707 must serialize the field it is keeping.
2. Part A's consumers: autopopulate, CSV, `PatientRecord`, admin, template.
3. Part B's permission delegation.
4. Part B's URLs, views, templates.
5. Part B's preview and search.

3–5 are views and templates on top of a settled model, so they parallelise with anything.

## Tests

- Extend `patients/tests/test_specimen_extraction.py`: the autopopulate truth table — reference+germline
  → `germline`, anything+somatic-only → `somatic`, affected+mixed → unset — asserted against
  `get_evidence_fields_for_sample_and_patient` rather than the helper, since that is where the
  derivation now lives.
- A migration test is not worth it, but a unit test that a `Specimen` created with no arguments comes
  out `UNKNOWN` is, because that default is the whole point of the issue.
- CSV round trip: the existing `process_record` tests in `test_specimen_extraction.py` cover the new
  column, and `_make_row` builds from `PatientColumns.COLUMNS` so it follows the rename automatically.
- `patients/tests/test_urls.py` already builds a specimen and an extraction in `setUpTestData`, so
  `view_specimen` and `view_extraction` go straight into `PRIVATE_OBJECT_URL_NAMES_AND_KWARGS`, which
  gives the 200-for-owner and 403-for-non-owner pair for free — and that 403 is the real test of the
  permission delegation.

## Done when

A specimen has its own page, reachable from search and from the patient, and `tissue_status` has
replaced `mutation_type` everywhere including the CSV — with a tumour block's classifications no longer
stamped "Germline" by default.

---

## Settled

- **The old CSV header is a hard break.** No alias — a spreadsheet carrying
  `'Specimen Mutation type (Germline/Somatic)'` fails the import naming the missing column, and people
  re-download the template. Germline → Unknown is a loss of meaning, so importing it silently would be
  worse than failing.
- **Whether `G` means anything is unknown**, so the migration keeps the `ManualOperation`. It is
  registered only where a `PatientRecord` actually carries `G` — evidence the CSV column was populated
  rather than defaulted — so on a deployment that never used it, no task appears and nothing is asked of
  anyone.
- **The pages are editable**, except `external_pk`, which is displayed read-only.
- **`tissue_status` is per-specimen only.** Both arms of one block share the material, so it does not go
  on `Extraction` the way `nucleic_acid_source` did.
- **`Tissue` gets no page.** It is a lookup table; admin covers it.
