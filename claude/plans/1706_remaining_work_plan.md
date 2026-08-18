# #1706 — remaining work: specimen and extraction grids

[#1706](https://github.com/SACGF/variantgrid/issues/1706) asked for four things: pages for `Specimen`
and `Extraction`, search, links from the model that contains them, and a decision on top-level grids.
The first three landed with TSO 500 Phase 3 (PR #1715) and the patient page work in `9c4354e5a`. What
is left is the grid half, already carried as part of Phase 6 in
[`tso500_overall_plan.md`](tso500_overall_plan.md).

This doc is the spec for it. Phase 6's other half (#1558 — CN/SM columns, fusion columns, fusion
search by a single partner) stays in the overall plan; nothing here depends on it, and it does not
depend on this. The issue's tissue comment became #1747 — see below.

---

## Where it stands

Done, and not to be redone:

| Piece | Where |
|---|---|
| Detail pages, editable, permission-gated | `patients/views.py:127` `view_specimen`, `patients/views.py:152` `view_extraction`, templates in `patients/templates/patients/` |
| Forms | `SpecimenForm` / `ExtractionForm`, `patients/forms.py:137` |
| URLs both grids link to | `view_specimen`, `view_extraction` in `patients/urls.py` |
| Search | `patients/signals/specimen_search.py` — `HAS_3_ANY`, since a TSO 500 reference is all digits; an extraction matches its own reference *or* its specimen's |
| Preview / hover cards | `Specimen.preview` (`patients/models.py:365`), `Extraction.preview` (`patients/models.py:453`) |
| Links down from the containing model | patient page Specimens / Extractions tabs; specimen page lists its extractions, their samples and its measures; extraction page lists samples and sequencing samples |
| Autocomplete | `SpecimenAutocompleteView`, `ExtractionAutocompleteView` (`patients/views_autocomplete.py`) |

Permissions are settled and both parts below inherit them: `Specimen` and `Extraction` take their
confidentiality from their `Patient` via `get_permission_class` / `_filter_from_permission_object_qs`,
so every queryset goes through `filter_for_user`. Samples are filtered separately, because a sample
carries its VCF's permissions rather than its patient's — `view_specimen` already does this with a
`Prefetch`, and the grids need the same care wherever they cross that line.

---

## The grid half

### A1. Top-level Specimen and Extraction grids

Nothing exists today: `patients/grids.py` has only the patient grids, and
`uicore/templates/uicore/menus/menu_bar_patients.html` has no entries.

**DataTables, not jQGrid** — new tables use `DatatableConfig` + `RichColumn`
(`snpdb/views/datatable_view.py`). `PatientRecordColumns` (`patients/grids.py:112`) is the closest
model to copy, including `renderer=self.view_primary_key` with
`client_renderer='TableFormat.linkUrl'` for the link column.

- `SpecimenColumns(DatatableConfig[Specimen])` and `ExtractionColumns(DatatableConfig[Extraction])`
  in `patients/grids.py`, each with `get_initial_queryset()` returning `filter_for_user(self.user)`.
- Listing views + templates alongside `patients`/`patients.html`, named `specimens` and `extractions`.
- URLs: the page plus a `DatabaseTableView.as_view(column_class=...)` datatables endpoint, following
  the `patient_records_datatables` pair in `patients/urls.py`.
- Menu: two `{% menu_item %}` entries in `menu_bar_patients.html`, under `patients`.
- Shariant has no patients menu, so add `"specimens": False` and `"extractions": False` to the
  `URLS_NAME_REGISTER` block at `variantgrid/settings/env/shariantcommon.py:185`, beside the
  `view_specimen` / `view_extraction` entries already there.

Columns worth having on the specimen grid: reference id (linked), patient (linked), `patient_code`,
tissue, tissue status, collection date, received date, external pk, extraction count. On the
extraction grid: reference id (linked), specimen (linked), patient, nucleic acid source, extraction
date, external pk, sample count. Follow `PatientListGrid`'s lead on identity — `patient_code` is the
de-identified column, so it is the one that belongs in a grid by default.

### A2. Links on the sample grid

`SamplesListGrid` (`snpdb/grids.py:115`) already pulls `extraction__specimen__reference_id`,
`extraction__specimen__tissue__name` and `extraction__specimen__collection_date`, but they render as
plain text — the grid has no extraction column at all, and nothing links anywhere.

- Add hidden `extraction__id` and `extraction__specimen__id` columns, then give the reference columns
  `'formatter': 'linkFormatter'` with `url_name` / `url_object_column` pointing at `view_specimen` and
  `view_extraction`. The `vcf__name` override a few lines up (`snpdb/grids.py:141`) is the pattern.
- Add an extraction reference column (`extraction__reference_id`), since the DNA and RNA arms of one
  block share a specimen and the specimen reference alone cannot tell those rows apart.

While in the neighbourhood: `snpdb/templates/snpdb/tags/related_data_for_patient.html:26` renders
`{{ sample.tissue }}`, and `Sample` has no `tissue` — it has a `specimen` property
(`snpdb/models/models_vcf.py:360`). That column has been rendering blank since Phase 1 moved the FK.
It wants to be the specimen, linked.

### A3. Specimen and patient IDs on the variant page's sample table

[SACGF/variantgrid_private#2837](https://github.com/SACGF/variantgrid_private/issues/2837) — *"be
able to see the specimen & patient IDs in the variant table on the variant page… it helps us work out
when the data come from the same patient/specimen"*. The question being asked is "which of these rows
are the same patient", so the columns want to be sortable and to read at a glance.

The table is DataTables built client-side from the `variant_sample_genotypes` API:

- `snpdb/variant_sample_information.py:33` — `COPY_SAMPLE_FIELDS` already fetches
  `patient__patient_code`, but not the specimen. Add the specimen (and extraction) reference and pks,
  and emit them in the row dict built at `snpdb/variant_sample_information.py:424`.
- `variantgrid/static_files/default_static/js/variant_sample_information.js:147` — add the columns.
  `renderSampleLink` (line 62) is the pattern for a linked cell; link the specimen at its
  `view_specimen` URL.

Both models are `filter_for_user`-gated on the patient, and this table is reachable by anyone who can
see the variant, so the specimen reference has to be filtered on the *patient's* permissions rather
than the sample's — a user who can see a sample cannot necessarily see its patient. Where the patient
is not visible, the row shows the sample and leaves the new columns blank, the same way the existing
phenotype columns behave.

---

## Tissue — moved to #1747

The issue's other loose end — "there's currently no way to create tissue, either" — turned out to be a
model and ontology change rather than a seed list, so it is
[#1747](https://github.com/SACGF/variantgrid/issues/1747), designed in
[`1747_specimen_tissue_ontology_plan.md`](1747_specimen_tissue_ontology_plan.md). Nothing here waits
on it.

What stays with #1706 is the visible half of it: the specimen page and the patient page's Specimens
tab both render a tissue `<select>` that can never have options, so `'tissue'` joins the excludes in
`SpecimenForm` (`patients/forms.py:137`) and `PatientSpecimenFormSet` (`patients/forms.py:122`) until
#1747 backs the field with something.

## Suggested order

1. **A2** — the smallest change, and the one that most directly answers "where did this sample come
   from".
2. **A3** — private#2837, a real user asking for it.
3. **A1** — the largest, and the one the issue itself was least sure was worth having ("think about
   whether it's worth…"). Doing A2 and A3 first is also the cheapest way to find out, since if the
   links off the sample and variant tables answer the question, the top-level grids may not need to be
   as rich as sketched above.

The two form excludes are a couple of lines and can ride with any of them.

## Done when

- Patients menu has Specimens and Extractions entries; both grids page, sort and filter, respect
  `filter_for_user`, and every reference cell links to its detail page. Shariant serves neither.
- The sample grid shows specimen and extraction references, both linked, and
  `related_data_for_patient.html` shows a linked specimen instead of a blank column.
- The variant page's sample table shows patient code and specimen reference, blank where the viewer
  cannot see the patient, and two rows from the same specimen are recognisable as such.
- Neither the specimen page nor the patient page's Specimens tab offers an empty tissue dropdown.
