# #1848 — AnnotSV ACMG class as a categorical enum

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-09
Status: landed e864c2869

[#1848](https://github.com/SACGF/variantgrid/issues/1848): `VariantAnnotation.annotsv_acmg_class` is a bare
`IntegerField` holding AnnotSV's `ACMG_class` (1..5). The variant details page maps it to a clinical significance
pill by hand, but the analysis grid shows the raw integer and the column summary treats it as a quantity and draws
a box plot instead of counting the five classes.

## The data

No new model and no column change. The field keeps its type and its stored values; it gains `choices` from an
enum that already exists with exactly AnnotSV's scale.

```python
# annotation/models/models_enums.py  (rename of ClinVarPathogenicity - same members, same values)
class Pathogenicity(models.IntegerChoices):
    """ The five-tier ACMG germline scale as an ordered integer: ClinVar.highest_pathogenicity (CLNSIG
        mapped onto it, 0 = none of the below) and VariantAnnotation.annotsv_acmg_class (AnnotSV's ACMG_class) """
    BENIGN = 1, "Benign"
    LIKELY_BENIGN = 2, "Likely benign"
    UNCERTAIN = 3, "Uncertain"
    LIKELY_PATHOGENIC = 4, "Likely pathogenic"
    PATHOGENIC = 5, "Pathogenic"
```

```python
# annotation/models/models.py  VariantAnnotation
annotsv_acmg_class = models.IntegerField(null=True, blank=True, choices=Pathogenicity.choices)
```

`ClinVarOncogenicity` stays as it is. The `ClinVar.highest_pathogenicity` field is out of scope: its `default=0`
is not in the enum and the classifications chip in `variantgrid/static_files/default_static/js/variantgrid_formats.js`
reads it as an integer, so giving it `choices` is a separate change (see Follow-ups).

Why this enum and not `classification/enums/classification_enums.py:ClinicalSignificance`: that one holds
`CharField` values `'0'..'5'` for classification records, so reusing it would mean an integer-to-text data migration
across `annotation_variantannotation` for no gain. `ClinVarPathogenicity` is an `IntegerChoices` whose values are
AnnotSV's values, so the change is `choices=` and nothing in the database moves.

## How the grid gets it from there

`snpdb/grid_columns/custom_columns.py:_model_field_column_kwargs` already gives any field with `choices` a
`select` column filter, which is what `analysis/views/views_node.py:node_column_summary` uses to decide categorical
(counts grid) versus quantitative (box plot). It only attaches the label renderer to a `CharField` with choices,
though, so an integer choice field would reach the client and the CSV as its number. Extending that renderer to
every field with `choices` is the one code change the grid needs; the composite cell and the summary grid both
read through it (`analysis/grids.py:NodeColumnSummaryConfig` reuses `rich_column.renderer` for its labels).

The only other `IntegerField(choices=...)` in the codebase are on `ClassificationGrouping`, which is never a
variant grid column, so nothing else changes behaviour.

## Steps

### §1 Rename the enum

- `annotation/models/models_enums.py`: rename `ClinVarPathogenicity` to `Pathogenicity` with the docstring above.
- Update the references: `analysis/models/nodes/filters/clinvar_node.py`, `analysis/models/nodes/node_counts.py`,
  `analysis/tests/test_clinvar_node.py`. `ClinVarOncogenicity` imports alongside stay as they are.
- Add `short_label` and `css_class` properties on the enum, and build `PATHOGENICITY_SHORT_LABELS` /
  `PATHOGENICITY_CSS_CLASSES` in `clinvar_node.py` from them plus the `NO_CLINVAR_CALL` entry, so the short forms
  (B, LB, VUS, LP, P) and the `.c-pill.cs-*` class names live in one place. `PATHOGENICITY_LABELS` becomes
  `{p: p.label for p in Pathogenicity} | {NO_CLINVAR_CALL: "Other"}`. The short labels are also the
  `clinical_significance` evidence key option values, which is what §3 relies on.

### §2 The field

- `annotation/models/models.py`: `annotsv_acmg_class = models.IntegerField(null=True, blank=True, choices=Pathogenicity.choices)`.
  Delete `ANNOTSV_ACMG_CLASS_CLINICAL_SIGNIFICANCE`.
- Migration annotation/migrations/0182_annotsv_acmg_class_choices.py (new): the `AlterField` makemigrations produces.
  Choices are not enforced by Postgres so this is metadata only and runs instantly on the 175 GB box.
- `annotation/vcf_files/bulk_annotsv_tsv_inserter.py:_parse_value`: an `ACMG_class` outside the enum's values is
  stored as `None`, the same as `NA`. AnnotSV documents 1..5 only, so this is a guard, and the comment on that
  line already says what the field is.

### §3 The variant details page

`variantopedia/templates/variantopedia/variant_details.html` renders the class through the
`clinical_significance` inclusion tag, which wants an evidence key option value. Keep that: the
`annotsv_acmg_clinical_significance` property on `VariantAnnotation` returns `Pathogenicity(value).short_label`
(or `None`), replacing the hand-written dict lookup. The page looks the same before and after.

### §4 The grid

- `snpdb/grid_columns/custom_columns.py:_model_field_column_kwargs`: attach `_make_choices_renderer` and
  `csv_rendered = True` for any field with `choices`, not only `CharField`. The composite `annotsv_acmg` headline
  then reads "Likely pathogenic", the CSV export writes the label, and the column summary counts per label.
- Draw the headline as a chip. Add to `snpdb/grids.py:get_standard_overrides` (where composite member entries in
  `snpdb/grid_columns/custom_columns.py:_composite_column_kwargs` read their `client_renderer` from) an override for `variantannotation__annotsv_acmg_class` naming `VariantGridFormat.pathogenicityChip`, and add
  that renderer to `variantgrid/static_files/default_static/js/variantgrid_formats.js`: the same `cs-chip cs-<short>`
  span the ClinVar and internal classification columns draw, showing the abbreviation (B, LB, VUS, LP, P) and keyed
  by the label the row now carries, in the style of `CLINVAR_SOMATIC_TIER_CHIPS` (a choice field's row carries the
  label, so the lookup is by label). The composite's hover title still uses the plain label text, and the chip
  leaves the title to the cell so the score and criteria stay reachable.
- Migration snpdb/migrations/0260_annotsv_acmg_class_description.py (new): update the `annotsv_acmg_class`
  `VariantGridColumn.description` so it names the classes as words (Benign … Pathogenic) rather than the 1..5 code
  table that `snpdb/migrations/0180_new_annotsv_variantgrid_columns.py` wrote. The annotation descriptions page
  draws the example cell from `snpdb/grid_columns/composite_examples.py`, whose `annotsv_acmg_class=4` needs no change.

### §5 Tests

Keep only tests of our branches:

- `annotation/tests/test_annotsv.py`: one case in `TestRowToUpdate` that an out-of-range `ACMG_class` (e.g. `"7"`)
  is absent from the update. The existing assertions that `"5"` becomes `5` still hold and stay.
- A test on `variant_column_rich_column("variantannotation__annotsv_acmg_class")` asserting the column filter
  type is `select` and the renderer maps `4` to "Likely pathogenic" - that is the branch §4 adds and the condition
  `node_column_summary` keys off. Put it in `snpdb/tests/test_composite_columns.py`, which already builds variant grid columns
  through `snpdb/grid_columns/custom_columns.py`.
- `analysis/tests/test_clinvar_node.py` keeps passing under the rename; no new test for the rename itself.

### §6 Docs and maps

- `scripts/vg map` after the enum rename and the field change; commit `claude/maps/*.md`.
- One line under `uicore/CLAUDE.md#grids`: a field with `choices` (any type) reaches the grid, CSV and column
  summary as its label, so a client renderer for it keys by label, not by stored value.
- `scripts/vg docs check` on this plan and the edited notes.

## Verification on vg-test2

- `vg page` a variant details page for an SV with AnnotSV annotation, before and after §3, and compare the pill.
- In an analysis with an SV VCF, add the AnnotSV ACMG column, confirm the cell reads as a pill, that the header sort
  still orders by class, and that "column summary" on it now shows a counts grid with the five labels instead of a
  box plot. Click a label to confirm the filter child node is created with the stored integer.
- Export the node as CSV and check the column carries the label.

## Follow-ups (out of scope)

- `ClinVar.highest_pathogenicity` has the same column summary problem (box plot over 0..5). Fixing it means either
  admitting 0 into the enum or giving the field `null=True`, and changing the classifications chip in
  `variantgrid/static_files/default_static/js/variantgrid_formats.js` to key by label. Separate issue.
