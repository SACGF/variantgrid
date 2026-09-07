# #1574 — `<CNV>` has no HGVS: record it, stop reporting it

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-07
Status: in progress

[#1574](https://github.com/SACGF/variantgrid/issues/1574): classifying a `<CNV>` variant raises
`ValueError: Unknown symbolic alt of '<CNV>'` from `ResolvedVariantInfo.set_variant_and_save`, which
sends it to Rollbar once per genome build. The variant is fine; there is simply no HGVS for it.

---

## 1. Diagnosis (confirmed against the code, reproduced on vg-test2)

* #1570 (`d43e9907c`) taught the VEP inserter to skip `as_external_explicit()` for alts that cannot be
  expanded, via `snpdb/models/models_variant.py:VariantCoordinate.can_be_made_explicit`.
* #1571 (`f7fcb259f`) gave `<DEL>`/`<DUP>`/`<INV>` a coordinate-only HGVS through
  `snpdb/models/models_variant.py:VariantCoordinate.symbolic_hgvs_interval`.
* `<CNV>` and `<INS>` were deliberately left out of both: "copy number changed somewhere in this range"
  has no unambiguous HGVS. `symbolic_hgvs_interval` returns None for them, so
  `genes/hgvs/hgvs_matcher.py:HGVSMatcher.variant_coordinate_to_hgvs_used_converter_type_and_method`
  falls through to `as_external_explicit()` and hits the `raise ValueError(f"Unknown symbolic alt ...")`.
  Both the g.HGVS (no transcript) and c.HGVS (transcript) paths raise the same way.
* `classification/models/classification_variant_info_models.py:ResolvedVariantInfo.set_variant_and_save`
  treats `NoTranscript` as a data issue (store in `error`, log a warning) and everything else as a bug
  (`report_exc_info`). A missing HGVS form is the first kind, not the second.

`Variant.can_make_g_hgvs` (`snpdb/models/models_variant.py:Variant`) hard-codes the same
`{DEL, DUP, INV}` set as `VariantCoordinate.can_be_made_explicit`; the template tag
`annotation/templatetags/clinvar_tags.py` is its only caller.

## 2. Change

No model or migration. Four edits.

### 2.1 A typed exception — `genes/hgvs/hgvs_converter.py`

```python
class HGVSNoRepresentationException(HGVSException):
    """ The variant is valid but HGVS has no way to write it - <CNV>, <INS> and other
        symbolic alts with neither a ranged form nor an explicit ref/alt expansion """
```

Export it from `genes/hgvs/__init__.py` next to the other `HGVSException` subclasses.

### 2.2 Raise it before any converter runs — `genes/hgvs/hgvs_matcher.py`

At the top of `variant_coordinate_to_hgvs_used_converter_type_and_method`, before the transcript loop:

```python
if variant_coordinate.symbolic_hgvs_interval is None and not variant_coordinate.can_be_made_explicit:
    raise HGVSNoRepresentationException(f"{variant_coordinate.alt} has no HGVS representation")
```

This is the one place both the g. and c. paths pass through, and it runs before the ClinGen Allele
Registry branch would be attempted for the coordinate. `variant_coordinate_to_g_hgvs` calls
`as_external_explicit()` itself first, so give it the same guard so it raises the typed exception
rather than the bare ValueError.

### 2.3 Treat it as data, not a bug — `classification/models/classification_variant_info_models.py`

In `ResolvedVariantInfo.set_variant_and_save`, catch `HGVSNoRepresentationException` alongside
`NoTranscript`: set `self.error = str(exception)` and `logging.info(...)`; leave `report_exc_info` to
the generic branch. Import it at the top of the file with the other `genes.hgvs` names.

The classification form already renders `ResolvedVariantInfo.error` in the c.HGVS field, so the user
sees "<CNV> has no HGVS representation" in place of a blank.

### 2.4 One source of truth for "expandable symbolic alt" — `snpdb/models/models_variant.py`

Make `Variant.can_make_g_hgvs` delegate to `self.coordinate.can_be_made_explicit` so the set of
expandable alts lives in one place.

## 3. Tests

* `genes/tests/test_hgvs.py:TestSymbolicHGVS` — add one test: a `<CNV>` coordinate raises
  `HGVSNoRepresentationException` from both `variant_coordinate_to_g_hgvs` and
  `variant_coordinate_to_hgvs_variant` with a transcript. The existing
  `test_symbolic_hgvs_interval_none_without_ranged_form` stays as is.
* `classification/tests/models/test_imported_allele_info.py` — add one test that
  `ResolvedVariantInfo.set_variant_and_save` on a `<CNV>` variant stores the message in `error`,
  leaves `c_hgvs` None, and calls `report_exc_info` zero times (patch `classification.models.classification_variant_info_models.report_exc_info`).

`scripts/vg tests --explain` will list the rest; run them with `--keepdb`.

## 4. Backfill

None needed. `ResolvedVariantInfo` rows already carrying `error = "Unknown symbolic alt of '<CNV>'"`
display correctly today; the new wording only applies from the next resolve. vg-test2 has zero `<CNV>`
variants and zero affected rows.

## 5. Docs

* `genes/CLAUDE.md`: one line under HGVS gotchas — `<CNV>`/`<INS>` have no HGVS; the matcher raises
  `HGVSNoRepresentationException` and classification records it as `ResolvedVariantInfo.error`.
* `scripts/vg docs check` after the edit.
