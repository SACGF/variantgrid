# Classification: reassign the sample (#1392)

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-11-09
Status: in progress

A classification is linked to a `Sample` at creation and the link cannot be changed afterwards, so a record created
against the wrong sample (the Classify & Report dialog lists samples the variant was not called in, #2737) stays
wrong. This adds a Sample card to the classification form from which a lab head reassigns the record to another
sample and, by default, re-populates the sample-derived evidence from it. Lab reassignment is a separate issue.

## Data

No new model. `Classification.sample` (`classification/models/classification.py:Classification`, nullable FK to
`Sample`) is re-pointed, and the evidence change is an ordinary `ClassificationModification` written through
`Classification.patch_value`, so it versions, diffs and publishes like any form edit.

```python
@dataclass(frozen=True)
class SampleEvidenceChange:
    key: str            # evidence key
    label: str          # EvidenceKey.pretty label for the dialog
    current: Any        # value on the record now (pretty_value), None when empty
    new: Any            # value the new sample supplies, None when it will be cleared


@dataclass(frozen=True)
class SampleReassignmentPreview:
    sample: Sample
    has_genotype: bool                           # sample.get_genotype(variant) found a call
    changes: list[SampleEvidenceChange]          # only the keys whose value would change
    allele_origin_bucket_before: AlleleOriginBucket
    allele_origin_bucket_after: AlleleOriginBucket
```

## Who may do it

`Classification.can_reassign_sample(user)`: the user is a superuser, or the record's lab is internal
(`Lab.external` False) and the user is one of its heads (`snpdb/models/models.py:Lab.is_head`, the `LabHead` rows
from #1470). The record also needs a resolved `variant`, since the genotype fields come from it. Everyone else sees
the Sample card read-only (the current sample, as today's Links card shows it) and no Change button.

## The sample-derived keys

`get_evidence_fields_for_sample_and_patient`
(`classification/autopopulate_evidence_keys/evidence_from_sample_and_patient.py`) decides what a sample populates.
Add a module constant there, `SAMPLE_EVIDENCE_KEYS`, listing every key that function and the two it calls
(`get_evidence_fields_for_extraction`, `get_evidence_fields_for_patient`) can set: sample id, allele origin, patient
and family id, specimen id, sample type, age, nucleic acid source, zygosity, allele depth, read depth, phred
likelihood, genotype quality, allele frequency, VCF filter, fold change, copy number, capture method, sequencing
platform, sex, ancestry. The re-populate step is defined over this set:

- a key the new sample supplies is overwritten with its value;
- a key in the set the new sample does not supply is cleared, so nothing from the old sample survives;
- allele origin is the exception: `get_allele_origin` deliberately returns None for a sample whose origin is
  per-variant (mixed tumour), leaving it to the curator, so allele origin is overwritten only when the new sample
  asserts one and otherwise kept.

The patch goes through `patch_value(patch, user=user, source=SubmissionSource.VARIANT_GRID,
leave_existing_values=False, save=True)`. VARIANT_GRID is the highest immutability level, so it overwrites values a
curator typed, which is the point: the preview is where they see what goes.

## Module: classification/sample_reassignment.py (new file)

Owns the two dataclasses above and:

- `sample_reassignment_preview(classification, sample, user) -> SampleReassignmentPreview`. Checks
  `Sample.check_can_view(user)`, builds the new values with `get_evidence_fields_for_sample_and_patient(variant, sample)`,
  diffs them against `classification.get(key)` over `SAMPLE_EVIDENCE_KEYS` (pretty values via
  `EvidenceKeyMap.cached()`), and computes the bucket after from the new allele origin with
  `EvidenceKeyMap.allele_origin_bucket`. `has_genotype` is `sample.get_genotype(variant) is not None`.
- `reassign_sample(classification, sample, user, repopulate: bool) -> ClassificationPatchResponse | None`. Permission
  check, sets `classification.sample`, saves, and when `repopulate` applies the patch described above. One
  `transaction.atomic`.

`get_evidence_fields_for_sample_and_patient` logs an error when the sample has no genotype for the variant; that
message is for the create path and stays.

## View: classification/views/views_classification_sample.py (new file)

One view, `classification_sample_panel(request, classification_id)`, following
`classification/views/views_gene_consensus.py:gene_consensus_panel` (the card, its dialog content and the POST on one
URL). URL `classification/<int:classification_id>/sample`, name `classification_sample`, in `classification/urls.py`.

- GET: renders the card body, classification/templates/classification/sample_card.html (new file): the current sample
  with its link (or "No sample"), and for a user who may reassign, the Change button.
- GET with `sample_id`: renders the preview partial, classification/templates/classification/sample_change_preview.html
  (new file): a warning line when `has_genotype` is False ("This sample has no genotype call for this variant. The
  genotype fields will be cleared."), an alert when the allele origin bucket changes (before and after), then one row
  per change with current and new value, "cleared" where new is None, and a line saying nothing changes when the
  list is empty. Any `PermissionDenied` from the sample check renders as the dialog's error line.
- POST `sample_id`, `repopulate` (checkbox), `confirm_no_genotype` (the override tick, required when the preview had
  no genotype): calls `reassign_sample`, returns JSON `{"reassigned": true}`. The page reloads.

Permission on the view is `classification.can_write(request.user)` for the card and `can_reassign_sample` for the
preview and the POST.

## Page

`classification/templates/classification/classification.html`:

- a `#vcform-sample` card in `#vc-extras` next to the Gene Content card, body loaded from the panel URL, when the view
  passes `sample_panel_enabled` (`settings.CLASSIFICATION_GRID_SHOW_SAMPLE`-independent: the card shows whenever the
  record has a variant and the user can write it);
- a `#sample-change-modal` holding a `SampleChoiceForm` (`snpdb/forms.py:SampleChoiceForm`, genome build forwarded
  from the page's `genome_build` so the autocomplete lists samples of the right build), the "Re-populate sample
  fields from the new sample" checkbox (checked), a preview `<div>`, the no-genotype override checkbox (hidden
  until the preview needs it) and the Confirm button, disabled until a preview has loaded and, with no genotype,
  until the override is ticked.

The autocomplete starts scoped to the current sample's patient: `snpdb/views/views_autocomplete.py:SampleAutocompleteView`
already filters on a forwarded `patient` (both the sample's own patient FK and the specimen chain), so the modal
carries a hidden `patient` input the form forwards, set to `classification.sample.patient_id` when the record has a
sample with a patient. A line under the field says "Showing samples for this patient" with a link "Show all samples",
which clears the hidden input and re-initialises the select2 so the next search is over every sample the user can
see. Without a patient the field starts unscoped and the line is absent. The wrong sample is usually a relative in
the same case (#2737), and a patient filter is a foreign key lookup; anything that asks which samples carry the
variant is a scan of the variant's genotypes and is far too slow for a dialog.

The JS lives in variantgrid/static_files/default_static/js/classification/classification_sample.js (new file):
load the card, open the modal, on sample select fetch the preview (`?sample_id=`), toggle the override and Confirm,
POST on confirm and reload. The template keeps a `<script src>` and one init call with the panel URL
(`variantgrid/static_files/AGENTS.md`).

`view_classification` in `classification/views/views.py` adds `sample_panel_enabled` to the context.

## Removing the old path

`classification/models/classification_evidence_normalise.py` has an `existing_sample_id` branch in
`normalise_evidence_cell` that sets `classification.sample` from a patch value. It is not a registered evidence key
and its re-populate line is commented out; this feature replaces it, so the branch goes, with the `Sample` import if
nothing else in the module uses it.

## Tests

classification/tests/test_sample_reassignment.py (new file), on the module rather than the view. Fixtures:
`ClassificationTestUtils` for lab and users (`claude/guides/testing.md`), `LabHead.objects.create` for the head,
`snpdb/fake_data.py:create_fake_cohort` for samples and `make_cohort_genotype` for a genotype on the record's
variant, `Classification.create` with sample-derived evidence already on the record.

- permission: head of the owning lab may, a plain lab member may not, superuser may, head of an external lab may not;
- preview: an overwritten key, a cleared key, allele origin kept when the new sample is mixed, the bucket change
  reported when the new sample is somatic-only, `has_genotype` False for a sample with no call;
- apply: the record points at the new sample, a new `ClassificationModification` carries the new values, the
  cleared keys are None, and with `repopulate=False` the evidence is untouched.

## Docs and close-out

- One line in `classification/AGENTS.md` under "Start with" naming the module and that it is the one way a record
  changes sample.
- A changelog entry in `variantgrid/templates/default_templates/changelog.html` following its header rules.
- Testing notes on #1392 (VariantGrid testing only; Shariant has no samples).
