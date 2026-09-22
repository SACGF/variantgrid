# Somatic classification: fill the fields a new record always starts by complaining about

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-22
Status: in progress - PR #1898

[sapath#453](https://github.com/SACGF/variantgrid_sapath/issues/453). A classification made from the Classify & Report
tab (and the create page) opens with four messages a scientist has to clear by hand every time:

| Message | Source | Why it appears |
|---|---|---|
| Classification requires a value | `classification/models/classification_variant_fields_validation.py:validate_variant_classification_significance` | Neither `clinical_significance` nor `somatic:clinical_significance` set. Wording is wrong for somatic ([#1893](https://github.com/SACGF/variantgrid/issues/1893), separate, Shariant testing) |
| Curated/reviewed by - Missing mandatory value | `Classification.patch_value` mandatory check | SA Path's lab/org `classification_config` marks `curated_by` mandatory; nothing fills it |
| Curation/review date - Missing mandatory value | as above | as above, for `curation_date` |
| Level D indicates Somatic Clinical Significance should be "Tier II" but it is blank | `classification_variant_fields_validation.py:validate_letter_to_tier` | Copy consensus brings the AMP level keys (`copy_scope=ALLELE`) but `somatic:clinical_significance` is `copy_scope=NONE`, so the tier is left behind |

Two changes, both applied once when the record is created and never again, so the scientist can change what
they disagree with:

1. **Dynamic defaults in the lab config** - `default_value` in `classification_config` learns two tokens, `$user`
   and `$today`, so the lab that makes a field mandatory says what it starts as.
2. **Tier derived from the AMP level** - after the consensus copy, an empty `somatic:clinical_significance` is set
   from the highest level present, using the same mapping the validator checks.

## Data

No model changes. The configuration lives in JSON already in the database:

`snpdb/models/models.py:Organization.classification_config` / `snpdb/models/models.py:Lab.classification_config`

```json
{
  "curated_by":    {"mandatory": true, "default_value": "$user"},
  "curation_date": {"mandatory": true, "default_value": "$today"}
}
```

`classification/models/evidence_key.py:EvidenceKeyOverrides.from_dict` accepts any EvidenceKey attribute per key and
`classification/models/evidence_key.py:EvidenceKeyMap` applies it with `setattr`, so `mandatory` and `default_value`
already reach the key. `EvidenceKey.default_value` is an instance attribute that only an override sets (it is not a
column), and `classification/models/classification.py:Classification.create_with_response` copies it into a new record
when `populate_with_defaults=True`. The one caller that passes that is the web create path
(`classification/views/views.py:create_classification_object`) - the API, file imports and sync never see a default, so
Shariant records are untouched.

SA Path applies the JSON above to the labs that want it (admin, Lab or Organization). That is the deploy step; there is
no migration.

## 1. Dynamic defaults

`classification/models/classification.py:Classification.create_with_response`, in the `populate_with_defaults` loop:

```python
if populate_with_defaults:
    for e_key in record.evidence_keys.all_keys:
        if e_key.default_value is not None and e_key.key not in data:
            data[e_key.key] = resolve_default_value(e_key.default_value, user)
```

`resolve_default_value` lives in `classification/models/evidence_key.py` beside `EvidenceKeyOverrides`:

| Token | Resolves to |
|---|---|
| `$user` | `user.get_full_name()`, falling back to `user.username` (the same choice `library/vg/inspect/user.py` makes) |
| `$today` | `date.today().isoformat()` - the `YYYY-MM-DD` string a `D` value type stores |
| anything else | returned unchanged, as today |

The record is created with the user, so no signature change. Values arrive as `SubmissionSource.VARIANT_GRID` like
every other initial value, editable in the form.

Document the tokens in the `EvidenceKeyOverrides.from_dict` docstring, which is the only description of the config
format there is.

## 2. Tier from the AMP level

`classification/autopopulate_evidence_keys/autopopulate_evidence_keys.py:classification_complete_web_create`, after the
copy loop:

```python
for source, copy_scopes in [(copy_from, COPY_SCOPES_ALL), (copy_gene_from, COPY_SCOPES_GENE)]:
    if source:
        ClassificationConsensus(modification=source, copy_scopes=copy_scopes).apply_to(classification, user)
apply_somatic_tier_from_amp_level(classification, user)
```

`apply_somatic_tier_from_amp_level` goes in `classification/models/classification_variant_fields_validation.py` next to
`validate_letter_to_tier`, because the level-to-tier table (`__LEVELS_TO_TIER`) it needs is there; the module-private
names become module-level constants the function and the validator share. It:

- returns without a write when `somatic:clinical_significance` already has a value, or no `amp:level_*` key does;
- takes the highest level present (A before B before C before D, the order `SpecialEKeys.AMP_LEVELS_TO_LEVEL` iterates);
- maps it through the table to a tier ("1" or "2") and picks the first option of `somatic:clinical_significance` whose
  `tier` attribute equals it (`tier_1` before `tier_1_or_2`, the option order in the key), via
  `EvidenceKey.option_dictionary_property("tier")`;
- patches that one value with `source=SubmissionSource.CONSENSUS`, `leave_existing_values=True`, then
  `publish_latest`, the way `ClassificationConsensus.apply_to` does.

Running it after both copies means a gene-level copy that carried the levels also gets its tier. A record with no copy
and no levels is left alone: the validator has nothing to say about it either.

The validator itself is unchanged. It still warns when a scientist later moves a level without moving the tier, which
is the point of it.

## Not in scope

- Refreshing `curation_date` on later saves or on publish. It is ClinVar's "date last evaluated", the reclassification
  timeline reads it and discordance ordering compares it, and every form field edit is a save. The initial value is
  the default; changing it afterwards is the curator's.
- The "Classification requires a value" wording - #1893.

## Tests

`classification/tests/models/test_classification_consensus.py` already covers the copy scopes. Add:

- `resolve_default_value`: `$user` with and without a full name, `$today`, a literal - one test method.
- `create_with_response(populate_with_defaults=True)` under a lab whose `classification_config` sets a `$user` default
  fills the key; the same call without `populate_with_defaults` leaves it empty.
- `apply_somatic_tier_from_amp_level`: Level D alone gives `tier_2`; Levels A and D give `tier_1`; an existing tier
  is kept; no levels writes nothing (no new ClassificationModification).

## Files

- `classification/models/evidence_key.py` - `resolve_default_value`, docstring for the tokens
- `classification/models/classification.py` - `create_with_response` resolves defaults
- `classification/models/classification_variant_fields_validation.py` - `apply_somatic_tier_from_amp_level`, shared constants
- `classification/autopopulate_evidence_keys/autopopulate_evidence_keys.py` - call it after the copies
- `classification/tests/models/test_classification_consensus.py` - tests above
- `classification/AGENTS.md` - one line under Patterns: lab config `default_value` accepts `$user` / `$today`, applied at web create only
