# Plan: Let g.HGVS-only classifications through (#1063)

**Issue**: A classification imported as g.HGVS (e.g. `NC_000009.12:g.99087902A>C`) has no transcript and no c.HGVS. `ImportedAlleleInfo._calculate_validation()` validates it as though a c.HGVS were expected, producing errors that set `include=False` and block export, plus bad UI links.

**Reported symptoms**:
- `General Transcript Type Not Supported` ERROR
- `Liftover c.nomen Change` ERROR
- `Liftover Transcript Version Change` ERROR
- `Builds Missing 37` / `Builds Missing 38` WARNINGs
- HGVS resolution tool link renders `?hgvs=None` (errors when visited)
- No allele link on the classification page (despite the allele existing)
- Record excluded from exports (`include=False`)

---

## Design decision

A g.HGVS-only submission is accepted on its own terms. The submitter may have a record of the genomic coordinate and not the gene or transcript; or there may have been no transcript at that position when they classified it, even though there is one now. Either way the record goes through.

Validation therefore keys off **what was supplied on import**, not off what the resolver later managed to find:

- **c.HGVS or transcript supplied on import** → validate exactly as today. A missing or unsupported c.HGVS is a genuine resolution failure and stays an error.
- **g.HGVS only** → the record is accepted. Where resolution does yield a transcript-based c.HGVS, surface the difference as a **warning** so it is visible without blocking.

This keeps the resolver's outcome as information for the curator, and keeps the accept/reject decision tied to the submission.

---

## Why each error fires today

**1. `transcript_type_not_supported` ERROR** — `classification_variant_info_models.py:601`
```python
if not ImportedAlleleInfo.is_supported_transcript(self.get_transcript):
    general["transcript_type_not_supported"] = ...
```
`get_transcript` (`:675`) returns `imported_transcript`, else the transcript parsed out of `imported_c_hgvs` — so `None` for a g.HGVS-only import. `is_supported_transcript(None)` returns `False`. The check exists to reject unsupported transcript prefixes (e.g. `NX_`); it fires on "no transcript at all" as a side effect.

**2. `Liftover c.nomen Change` + `Liftover Transcript Version Change`** — `:587-589`
`ResolvedVariantInfo.recalc_c_hgvs()` (`:207-232`) calls the matcher with `imported_transcript=None`, which falls back to the genomic representation, and that g.HGVS is stored in the `c_hgvs` field. `_calculate_validation()` then diffs the GRCh37 and GRCh38 values. The contig accessions differ by build (`NC_000009.11` vs `NC_000009.12`), so the diff always reports `transcript_version_change` (W) and `c_nomen_change` (E). The `c_nomen_change` error is what sets `include=False`.

**3. `Builds Missing 37` / `Builds Missing 38`** — `:591-595`
The check requires `c_hgvs_obj` on each `ResolvedVariantInfo`. For the partially-overlapping deletion in the issue comments (`NC_000009.12:g.99044538_99087902del`), liftover succeeds and both `ResolvedVariantInfo` rows exist, but `recalc_c_hgvs()` errors because the variant runs off the transcript, leaving `c_hgvs` empty. Both severities are `W` (`:278-279`), so these are cosmetic rather than export-blocking.

**4. `?hgvs=None` link** — `imported_allele_info_detail.html:60` and `imported_allele_info_check.py:113`
Both build the resolution-tool URL from `imported_c_hgvs`; Django's template engine renders `None` as the literal string `"None"`, and `urlencode()` stringifies it the same way.

---

## Step 1: Add `imported_as_c_hgvs` to `ImportedAlleleInfo`

**File**: `classification/models/classification_variant_info_models.py`, near `get_transcript` (`:675`)

```python
@property
def imported_as_c_hgvs(self) -> bool:
    """True when the submission supplied a transcript-based HGVS, so a resolved c.HGVS is expected."""
    return bool(self.imported_c_hgvs or self.imported_transcript)
```

Read the fields directly rather than going via `get_transcript`, so a malformed `imported_c_hgvs` that parses to no transcript still counts as a c.HGVS submission and keeps its errors.

---

## Step 2: Key validation off the submission

**File**: `classification/models/classification_variant_info_models.py`, `_calculate_validation()` (`:565-609`)

### 2a. Gate `transcript_type_not_supported`
```python
if self.imported_as_c_hgvs and not ImportedAlleleInfo.is_supported_transcript(self.get_transcript):
    general["transcript_type_not_supported"] = _VALIDATION_TO_SEVERITY.get("transcript_type_not_supported", "E")
```

### 2b. Run the normalize diff for c.HGVS submissions
The `normalize` diff compares `imported_c_hgvs` against the resolved value, so it applies where a c.HGVS was imported. The existing `if imported_c_hgvs and normalized_c_hgvs:` guard already covers this — leave it as is.

### 2c. Report the liftover diff as a warning for g.HGVS-only submissions
Give `calculate_diff_dict` an optional severity override:
```python
def calculate_diff_dict(c_hgvs_diff: CHGVSDiff, severity: Optional[ALLELE_INFO_VALIDATION_SEVERITY] = None) -> ImportedAlleleValidationTagsDiff:
    diff_dict: ImportedAlleleValidationTagsDiff = {}
    for diff_flag, field_name in _DIFF_TO_VALIDATION_KEY.items():
        if c_hgvs_diff & diff_flag:
            diff_dict[field_name] = severity or _VALIDATION_TO_SEVERITY.get(field_name, "E")
    return diff_dict
```

Run the liftover diff when both builds resolved to a genuine transcript, and downgrade to `"W"` when the submission was g.HGVS-only:
```python
both_builds_have_transcript = all(
    build_info and build_info.transcript_version_id
    for build_info in (self.variant_info_for_imported_genome_build,
                       self.variant_info_for_lifted_over_genome_build)
)
if normalized_c_hgvs and lifted_c_hgvs and (self.imported_as_c_hgvs or both_builds_have_transcript):
    diff_severity = None if self.imported_as_c_hgvs else "W"
    if lifted_diff_dict := calculate_diff_dict(normalized_c_hgvs.diff(lifted_c_hgvs), diff_severity):
        validation_dict["liftover"] = lifted_diff_dict
```

`transcript_version_id` distinguishes a real c.HGVS from the g.HGVS fallback: `recalc_c_hgvs()` sets `transcript_version` from `c_hgvs_obj.transcript_version_model()` (`:219`), which resolves for `NM_`/`ENST` accessions and stays null for the `NC_` genomic form. So a g.HGVS-only record that resolves to a transcript in both builds gets its liftover difference reported as a warning, and one that resolves to genomic form in both builds gets nothing.

### 2d. Check build coverage against the variant for g.HGVS-only submissions
```python
builds: ImportedAlleleValidationTagsBuilds = {}
for build_str, build_info in (("37", self.grch37), ("38", self.grch38)):
    if self.imported_as_c_hgvs:
        resolved = build_info and build_info.c_hgvs_obj
    else:
        resolved = build_info and build_info.variant_id
    if not resolved:
        builds[f"missing_{build_str}"] = _VALIDATION_TO_SEVERITY.get(f"missing_{build_str}", "E")
```
For a g.HGVS-only submission the variant coordinate is the thing that has to exist in each build; the c.HGVS is a bonus.

---

## Step 3: Fall back to g.HGVS in the resolution-tool links

**File**: `classification/templates/classification/imported_allele_info_detail.html:60`
```html
<a class="hover-link" href="{% url 'hgvs_resolution_tool' %}?genome_build={{ allele_info.imported_genome_build|urlencode }}&hgvs={{ allele_info.imported_c_hgvs|default:allele_info.imported_g_hgvs|urlencode }}">Click here to test</a>
```

**File**: `classification/management/commands/imported_allele_info_check.py:113` (`_compare_url`)
```python
"hgvs": self.imported_allele_info.imported_c_hgvs or self.imported_allele_info.imported_g_hgvs
```

---

## Step 4: Allele link on the classification page

**Files**: `classification/templates/classification/classification.html:222`, `classification/views/classification_view.py`

The template renders the link under `{% if record == 'view_allele' and value.id %}`. A search of `classification_view.py` found no gating on `latest_validation.include`, so this symptom is likely downstream of `include=False` and should clear once Steps 1–2 land. Confirm against a live g.HGVS-only record after the fix; if the link is still absent, trace where the `view_allele` record value is populated and make it depend on the allele existing.

---

## Step 5: ClinVar export message

**File**: `classification/models/clinvar_export_convertor.py:518`

ClinVar requires a transcript-based HGVS, so excluding these submissions there is correct. Make the message say why:
```python
if allele_info and not allele_info.imported_as_c_hgvs:
    return ValidatedJson(None, JsonMessages.error(
        f"No c.HGVS available for this variant in {genome_build} — ClinVar requires a transcript-based HGVS."))
return ValidatedJson(None, JsonMessages.error(f"No normalised c.hgvs in genome build {genome_build}"))
```

---

## Step 6: Re-validate existing records

Existing g.HGVS-only records carry the old validation tags and stay excluded until re-validated. Add a management command, and register it for deployments via a `ManualOperation` migration with a `test=` callable that checks for affected rows:

```python
for allele_info in ImportedAlleleInfo.objects.filter(imported_c_hgvs__isnull=True, imported_transcript__isnull=True):
    allele_info.apply_validation(force_update=True)
    allele_info.save()
```

`apply_validation` preserves `include` on records a curator has confirmed (`:638`), so manual decisions survive the sweep.

---

## Files to change

| File | Change |
|------|--------|
| `classification/models/classification_variant_info_models.py` | Add `imported_as_c_hgvs`; Steps 2a, 2c, 2d in `_calculate_validation()` |
| `classification/templates/classification/imported_allele_info_detail.html` | g.HGVS fallback in resolution-tool link |
| `classification/management/commands/imported_allele_info_check.py` | g.HGVS fallback in `_compare_url` |
| `classification/models/clinvar_export_convertor.py` | Clearer exclusion message |
| `classification/views/classification_view.py` | Confirm allele link appears (Step 4) |
| New management command + `ManualOperation` migration | Re-validate existing g.HGVS-only records |

---

## Display check

`Classification.c_hgvs_best()` (`classification/models/classification.py:2175-2200`) falls through to `is_normalised = False` when `chgvs_grch37`/`chgvs_grch38` are empty, which is what renders *"not resolved, showing imported GRCh38.p14"* in the c.HGVS column (`vc_form.js:2392`). After the fix, check both cases from the issue:

- `NC_000009.12:g.99087902A>C` — resolves to the genomic form, so `chgvs_grch38` holds the g.HGVS and the column should show it.
- `NC_000009.12:g.99044538_99087902del` — runs off the transcript, so `recalc_c_hgvs()` errors and both cached columns stay empty. Confirm what the column shows and whether falling back to `imported_g_hgvs` reads better than the imported-c.HGVS path.

---

## Tests

In `classification/tests/`:

1. **`imported_as_c_hgvs`** — `True` for `imported_c_hgvs` set, `True` for `imported_transcript` set, `False` for g.HGVS only.
2. **g.HGVS-only, resolves to genomic form in both builds** — `_calculate_validation()` returns no `transcript_type_not_supported`, no `liftover` diff, and `include=True`.
3. **g.HGVS-only, resolves to a transcript in both builds with differing c.HGVS** — `liftover` diff present with every severity `"W"`, and `include=True`.
4. **c.HGVS import with an unsupported transcript prefix** — `transcript_type_not_supported` still fires as `"E"` and `include=False`.
5. **Build coverage** — for a g.HGVS-only record, `missing_37`/`missing_38` track variant presence per build.
