# Plan: Fix Intergenic Classifications (#1063)

**Issue**: Intergenic variants (e.g. `NC_000009.12:g.99087902A>C`) have no transcript and no c.HGVS by design. The validation system incorrectly treats them as having resolution errors, causing false "Resolution Differences" that mark them as un-exportable, generate bad UI links, and flood the classification record with spurious errors.

**Reported symptoms**:
- `General Transcript Type Not Supported` ERROR
- `Liftover c.nomen Change` ERROR
- `Liftover Transcript Version Change` ERROR
- `Builds Missing 37` / `Builds Missing 38` WARNINGs
- HGVS resolution tool link shows `?hgvs=None` (causes error when visited)
- No link to allele on classification page (despite allele existing)
- Variant excluded from exports (`include=False`)

---

## Root Cause Analysis

### What identifies a truly intergenic variant?
An `ImportedAlleleInfo` is intergenic **by design** when:
- `imported_g_hgvs` is set (e.g. `NC_000009.12:g.99087902A>C`)
- `imported_c_hgvs` is None and `imported_transcript` is None
- **AND** the HGVSMatcher found no overlapping transcript when resolving the variant (i.e. `ResolvedVariantInfo.transcript_version` is null for all resolved builds)

If the resolver DID find a transcript (i.e. `transcript_version` is set on a `ResolvedVariantInfo`), the missing c.HGVS is **accidental** — it's a resolution failure, not expected behaviour, and should still be flagged as an error.

### Why each error fires

**1. `transcript_type_not_supported` ERROR**
`classification_variant_info_models.py:587`
```python
if not ImportedAlleleInfo.is_supported_transcript(self.get_transcript):
    general["transcript_type_not_supported"] = "E"
```
`get_transcript` returns `None` for intergenic (no `imported_transcript`, no `imported_c_hgvs` to extract from).
`is_supported_transcript(None)` returns `False` at line 712 (`if not transcript_or_hgvs: return False`).
**This check is meant for unsupported transcript prefixes (e.g. `NX_`), not "no transcript at all".**

**2. `Liftover c.nomen Change` + `Liftover Transcript Version Change` ERRORs / WARNINGs**
`classification_variant_info_models.py:573-575`
`ResolvedVariantInfo.recalc_c_hgvs()` (lines 192–217) calls:
```python
result = hgvs_matcher.variant_to_hgvs_variant_used_converter_type_and_method(variant, imported_transcript=None)
```
When `imported_transcript` is None for intergenic variants, the matcher falls back to the genomic (g.) HGVS representation (e.g. `NC_000009.12:g.99087902A>C`). This **g.HGVS gets stored in the `c_hgvs` field** of `ResolvedVariantInfo`.

Then `_calculate_validation()` computes `c_hgvs_obj` from `ResolvedVariantInfo.c_hgvs` for both builds, and diffs them:
```python
if normalized_c_hgvs and lifted_c_hgvs:
    if lifted_diff_dict := calculate_diff_dict(normalized_c_hgvs.diff(lifted_c_hgvs)):
        validation_dict["liftover"] = lifted_diff_dict
```
Since GRCh37 and GRCh38 have **different contig versions** (`NC_000009.11` vs `NC_000009.12`), the diff always flags changes. These are spurious for intergenic variants.

**3. `Builds Missing 37` / `Builds Missing 38` WARNINGs**
`classification_variant_info_models.py:578-581`
```python
if not self.grch37 or not self.grch37.c_hgvs_obj:
    builds["missing_37"] = "W"
if not self.grch38 or not self.grch38.c_hgvs_obj:
    builds["missing_38"] = "W"
```
For the partially-overlapping deletion case, even though liftover succeeds and a `ResolvedVariantInfo` exists for both builds, `c_hgvs_obj` may be `None` because no transcript representation could be generated. The "builds" check is `c_hgvs_obj`-specific, so it triggers even when the variant IS resolved in both builds.

**4. HGVS resolution tool link shows `?hgvs=None`**
`imported_allele_info_detail.html:60`
```html
<a href="{% url 'hgvs_resolution_tool' %}?genome_build={{ allele_info.imported_genome_build|urlencode }}&hgvs={{ allele_info.imported_c_hgvs|urlencode }}">
```
`imported_c_hgvs` is `None` for intergenic, so Django's template engine renders the literal string `"None"`.
Same bug exists in `management/commands/imported_allele_info_check.py:116`.

**5. No allele link on classification page**
This is a downstream effect: because validation errors mark `include=False`, and possibly because the classification view omits the allele link when the variant isn't resolved to a "supported" state. Needs verification in `classification_view.py` — the allele link context variable may be gated on `latest_validation.include`.

---

## Proposed Changes

### Step 1: Add `has_c_hgvs` property to `ImportedAlleleInfo`

**File**: `classification/models/classification_variant_info_models.py`
Add after the `get_transcript` property (~line 665):

```python
@property
def has_c_hgvs(self) -> bool:
    """True when a c.HGVS (transcript-based) representation is expected for this record.

    If c.HGVS or a transcript was explicitly provided on import, c.HGVS is always expected.

    If only g.HGVS was imported, we ask whether the resolver found a transcript when
    resolving the variant. If it did (transcript_version is set on any ResolvedVariantInfo),
    the variant overlaps a known transcript and c.HGVS should have been produced — its
    absence is a resolution failure (missing by accident), not expected behaviour.
    Only if no transcript was found in any build is the variant confirmed as truly intergenic
    (missing by design)."""
    if self.imported_c_hgvs or self.imported_transcript:
        return True
    # g.HGVS-only import: check whether resolution found a transcript in any build
    for build_info in [self.grch37, self.grch38]:
        if build_info and build_info.transcript_version_id:
            return True
    return False
```

---

### Step 2: Fix `_calculate_validation()` for intergenic variants

**File**: `classification/models/classification_variant_info_models.py`, `_calculate_validation()` (~lines 551-595)

#### 2a. Skip `transcript_type_not_supported` for variants without c.HGVS
Change:
```python
if not ImportedAlleleInfo.is_supported_transcript(self.get_transcript):
    general["transcript_type_not_supported"] = _VALIDATION_TO_SEVERITY.get("transcript_type_not_supported", "E")
```
To:
```python
if self.has_c_hgvs and not ImportedAlleleInfo.is_supported_transcript(self.get_transcript):
    general["transcript_type_not_supported"] = _VALIDATION_TO_SEVERITY.get("transcript_type_not_supported", "E")
```

#### 2b. Skip liftover/normalize diffs when no c.HGVS is expected
The diff comparisons at lines 569-575 use `c_hgvs_obj` from `ResolvedVariantInfo`. For intergenic variants this field may hold a g.HGVS (different contig versions across builds). Wrap the entire diff section:
```python
if self.has_c_hgvs:
    if imported_c_hgvs and normalized_c_hgvs:
        if normal_diff_dict := calculate_diff_dict(imported_c_hgvs.diff(normalized_c_hgvs)):
            validation_dict["normalize"] = normal_diff_dict
    if normalized_c_hgvs and lifted_c_hgvs:
        if lifted_diff_dict := calculate_diff_dict(normalized_c_hgvs.diff(lifted_c_hgvs)):
            validation_dict["liftover"] = lifted_diff_dict
```

#### 2c. Fix `builds` check when no c.HGVS is expected
For intergenic variants the absence of `c_hgvs_obj` is expected; check for variant existence instead:
```python
builds: ImportedAlleleValidationTagsBuilds = {}
if self.has_c_hgvs:
    if not self.grch37 or not self.grch37.c_hgvs_obj:
        builds["missing_37"] = _VALIDATION_TO_SEVERITY.get("missing_37", "E")
    if not self.grch38 or not self.grch38.c_hgvs_obj:
        builds["missing_38"] = _VALIDATION_TO_SEVERITY.get("missing_38", "E")
else:
    # For intergenic, check the variant itself exists in each build (no c.HGVS expected)
    if not self.grch37 or not self.grch37.variant_id:
        builds["missing_37"] = _VALIDATION_TO_SEVERITY.get("missing_37", "E")
    if not self.grch38 or not self.grch38.variant_id:
        builds["missing_38"] = _VALIDATION_TO_SEVERITY.get("missing_38", "E")
```

---

### Step 3: Fix HGVS resolution tool link in template

**File**: `classification/templates/classification/imported_allele_info_detail.html`, line 60

Change:
```html
<a class="hover-link" href="{% url 'hgvs_resolution_tool' %}?genome_build={{ allele_info.imported_genome_build|urlencode }}&hgvs={{ allele_info.imported_c_hgvs|urlencode }}">Click here to test</a>
```
To (use `imported_c_hgvs` if set, else `imported_g_hgvs`):
```html
<a class="hover-link" href="{% url 'hgvs_resolution_tool' %}?genome_build={{ allele_info.imported_genome_build|urlencode }}&hgvs={{ allele_info.imported_c_hgvs|default:allele_info.imported_g_hgvs|urlencode }}">Click here to test</a>
```

Also fix the same bug in the management command:
**File**: `classification/management/commands/imported_allele_info_check.py`, line 116
Change:
```python
"hgvs": self.imported_allele_info.imported_c_hgvs
```
To:
```python
"hgvs": self.imported_allele_info.imported_c_hgvs or self.imported_allele_info.imported_g_hgvs
```

---

### Step 4: Investigate and fix "no allele link" on classification page

**File to investigate**: `classification/views/classification_view.py`
Check whether the context variable that drives the allele link is gated on `latest_validation.include`. If so, it should be gated on `allele_info.allele` being non-null instead (the allele exists even if exports are excluded).

Look for pattern like:
```python
context['allele'] = classification.allele if classification.variant_info.latest_validation.include else None
```
If found, change to always pass the allele when it exists:
```python
context['allele'] = classification.allele  # allele link is independent of export eligibility
```

---

### Step 5: ClinVar export — improve error message for intergenic

**File**: `classification/models/clinvar_export_convertor.py`, ~line 476
The current code already excludes intergenic variants from ClinVar (no c.HGVS → error). This is correct behaviour since ClinVar requires c.HGVS. However, the error message should be clearer:

Change:
```python
return ValidatedJson(None, JsonMessages.error(f"No normalised c.hgvs in genome build {genome_build}"))
```
To:
```python
if allele_info and not allele_info.has_c_hgvs:
    return ValidatedJson(None, JsonMessages.error("Intergenic variant — no c.HGVS available. ClinVar requires a transcript-based HGVS."))
else:
    return ValidatedJson(None, JsonMessages.error(f"No normalised c.hgvs in genome build {genome_build}"))
```

---

### Step 6: Re-validate existing intergenic ImportedAlleleInfo records

After deploying the above code changes, existing records that were incorrectly validated need to be re-processed. This requires a **data migration** or a **management command run**:

```python
# One-off migration or manage.py command:
from classification.models import ImportedAlleleInfo
for iai in ImportedAlleleInfo.objects.filter(imported_g_hgvs__isnull=False, imported_c_hgvs__isnull=True):
    iai.apply_validation(force_update=True)
    iai.save()
```

This will clear the false error tags and set `include=True` for intergenic variants that successfully resolved to a variant coordinate in both builds.

---

## Files to Change (Summary)

| File | Change |
|------|--------|
| `classification/models/classification_variant_info_models.py` | Add `has_c_hgvs` property; fix `_calculate_validation()` (3 sub-changes) |
| `classification/templates/classification/imported_allele_info_detail.html` | Fix resolution tool link to use `imported_g_hgvs` as fallback |
| `classification/management/commands/imported_allele_info_check.py` | Fix `hgvs` param in resolution tool URL |
| `classification/models/clinvar_export_convertor.py` | Better error message for intergenic variants |
| `classification/views/classification_view.py` | Investigate and fix allele link visibility (Step 4) |
| New migration or management command | Re-validate existing intergenic records |

---

## Test Cases to Add

In `classification/tests/models/test_classification_utils.py` or a new test file:

1. **`has_c_hgvs` property**:
   - `True` when `imported_c_hgvs` or `imported_transcript` is set (always expects c.HGVS)
   - `True` when g.HGVS-only import AND resolver found a transcript (`transcript_version_id` set on a `ResolvedVariantInfo`) — missing c.HGVS is an error
   - `False` when g.HGVS-only import AND resolver found no transcript in any build — truly intergenic, no c.HGVS expected
2. **Validation tags for intergenic**: Create an `ImportedAlleleInfo` with `imported_g_hgvs` set, assert `_calculate_validation()` does NOT produce `transcript_type_not_supported`, liftover diffs, or normalize diffs.
3. **Builds missing check**: Verify `missing_37`/`missing_38` only fires when the variant itself is absent in that build for intergenic records.
4. **`include=True` for resolved intergenic**: Verify a fully-resolved intergenic record (variant in both builds) gets `include=True`.
