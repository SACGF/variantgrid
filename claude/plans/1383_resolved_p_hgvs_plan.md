# #1383 — Compute and store canonical p.HGVS on ResolvedVariantInfo

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-10-15
Status: in progress

[#1383](https://github.com/SACGF/variantgrid/issues/1383). Split off from #19 (display), which this unblocks along
with variantgrid_private#3557 (mismatch warnings). The measurements behind the source order are in the issue comments.

## The problem

`ResolvedVariantInfo.resolved_hgvs` is the one resolved c.HGVS every grid, export and sort reads. There is no
equivalent for p.HGVS: the grids read `published_evidence__p_hgvs__value` - whatever the lab typed, with mixed
brackets, 1-letter codes, `p.(?)`, typos and sometimes nothing at all - which is why Shariant sets
`CLASSIFICATION_GRID_SHOW_PHGVS = False` (`variantgrid/settings/env/shariantcommon.py`). We can compute a p.HGVS
ourselves from the resolved c.HGVS with the biocommons mapper we already run (`BioCommonsHGVSConverter.am.c_to_p`),
using cdot's protein accession and the transcript sequences we hold. On a dev database it produced a value for
1,096 of 1,098 transcript-bearing GRCh38 rows in ~6 s; on vg-test2 a spot check of the first rows gave
`NM_000492.4:c.1438G>T` → `p.(Gly480Cys)` and an intronic `c.234+83_234+84del` → `p.?`.

## Data

`classification/models/classification_variant_info_models.py:ResolvedVariantInfo` gains two columns
(migration classification/migrations/0228_resolvedvariantinfo_p_hgvs.py (new file), after `classification/migrations/0227_data_evidencekeys.py`):

```python
class ResolvedVariantInfo(TimeStampedModel):
    # ... existing fields ...
    p_hgvs = TextField(null=True, blank=True)
    """ Canonical p.HGVS for resolved_hgvs: p-dot only (no protein accession), 3-letter, Ter, always in predicted
        brackets - p.(Gly12Asp), p.(Met246IlefsTer14), p.?, p.(=). Null when the row has no transcript_version
        (gene-level, mito, g.HGVS-only) or no source could translate it (see p_hgvs_error) """

    p_hgvs_source = TextField(null=True, blank=True, choices=PHGVSSource.choices)
    """ Which source p_hgvs came from """

    p_hgvs_error = TextField(null=True, blank=True)
    """ Why p_hgvs is null when a transcript_version is set - the computed translation's error text
        (e.g. an Ensembl CDS not divisible by 3) """


class PHGVSSource(TextChoices):
    VEP = "V", "VEP"             # VariantTranscriptAnnotation.hgvs_p, exact transcript version, latest VAV
    COMPUTED = "C", "Computed"   # biocommons c_to_p from resolved_hgvs
    CLINGEN = "G", "ClinGen"     # cached ClinGenAllele proteinEffect
    EVIDENCE = "E", "Evidence"   # a classification's p_hgvs evidence value, parsed
```

No converter-version column: a value legitimately changes when a new VEP version lands, a transcript sequence gets
cached or ClinGen is refreshed, and the source column plus an idempotent, rerunnable backfill is enough. The short form
(`p.(G12D)`, `p.(M246Ifs*14)`) is derived on the fly for display; nothing sorts or searches on it.

The dataclass the resolver returns (one value, one source, or an error):

```python
@dataclass(frozen=True)
class PHGVSResolution:
    p_hgvs: Optional[str]
    source: Optional[PHGVSSource]
    error: Optional[str]
```

## The canonical form and its parser

New module genes/hgvs/phgvs_canonical.py (new file), owning "turn any p.HGVS string, from any source, into the one
canonical string" and the short form. It wraps biocommons `hgvs.parsers.Parser` (the regex
`genes/hgvs/phgvs.py:PHGVS` only handles single-residue changes and cannot parse ranges, dups or delins, which are
2-12% of real values; it stays for its existing callers `classification/models/evidence_mixin.py:EvidenceMixin.p_parts`
and `ClassificationGrouping.p_hgvses`).

Entry points:

- `canonical_p_hgvs(raw: str) -> Optional[str]` - returns the canonical string or None when the string is not a
  p.HGVS biocommons can parse. Steps: strip whitespace; split off anything before `p.` (protein accession, with or
  without a gene in brackets) and parse only the p-dot; map legacy `X` for stop to `Ter` before parsing (biocommons
  reads `X` as `Xaa`, unknown residue, and `p.Trp26X` would come out as `p.(Trp26Xaa)`); url-decode `%3D` (VEP
  writes `p.%3D` for `p.=`); a bare p-dot has no accession so it goes through `Parser.parse_p_posedit`;
  `p.?` has `posedit is None` and is returned as `p.?` directly. Format with
  `format(conf={"p_3_letter": True, "p_term_asterisk": False})`, then force predicted brackets - `p.(…)` - on anything
  that is not `p.?` or `p.=`/`p.(=)` (write `p.(=)`), matching `p_parts`' `override_is_confirmed_to=False`.
- `short_p_hgvs(canonical: str) -> str` - the same parse formatted with `p_3_letter=False, p_term_asterisk=True`.
- `canonical_from_sequence_variant(variant: hgvs.sequencevariant.SequenceVariant) -> str` - what the computed path
  calls; same formatting rules, no reparse.
- `is_frameshift_in_disguise(variant) -> bool` - a computed `delins` whose range ends at the protein's last residue
  is a frameshift biocommons could not name (it writes `p.(Ser1321_Val2843delinsIleLeuTer)` where VEP writes
  `p.Ser1321IlefsTer3`); the resolver treats it as "no computed value".

## Resolving

`ResolvedVariantInfo.recalc_p_hgvs() -> PHGVSResolution`, next to `recalc_resolved_hgvs()`, tries the sources in
order and returns the first canonical value:

1. **VEP** - `VariantTranscriptAnnotation` for `self.variant`, `version=AnnotationVersion.latest(genome_build).variant_annotation_version`
   and `transcript_version=self.transcript_version` (exact version only: another version of the transcript can have a
   different protein). Use `AnnotationVersion.latest_or_none` so a build with no annotation version falls through
   rather than raises. `hgvs_p` null/empty falls through.
2. **Computed** - `HGVSMatcher.instance(genome_build).hgvs_converter` is the matcher's converter and may be the
   ClinGen one by configuration, so the computed path constructs its own `BioCommonsHGVSConverter(genome_build)`
   (cached per build the way `HGVSMatcher.instance` is) and calls
   `converter.am.c_to_p(Parser().parse(converter.c_hgvs_remove_gene_symbol(self.resolved_hgvs)))`
   (`resolved_hgvs` carries the gene in brackets, `NM_x.y(GENE):c.…`). Catch `Exception`, not just `HGVSError`: the
   Ensembl CDS-length case is a `NotImplementedError`
   (`Transcript ENST00000361381.2 is not supported because its sequence length of 1378 is not divisible by 3`). The
   error text is kept for `p_hgvs_error` if nothing later succeeds. A result that `is_frameshift_in_disguise` falls
   through.
3. **ClinGen** - `self.allele_info.allele.clingen_allele.get_p_hgvs(self.transcript_version.accession, match_version=True)`
   from the cached JSON only (the model method reads `api_response`; this path never calls the registry). Reject
   anything whose p-dot does not start with `p.` (ClinGen returns `NP_x:n.123+1G>A` for intronic variants, the same
   bug class `ClinGenAllele.get_c_hgvs_variant` already works around; the right answer there is `p.?`, which the
   computed path already gave if it ran). A value ending `fsTer?` is accepted only when no computed value exists.
4. **Evidence** - `published_evidence__p_hgvs__value` of any classification with `allele_info=self.allele_info`
   (`Classification.allele_info`), latest published modification first; the first that `canonical_p_hgvs` can parse
   wins. One query, `values_list`, no loop over Classification objects.

Nothing translates → `PHGVSResolution(None, None, error)` where `error` is the computed path's message, or
"No p.HGVS source" when the computed path did not run.

### Hook

In `set_variant_and_save()`: reset `p_hgvs`, `p_hgvs_source` and `p_hgvs_error` with the other derived fields at
the top; after `recalc_resolved_hgvs()` succeeds and `self.transcript_version` is set, call `recalc_p_hgvs()` and
copy the three values. The gene-level early return and the mito / g.HGVS-only cases (no `transcript_version`) leave
all three null. A failure inside `recalc_p_hgvs` is recorded in `p_hgvs_error`, never raised - a p.HGVS problem
must not stop the c.HGVS from saving. The existing refresh paths
(`ImportedAlleleInfo._update_variant(force_update=True)`, `ImportedAlleleInfo.relink_variants`) then pick it up.

### Autopopulate agreement

`classification/autopopulate_evidence_keys/evidence_from_variant.py:get_evidence_fields_from_preferred_transcript`
already does "VEP, else ClinGen" for the `p_hgvs` evidence key (line ~437). It keeps its own flow (the evidence value
keeps the protein accession and the `=`→`?` splice rule under `CLASSIFICATION_AUTO_POPULATE_P_HGVS_SYNONYMOUS_SPLICE_CHANGE_TO_UNKNOWN`,
because that is what the lab is asked to confirm), and is out of scope here.

## Backfill

New command classification/management/commands/resolved_variant_info_p_hgvs.py (new file), `category = "maintenance"`:

- Pages `ResolvedVariantInfo.objects.filter(transcript_version__isnull=False)` with
  `library/django_utils/__init__.py:chunked_queryset` (pk range, default chunk 1,000, `--chunk-size`), running
  `recalc_p_hgvs()` per row and writing with `save(update_fields=["p_hgvs", "p_hgvs_source", "p_hgvs_error"])`
  - the two c.HGVS columns are left alone, so this is cheap and safe to rerun. `--genome-build` and `--only-missing`
  (rows with `p_hgvs IS NULL`) narrow it; `--retrieve-transcripts` builds the converter with
  `DjangoTranscriptDataProvider(genome_build, retrieve_transcripts=True)` for a second, slower pass that fetches
  missing transcript sequences from NCBI (default off, so the first pass cannot fan out into thousands of fetches).
- Prints a summary: rows per source, rows with an error and the top error messages by count (aggregated in SQL
  from `p_hgvs_error` after the run, not accumulated in Python).
- Registered as a `ManualOperation` in migration 0228
  (`manual/operations/manual_operations.py`, example `snpdb/migrations/0188_one_off_migrate_common_filter_gnomad_versions.py`)
  with `test=` "any ResolvedVariantInfo with a transcript_version and p_hgvs null and p_hgvs_error null".

Rows are in the tens of thousands at most (vg-test2: 50 k classifications; a dev box's 1,098 rows took ~6 s), so a
single process is fine; the chunking is for memory and for resuming, not throughput.

## Display (#19 reads it; this issue only makes it available)

Add `p_hgvs_short` as a property on `ResolvedVariantInfo` (`short_p_hgvs(self.p_hgvs)`), and show both on the
admin / `vg inspect classification` output so the backfill can be eyeballed. The grid and template-tag consumers
(`classification/views/classification_datatables.py` line ~104, `classification/views/classification_grouping_datatables.py`
lines ~178 and ~487, `classification/templatetags/classification_tags.py` line ~638) switch to the new column in #19.

## Tests

`classification/tests/utils/test_phgvs.py` and `classification/tests/models/test_imported_allele_info.py` are the
neighbours; `snpdb/tests/utils/mock_clingen_api.py` has `proteinEffect` payloads. Fixtures per
`claude/guides/testing.md`. Worth keeping:

- `canonical_p_hgvs`: 1-letter → 3-letter, `*` → `Ter`, `X` → `Ter`, missing brackets added, accession stripped,
  `p.%3D` → `p.(=)`, `p.?` unchanged, a range/dup/delins (`p.(Lys550_Lys558del)`) round-trips, junk returns None.
- `short_p_hgvs`: `p.(Met246IlefsTer14)` → `p.(M246Ifs*14)`.
- `recalc_p_hgvs` source ordering: VEP value wins over a computed one; computed wins when VEP has no row for the
  exact transcript version; a ClinGen `NP_x:n.…` is rejected; an evidence typo is skipped and the next classification's
  value used; the dup-frameshift guard falls through.
- `set_variant_and_save`: a row without a transcript_version stays null in all three columns; a computed failure
  lands in `p_hgvs_error` and the c.HGVS still saves.

Tests that only restate a field declaration or check that biocommons formats what biocommons formats are dropped at the
audit.
