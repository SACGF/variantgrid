# Multi variant report: tick the assay flags from the measures

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-22; phase 1 implemented by Claude Opus 5 (claude-opus-5), revised by Claude Fable 5.1 the same day
Status: landed - phase 1 in core ee184f8ba (PR 1895) and sapath PR 456, with the lab's thresholds and derived Confirmed as a follow-up straight on master; phase 2 is [sapath#455](https://github.com/SACGF/variantgrid_sapath/issues/455), planned in `claude/plans/sapath455_library_qc_metrics_output_plan.md`

[sapath#454](https://github.com/SACGF/variantgrid_sapath/issues/454). The TSO 500 Build report form asks the scientist
to tick seven "Assay Success" boxes and three "Caveats" by hand. Most of those answers are already in the database, or
in a DRAGEN file we are not yet importing. This plan makes the form show each measure against its checkbox and start
the tick from a rule, leaving the scientist to adjust.

## Where the checkboxes come from

- *sapath/tso500_case_template.py* declares them in `CASE_FIELDS`; the sapath migration `0014_tso500_report_template`
  loads that list into the "TSO 500" row of `classification/models/classification_report_models.py:ClassificationReportTemplate`
  (`case_fields`). A change to the list lands as a new sapath migration that updates the row.
- They were literals in the legacy PDF parser (`Assay Success` all True, `Caveats` all False). Nobody ever chose them.
  They only reach the Omico JSON (*sapath/tso500_report.py* `_booleans`); the printed document ignores them.
- `analysis/views/views_classify_report.py:_case_values_for_form` already prefills a text field from an evidence key
  (`prefill_key`). Nothing prefills from a measure.
- `classification/report/case_report_context.py:build_report_context` already gives the report the specimen's
  `patients/models.py:SpecimenMeasure` rows (`measures.tmb/msi/gis/tumour_fraction/ploidy`), each with value, unit,
  call, threshold and the DRAGEN section it was read from in `source_payload` (so `Usable MSI Sites` is there).
- The CombinedVariantOutput import (`upload/tso500/dragen_combined_variant_output_records.py:write_specimen_measures`)
  writes the value and never a `call`. So today every DRAGEN case sends Omico `MSI: Not able to be determined`.
  Mocha's tumour fraction does carry the pathologist's call (`Insufficient`, `No tumour`).

What each flag actually means, and what can answer it:

| Flag | Meaning | Source today | Source after phase 2 |
|---|---|---|---|
| TMB | a TMB was reportable | `tmb` measure has a value and a call | small-variant/TMB library QC |
| MSI | an MSI status was reportable | `msi` usable sites reach the lab minimum | MSI library QC |
| Amplifications | CNV calling passed QC | none | CNV library QC |
| Variants | small variant calling passed QC | none | small-variant/TMB library QC |
| Fusions | RNA arm passed QC | none | RNA library QC |
| Deletions | (assay is not validated for loss) | none | none - ask the lab |
| Confirmed | orthogonal confirmation? | none | none - ask the lab |
| Purity caveat | tumour content too low | `tumour_fraction` call or value | same |
| Quality caveat | library QC marginal | none | any library QC failed |
| Fail caveat | assay failed | none | DNA library QC / contamination failed |

## Data

### Phase 1: no new model. Two case_field keys and the calls the import stops leaving blank

A `case_fields` entry (JSON on the template row, editable in admin) gains two optional keys:

```python
# ClassificationReportTemplate.case_fields - one bool field
{"key": "assay_success_msi", "label": "MSI", "type": "bool", "default": True, "group": "assay_success",
 "measure": "msi",                       # context key (MEASURE_CONTEXT_KEYS) shown beside the checkbox
 "tick_when": {"called": True}}          # rule the tick starts from; absent = the field's default

# tick_when, one of:
{"called": True}                         # the measure has a call
{"call_in": ["Insufficient", "No tumour"]}
{"value_below": 20}                      # in the measure's own unit
```

`SpecimenMeasure` is unchanged. The import starts filling the three columns it already has:

| measure | `call` | `threshold` | `threshold_source` |
|---|---|---|---|
| MSI | `Stable` / `Unstable`, or `None` when usable sites are below the minimum | `>= {min_sites} usable sites, >= {pct}% unstable` | `settings.TSO500_MSI_*` |
| TMB | `High` / `Low` | `>= {mut_per_mb} mut/Mb` | `settings.TSO500_TMB_HIGH_MUT_PER_MB` |

Three settings, default `None` (no call written - an installation that has not set its policy sends
`Not able to be determined`, as now): `TSO500_MSI_MIN_USABLE_SITES`, `TSO500_MSI_UNSTABLE_PERCENT`,
`TSO500_TMB_HIGH_MUT_PER_MB`. SA Path's values go in its env settings - Illumina's guideline is 20 sites and 20% for MSI,
10 mut/Mb for TMB, but they are the lab's numbers to confirm.

### Phase 2 ([sapath#455](https://github.com/SACGF/variantgrid_sapath/issues/455)): library QC, so Amplifications / Variants / Fusions / Quality / Fail have a source

DRAGEN writes *MetricsOutput.tsv* beside the CombinedVariantOutput. The fixture is in - `upload/test_data/tso500/ExampleSample_2600000001/ExampleSample_2600000001_MetricsOutput.tsv`
(quirks in `upload/test_data/tso500/README.md`) - and it shows exactly what there is to store:

| Section | Shape | Use |
|---|---|---|
| `[Header]` | key/value: Output Date / Time, Workflow Version | `measured_date`, `method` - as `measured_date()` reads the CVO's `[Analysis Details]` |
| `[Run QC Metrics]` | metric, LSL, USL, one `Value` column | the run's, not a library's - out of scope here (NA when analysis started from FASTQ) |
| `[Analysis Status]` | one column per sample: `COMPLETED_ALL_STEPS`, `FAILED_STEPS`, `STEPS_NOT_EXECUTED` | a sample whose steps failed is the **Fail** caveat |
| `[DNA Library QC Metrics]` | `CONTAMINATION_SCORE` USL 1457 | category `DNA` |
| `[... for Small Variant Calling and TMB]` | `MEDIAN_INSERT_SIZE` >=70, `MEDIAN_EXON_COVERAGE` >=150, `PCT_EXON_50X` >=90 | category `SMALL_VARIANT_TMB` - the **Variants** and **TMB** flags |
| `[... for MSI]` | `USABLE_MSI_SITES` >=40 | category `MSI` - the same number `TSO500_MSI_MIN_USABLE_SITES` gates the call on |
| `[... for CNV]` | `GENE_SCALED_MAD` <=0.134, `MEDIAN_BIN_COUNT_CNV_TARGET` >=1 | category `CNV` - the **Amplifications** flag |
| `[... for GIS]` | `PCT_TARGET_HRD_50X` >=50 | category `GIS` - not a flag today, but the HRD block's reliability |
| `[RNA Library QC Metrics]` | `MEDIAN_CV_GENE_500X` <=0.93, `TOTAL_ON_TARGET_READS` >=9M, `MEDIAN_INSERT_SIZE` >=80 | category `RNA` - the **Fusions** flag and the **Quality** caveat |
| `[DNA Expanded Metrics]`, `[RNA Expanded Metrics]` | no guidelines | kept in the row's `metrics` for the record, judged by nothing |

The guideline columns are the file's own LSL/USL and match the lab's methods paragraph (40 sites, 1457, 0.134, 150x,
70bp, 90%; RNA 0.93, 9M, 80bp), so `passed` is "every metric in the section within the file's guideline": `NA` as a
guideline is no bound, `NA` as a value is an arm the section does not apply to (the RNA column of a DNA section) and
counts as neither. Every library QC section carries a column per sample, so one file writes rows for both arms.

```python
class LibraryQCCategory(models.TextChoices):
    DNA = 'D', 'DNA library'                  # contamination
    SMALL_VARIANT_TMB = 'V', 'Small variants and TMB'
    MSI = 'M', 'MSI'
    CNV = 'C', 'CNV'
    GIS = 'G', 'GIS'
    RNA = 'R', 'RNA library'                  # fusions and splice variants


class LibraryQC(GuardianPermissionsMixin, TimeStampedModel):
    """ One DRAGEN QC category for one sequenced library - what 'the assay succeeded for X' means """
    extraction = models.ForeignKey(Extraction, on_delete=CASCADE)
    category = models.CharField(max_length=1, choices=LibraryQCCategory.choices)
    passed = models.BooleanField(null=True)    # None: no metric in the section applied to this arm
    completed = models.BooleanField(null=True)  # [Analysis Status] COMPLETED_ALL_STEPS for the arm, same on every row
    metrics = models.JSONField(default=dict)   # {metric: {"value", "unit", "lsl", "usl"}} - the section, so 'why' stays answerable
    method = models.TextField(blank=True)      # Workflow Version off [Header]
    measured_date = models.DateTimeField(null=True, blank=True)
    user = models.ForeignKey(User, null=True, on_delete=SET_NULL)

    class Meta:
        unique_together = ("extraction", "category")   # a re-analysis replaces, as SpecimenMeasure does
```

Implementation, when the lab wants the five flags derived:

1. **Parser.** The file is the CVO's layout (banner line, `[Section]`, tab padding), so the section reader in
   `upload/tso500/dragen_combined_variant_output_parser.py` is shared and a `dragen_metrics_output_parser.py` reads the
   sections above by name; `can_process_file` off the banner "Metrics Output". A new `UploadedFileTypes` value and task
   beside `upload/tasks/import_dragen_tso500_combined_variant_output_task.py`.
2. **Records.** The sample columns are the CVO's DNA / RNA Sample IDs, so each resolves to its extraction the way
   `upload/tso500/dragen_combined_variant_output_records.py:resolve_pair` does (the sample ID's accession and container
   suffix); a column naming no extraction is the same needs-attention path as a CVO arm that names none. One
   `LibraryQC` per (extraction, category), upserted.
3. **Form.** A case_field says `"qc": "CNV"` (a `LibraryQCCategory` name) and ticks when that category passed on any
   of the case's extractions; `"qc_failed": ...` for the Quality and Fail caveats. The note under the group lists
   the category's metrics against their guideline the way phase 1 lists the measure's threshold, so a scientist
   sees which number failed.
4. **sapath.** Amplifications → CNV, Variants → SMALL_VARIANT_TMB, Fusions → RNA, Quality → RNA failed,
   Fail → DNA failed or not `completed`; a migration updating the template row.
5. **Pipeline.** *vg_api_full.py* (../NGS-pipelines) uploads the file with the run, as #443 wired the CVO.

## Phase 1 implementation

1. **Calls at import.** `upload/tso500/dragen_combined_variant_output_records.py`: `MeasureSource` gains an optional
   `call` callable taking the section values and returning `(call, threshold, threshold_source)` or `None`; MSI and TMB
   get one, reading the settings above. `write_specimen_measures` passes the three columns through. Settings and their
   comment in `variantgrid/settings/components/default_settings.py` next to `TSO500_PAIR_ID_PATIENT_CODE_REGEX`.
2. **Measures reach the form.** `classification/report/case_report_context.py`: extract the specimen-for-source lines of
   `build_report_context` into `case_specimen(source_level, source)` and reuse it in
   `analysis/views/views_classify_report.py:case_report_build_dialog`, which puts `measures=_specimen_measures(specimen)`
   into the template context and passes them to `_case_values_for_form`.
3. **Rules.** `_case_values_for_form` gains `measures`; a bool field with `tick_when` starts from the rule against
   `measures[field["measure"]]` when the measure exists, else from `default`. A draft's own answers still win. The rule
   evaluation is one small function in `classification/models/classification_report_models.py` beside
   `case_values_from_form`, so a template's JSON is validated in `clean()` (unknown rule or measure key is an error).
4. **Form.** `analysis/templates/analysis/case_report_build_dialog.html`: a bool field with a `measure` renders the
   measure after its label, `<span class="text-muted small">2.48% (Stable)</span>`, using the same
   value/unit/call wording as `patients/grids.py` renders on the specimen grid - pull that into a `SpecimenMeasure`
   `__str__`-style property both use. A field whose measure is missing says `no measure`.
5. **sapath.** *sapath/tso500_case_template.py* sets on the fields:
   - `assay_success_tmb`: `measure: tmb`, `tick_when: {called}`
   - `assay_success_msi`: `measure: msi`, `tick_when: {called}`
   - `caveat_purity`: `measure: tumour_fraction`, `tick_when: {call_in: [Insufficient, No tumour]}` - or `value_below`
     once the lab names the percentage
   - the rest unchanged (Amplifications, Variants, Fusions wait for phase 2; Deletions and Confirmed are the lab's call).
   Migration `0016` updates the existing row's `case_fields` from `CASE_FIELDS`. SA Path env settings gain the three
   thresholds.
6. **Tests.** Core: `analysis/tests/test_case_report.py` - the three rules, a missing measure falling back to the
   default, a draft answer winning over the rule; `upload/tests/test_import_dragen_tso500_combined_variant_output.py` -
   MSI call from usable sites / percent, no call with settings unset. sapath: *tests/test_tso500_report_json.py* is
   unaffected (the JSON reads the answers, not the rules).
7. **Docs.** `classification/AGENTS.md` case report section: the two keys and the rule vocabulary; `scripts/vg docs check`.

## What the lab's own reports say (found after phase 1 landed)

The legacy PDF parser in *../tso500_reports* and its `documents/reporting_method.md` answered most of the questions
phase 1 left open, and the follow-up on master applied them:

- **Thresholds** (the report's methods paragraph): MSI High >= 30% sites unstable, MSI Low >= 10% and < 30%,
  MS Stable < 10%, off > 40 usable MSI sites; TMB high >= 10 mutations/Mb; minimum tumour content 20%. The
  vocabulary is `MSS` / `MSI-Low` / `MSI-High`. The settings are now band lists (`TSO500_MSI_CALL_BANDS`,
  `TSO500_TMB_CALL_BANDS`, `TSO500_MSI_MIN_USABLE_SITES`) so the words and the cut points are both the lab's.
- **Confirmed** is every other Assay Success flag being true (`parse_assay_success`), so the JSON derives it and
  the form no longer asks.
- **Deletions** was never false in any issued report; it stays a plain checkbox.
- **Purity caveat** ticks off the pathologist's call or a tumour content under 20% (`tick_when` takes a list).
- **The policy is on the form**: under each group the build form lists the measure, the rule in words and the
  threshold the import applied with the setting it came from, so a scientist who disagrees asks for the setting
  to change rather than silently unticking.

For phase 2 (#455), the same paragraph gives the QC minimums the flags mean: DNA median exon coverage >150x, median
insert size >70bp, % exon 50x >90%, >40 usable MSI sites, gene scaled MAD <=0.134, contamination score <=1457; RNA
median CV <93% for genes with median coverage >500x, >9M on target reads, median insert size >80bp. And the legacy's
derivations: Fusions false on "Gene Fusions Fail"; Quality caveat on RNA QC not met or reduced sensitivity; Fail
caveat on "sequencing failed for this specimen", which also sets every Assay Success flag false.
