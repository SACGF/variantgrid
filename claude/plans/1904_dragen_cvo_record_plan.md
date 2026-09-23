# Store the DRAGEN TSO500 CombinedVariantOutput as a seqauto record

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-23
Status: draft

[#1904](https://github.com/SACGF/variantgrid/issues/1904). TMB, MSI, GIS, tumour fraction and ploidy are results of one
DRAGEN analysis of one sequenced pair. Today `upload/tso500/dragen_combined_variant_output_records.py:write_specimen_measures`
writes them as `patients/models.py:SpecimenMeasure` rows keyed on (specimen, measure type), replaced on resend, so which
analysis produced a report's TMB is not recorded, two analyses of one specimen cannot coexist, and Mocha's pathology
tumour percent (`sapath/models/sapath_mocha_vg_sync.py` in variantgrid_sapath) and DRAGEN's `Tumor Fraction` overwrite
each other under `SpecimenMeasureType.TUMOUR_FRACTION`. Nothing stores the file itself.

This plan gives the CVO a record of its own in seqauto, shaped like `seqauto/models/models_seqauto.py:LibraryQC` (its
MetricsOutput sibling, [sapath#455](https://github.com/SACGF/variantgrid_sapath/issues/455)): one row per (run, pair),
typed columns for the file's scalars, the standard nullable claims on the seqauto and snpdb objects, and the case report
reading the row for its own samples. `SpecimenMeasure` keeps only what is genuinely of the specimen.

It assumes #1903 step 3 has landed: the CVO import writes no VCF, so it becomes a single-shot import like MetricsOutput.
Until then the row is written from the existing post-header task and the VCF pipeline stays.

## Data

### `DragenTSO500CombinedVariantOutput` - one analysis of one pair

`seqauto/models/models_seqauto.py`, beside `LibraryQC`. Vendor-named, as `IlluminaFlowcellQC` and `FastQC` are; a
second caller's pair summary would be its own record, not rows in this one. A plain `TimeStampedModel`, not a
`SeqAutoRecord`: it arrives by upload, not the disk scanner, and has no path.

```python
class DragenTSO500CombinedVariantOutput(TimeStampedModel):
    # Key - the file names its run 'NA', so the run comes in as upload metadata, as MetricsOutput's does
    sequencing_run_name = models.TextField()
    pair_id = models.TextField()                  # [Analysis Details] Pair ID
    # [Analysis Details]
    dna_sample_id = models.TextField(blank=True)  # a Sample.vcf_sample_name and a SequencingSample.sample_name
    rna_sample_id = models.TextField(blank=True)
    output_datetime = models.DateTimeField(null=True, blank=True)   # Output Date + Output Time, made aware
    module_version = models.TextField(blank=True)
    pipeline_version = models.TextField(blank=True)
    # [TMB]
    total_tmb = models.FloatField(null=True, blank=True)            # mut/Mb
    coding_region_size_mb = models.FloatField(null=True, blank=True)
    passing_eligible_variants = models.IntegerField(null=True, blank=True)
    # [MSI]
    usable_msi_sites = models.IntegerField(null=True, blank=True)
    total_msi_sites_unstable = models.IntegerField(null=True, blank=True)
    percent_unstable_msi_sites = models.FloatField(null=True, blank=True)
    # [GIS]
    genomic_instability_score = models.FloatField(null=True, blank=True)
    tumor_fraction = models.FloatField(null=True, blank=True)       # a fraction, 0.62 - Mocha's is a percent
    ploidy = models.FloatField(null=True, blank=True)
    # Links, all nullable claims filled where the object exists and reconciled after (as LibraryQC)
    sequencing_run = models.ForeignKey(SequencingRun, null=True, blank=True, on_delete=SET_NULL)
    dna_sequencing_sample = models.ForeignKey(SequencingSample, null=True, blank=True, on_delete=SET_NULL,
                                              related_name="cvo_dna_set")
    rna_sequencing_sample = models.ForeignKey(SequencingSample, null=True, blank=True, on_delete=SET_NULL,
                                              related_name="cvo_rna_set")
    dna_sample = models.ForeignKey(Sample, null=True, blank=True, on_delete=SET_NULL, related_name="cvo_dna_set")
    rna_sample = models.ForeignKey(Sample, null=True, blank=True, on_delete=SET_NULL, related_name="cvo_rna_set")
    specimen_reference = models.TextField()       # the accession inside the sample IDs
    specimen = models.ForeignKey(Specimen, null=True, blank=True, on_delete=SET_NULL)
    specimen_match_status = models.CharField(max_length=1, choices=MatchStatus.choices, null=True, blank=True)
    specimen_match_error = models.TextField(null=True, blank=True)
    specimen_match_date = models.DateTimeField(null=True, blank=True)
    # Provenance
    file_upload = models.ForeignKey("upload.FileUpload", null=True, on_delete=SET_NULL)
    user = models.ForeignKey(User, null=True, on_delete=SET_NULL)

    class Meta:
        # A re-analysis of the same run replaces its row; the same pair on a later run is a new one
        unique_together = ("sequencing_run_name", "pair_id")
```

The specimen claim fields and `apply_specimen_match` / `park_specimen_claim` are shared with `LibraryQC` through an
abstract `SpecimenClaimMixin` in seqauto, lifted out of `LibraryQC` unchanged. The four sample-side FKs are what makes
the row findable from a case: a `Sample` by `vcf_sample_name`, a `SequencingSample` by the sheet's `Pair_ID` /
`Sample_Type` data through `seqauto/models/models_seqauto.py:sequencing_sample_for_pair`.

The lab's calls are computed, not stored. `msi_call` and `tmb_call` in
`upload/tso500/dragen_combined_variant_output_records.py` take the row instead of a section dict and move to the model
as properties, so the `TSO500_MSI_MIN_USABLE_SITES`, `TSO500_MSI_CALL_BANDS` and `TSO500_TMB_CALL_BANDS` settings are
applied when the number is read. The threshold and its source were stored on `SpecimenMeasure` so a re-interpretation
survived the row being overwritten; this row is never overwritten by a later analysis, and
`classification/models/classification_report_models.py:CaseReport` snapshots what it printed.

### `SpecimenMeasure` - what is of the specimen

`patients/models.py:SpecimenMeasure` keeps its fields. `patients/models_enums.py:SpecimenMeasureType` loses TMB, MSI, GIS
and PLOIDY; `TUMOUR_FRACTION` ('F') stays with the label "Tumour content (pathology)" so Mocha's rows, its sync and the
sapath template's `tumour_fraction` measure key are untouched. The DRAGEN rows are deleted by the migration below.

### `Measure` - what the report and the build form see

`classification/report/case_report_context.py` gains a frozen dataclass so the case-field `measure` / `tick_when`
machinery, the templates and sapath's `_measure(context, key)` keep reading `context["measures"][key]` with the same
fields whatever the number came from:

```python
@dataclass(frozen=True)
class Measure:
    value: Optional[float]
    unit: Optional[str]
    call: Optional[str]
    threshold: Optional[str]      # the policy in words, describe_bands(...)
    method: str
```

Context keys (`patients/models_enums.py:MEASURE_CONTEXT_KEYS` moves to seqauto's enums as `CVO_MEASURE_CONTEXT_KEYS`
plus the one pathology key):

| key | from |
|---|---|
| `tmb` | row `total_tmb`, `tmb_call` |
| `msi` | row `percent_unstable_msi_sites`, `msi_call` |
| `gis` | row `genomic_instability_score` |
| `ploidy` | row `ploidy` |
| `tumour_fraction_sequencing` | row `tumor_fraction`, unit `fraction` |
| `tumour_fraction` | `SpecimenMeasure` TUMOUR_FRACTION, unit `%` - the pathology estimate, as sapath's purity caveat and `_purity` are written for |

`measure_tick` and `_rule_holds` take a `Measure`; `validate_case_fields` validates a `measure` key against the union.

### `UploadedDragenTSO500CombinedVariantOutput`

`upload/models/models_uploaded_files.py`, mirroring `UploadedDragenTSO500MetricsOutput`:

```python
class UploadedDragenTSO500CombinedVariantOutput(UploadData):
    file_upload = models.OneToOneField(FileUpload, on_delete=CASCADE)
    combined_variant_output = models.OneToOneField(DragenTSO500CombinedVariantOutput, null=True, on_delete=SET_NULL)
```

## Import

1. **Factory.** `upload/import_task_factories/import_task_factories.py:DragenTSO500CombinedVariantOutputImportTaskFactory`
   becomes an `ImportTaskFactory` like `DragenTSO500MetricsOutputImportTaskFactory`: metadata keys
   `{SEQUENCING_RUN}` (required, part of the key), data class the new `UploadData`, one celery task. The pipeline
   already sends `sequencing_run` with the MetricsOutput, so the client and NGS-pipelines changes are the same one
   line for this file type (variantgrid_api and the NGS-pipelines plan alongside; the CVO's `[Sequencing Run Details]`
   is `NA` in the file).
2. **Task.** `upload/tasks/import_dragen_tso500_combined_variant_output_task.py` parses the header and the three
   sections into the row with `update_or_create` on the key, then, as `write_library_qc` does: resolves the specimen
   from the accession in the sample IDs (`SAMPLE_ID_ACCESSION_PATTERN`) or parks the claim, fills `sequencing_run`,
   the two arm `SequencingSample`s through `sequencing_sample_for_pair`, and the two `Sample`s by
   `vcf_sample_name` the way `link_samples_to_extractions` finds them. `resolve_pair` / `link_samples_to_extractions`
   / `link_to_sequencing_run` keep doing the patient chain and the `VCFFromSequencingRun` /
   `SampleFromSequencingSample` links; `write_specimen_measures`, `MeasureSource`, `MEASURE_SOURCES` and the
   `measured_date` helper go.
3. **Reconciliation.** `patients/tasks/extraction_matching_tasks.py:reconcile_pending_extractions` re-resolves the
   row's specimen claim beside `LibraryQC`'s, and `link_library_qc_to_sequencing_samples` generalises to both models
   (both arms for this one). A new pass links `dna_sample` / `rna_sample` for rows whose sample arrived after the file:
   `Sample.vcf_sample_name` equal to the row's sample ID, restricted to VCFs linked to the row's run where it is known.
4. **Signals.** `link_samples_and_vcfs_to_sequencing` in `upload/vcf/vcf_import.py` already runs when an arm VCF
   lands; it gains the one query that fills a waiting row's sample FK.

## Report

`classification/report/case_report_context.py`:

- `case_combined_variant_output(samples, specimen)` picks the row whose `dna_sample` or `rna_sample` is one of the
  case's samples; where none is linked yet, the newest `output_datetime` row claiming the specimen (the
  `specimen_library_qc` rule). One row, so the report's TMB, MSI and GIS come off the same analysis.
- `specimen_measures` becomes `case_measures(cvo, specimen) -> dict[str, Measure]`, the table above.
- `ReportContext` gains `combined_variant_output` beside `measures`; `context_as_dict` snapshots the row's key,
  versions and `output_datetime` under `"analysis"` so a report can say which analysis it printed.
- `LIBRARY_QC_CONTEXT_KEYS` and `specimen_library_qc` are unchanged; a follow-up can pick `LibraryQC` by the same
  sample rule.

`analysis/views/views_classify_report.py` and the build form template read `Measure` fields, which are the names they
already use.

## Pages

- `patients/views.py:view_specimen` lists the specimen's CVO rows under the library QC table (run, pair, versions,
  output date, the five numbers with their calls); the measures table shrinks to the pathology row.
- `patients/grids.py:SpecimenColumns` `measures` column reads only `SpecimenMeasure`, which is now one row at most;
  the "# Measures" column goes.
- `seqauto/views.py` sequencing run page lists the run's CVO rows beside its `LibraryQC`; `seqauto/admin.py`
  registers the model.
- The specimen-measure REST endpoints go: `SpecimenMeasureViewSet`, `SpecimenMeasureBulkCreateView`, their
  serializers, `patients/urls.py` routes and the `specimen_measures` entry in `variantgrid/views_rest.py:API_FEATURES`. #1559's API was
  superseded by the CVO import and has no client.

## Migrations

- seqauto 0049 (new): the model, the mixin refactor of `LibraryQC` (no
  column change).
- upload 0048 (new): the `UploadData` record.
- patients 0019 (new): `RunPython` deleting `SpecimenMeasure` rows whose
  `measure_type` is not TUMOUR_FRACTION or whose `method` starts with `DRAGEN` (the CVO import's `method`), then the
  choices change. Existing rows are recreated by re-sending the CVOs through the pipeline rather than by re-parsing
  the stored uploads - the old uploads carry no `sequencing_run` metadata, which is the row's key.

## Tests

- `upload/tests/test_import_dragen_tso500_combined_variant_output.py`: the row's columns off the test file, the key
  replacing on re-import of the same run and adding for a new run, the sample and sequencing-sample links landing
  when the arms arrive before and after the file, the parked specimen claim, computed `msi_call` / `tmb_call` with and
  without the settings.
- `analysis/tests/test_case_report.py`: the report picks the row for its samples over a newer row claiming the
  specimen; `measures` keys and the `analysis` snapshot; `tumour_fraction` stays the pathology row.
- `patients/tests/test_specimen_measures.py`: reduced to the Mocha path; the REST tests go with the endpoints.
- `upload/test_data/tso500/README.md` and `upload/AGENTS.md`, `seqauto/AGENTS.md`, `classification/AGENTS.md`
  (the `measure` case-field paragraph) updated; `scripts/vg docs check`.

## Out of scope

- Picking `LibraryQC` by the case's samples rather than the specimen's newest run (same rule, separate change).
- MetricsOutput's `[DNA Expanded Metrics]` / `[RNA Expanded Metrics]`, which `LibraryQC` does not keep either.
- The sapath report template: its `measure` keys keep working unchanged; `tumour_fraction_sequencing` is available
  to it when the lab wants the sequenced estimate on the form.
