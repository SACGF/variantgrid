# Library QC from DRAGEN's MetricsOutput.tsv: Amplifications / Variants / Fusions / Quality / Fail have a source

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-22
Status: landed 181b8df16 - sapath settings and the NGS-pipelines plan alongside; the sapath report-template step 6 remains

[sapath#455](https://github.com/SACGF/variantgrid_sapath/issues/455), phase 2 of #454
(`claude/plans/sapath454_multi_report_tick_measures_plan.md`, landed). Phase 1 put the measure behind the TMB, MSI and
Purity checkboxes on the TSO 500 Build report form. The other five flags mean "did DRAGEN's library QC pass for that
category", which lives in a file we do not import: *MetricsOutput.tsv*, written per pair beside the
CombinedVariantOutput. This plan imports it, stores one row per library per QC category, and lets a `case_field`
start its tick from that row. The client half is done - [variantgrid_api#24](https://github.com/SACGF/variantgrid_api/issues/24), released as variantgrid_api 1.8.0;
the pipeline half is *claude/plans/sapath455_metrics_output_upload_plan.md* in ../NGS-pipelines.

## What the file says

The fixture is `upload/test_data/tso500/ExampleSample_2600000001/ExampleSample_2600000001_MetricsOutput.tsv`, its
quirks in `upload/test_data/tso500/README.md`. It is the CVO's layout - banner line, `[Section]`, tab padding - so
`upload/tso500/dragen_combined_variant_output_parser.py:read_combined_variant_output` already reads it whole.

| Section | Shape | Becomes |
|---|---|---|
| `[Header]` | key/value: Output Date / Time, Workflow Version | `measured_date`, `method` - as `upload/tso500/dragen_combined_variant_output_records.py:measured_date` reads the CVO |
| `[Run QC Metrics]` | metric, LSL, USL, one `Value` column | the run's, not a library's - not stored |
| `[Analysis Status]` | first header cell empty, then one column per sample: `COMPLETED_ALL_STEPS` TRUE/FALSE, `FAILED_STEPS`, `STEPS_NOT_EXECUTED` | `completed` on every row of that sample |
| `[DNA Library QC Metrics]` | `CONTAMINATION_SCORE` USL 1457 | category `DNA` |
| `[... for Small Variant Calling and TMB]` | `MEDIAN_INSERT_SIZE` >=70, `MEDIAN_EXON_COVERAGE` >=150, `PCT_EXON_50X` >=90 | `SMALL_VARIANT_TMB` - the **Variants** flag |
| `[... for MSI]` | `USABLE_MSI_SITES` >=40 | `MSI` |
| `[... for CNV]` | `GENE_SCALED_MAD` <=0.134, `MEDIAN_BIN_COUNT_CNV_TARGET` >=1 | `CNV` - the **Amplifications** flag |
| `[... for GIS]` | `PCT_TARGET_HRD_50X` >=50 | `GIS` |
| `[RNA Library QC Metrics]` | `MEDIAN_CV_GENE_500X` <=0.93, `TOTAL_ON_TARGET_READS` >=9M, `MEDIAN_INSERT_SIZE` >=80 | `RNA` - the **Fusions** flag and the **Quality** caveat |
| `[DNA Expanded Metrics]`, `[RNA Expanded Metrics]`, 2.6's `[DNA Library Sex Metrics]` | no guidelines | not stored |

Every QC section has a column per sample and `NA` down the arm it does not apply to, so one file writes rows for both
arms. A metric passes when `LSL <= value <= USL` with `NA` as no bound - the same rule the lab's own release summary
applies (../NGS-pipelines *scripts/tso500_metrics_parser.py* `run_tso_qc_checks`), and the guideline numbers are the
ones the report's methods paragraph quotes (40 sites, 1457, 0.134, 150x, 70bp, 90%; RNA 0.93, 9M, 80bp). So the file's
own guideline is the policy and there is no threshold setting: a lab whose number differs from Illumina's asks for one
then, the way `TSO500_MSI_CALL_BANDS` arrived.

The sample columns are DRAGEN's Sample_IDs - exactly the CVO's `DNA Sample ID` / `RNA Sample ID`, a `Sample.vcf_sample_name`
and a `SequencingSample.sample_name`. The file names no pair, patient or specimen, so it cannot create the chain the way
the CVO does; it claims an extraction by the accession inside the sample ID
(`upload/tso500/dragen_combined_variant_output_records.py:SAMPLE_ID_ACCESSION_PATTERN`) and, where that extraction is
not there yet, parks the claim exactly as a VCF sample does.

## Data

### `LibraryQC` - one library, one QC category

`seqauto/models/models_seqauto.py`, beside the other QC models; the enum in `seqauto/models/models_enums.py`. It started
in patients beside `SpecimenMeasure` and moved (2026-09-22): its identity is a sequencing library (run + pair), the specimen
is a nullable claim the way `SequencingSample.extraction` is, and seqauto already imports `patients.models`.

```python
class LibraryQCCategory(models.TextChoices):
    """ What a caller's library QC section vouches for - vendor neutral, DRAGEN's sections map onto it in the parser """
    DNA = 'D', 'DNA library'                      # contamination
    SMALL_VARIANT_TMB = 'V', 'Small variants and TMB'
    MSI = 'M', 'MSI'
    CNV = 'C', 'CNV'
    GIS = 'G', 'GIS'
    RNA = 'R', 'RNA library'                      # fusions and splice variants


# The key a report template's case_fields names a category by, as MEASURE_CONTEXT_KEYS does for measures
LIBRARY_QC_CONTEXT_KEYS = {
    LibraryQCCategory.DNA: "dna",
    LibraryQCCategory.SMALL_VARIANT_TMB: "small_variant_tmb",
    LibraryQCCategory.MSI: "msi",
    LibraryQCCategory.CNV: "cnv",
    LibraryQCCategory.GIS: "gis",
    LibraryQCCategory.RNA: "rna",
}


class LibraryQC(ExtractionMatchMixin, TimeStampedModel):
    """ One QC category of one sequenced library, as the caller judged it - what 'the assay succeeded for X' means.
        Keyed on the library (the caller's sample name) rather than the extraction: a repeat sequencing is a new
        library with its own QC, and the extraction it belongs to may not be accessioned yet """
    sample_name = models.TextField()            # the file's column - a Sample.vcf_sample_name / SequencingSample.sample_name
    category = models.CharField(max_length=1, choices=LibraryQCCategory.choices)
    passed = models.BooleanField(null=True)     # every metric within its guideline; None where the section is all NA for this arm
    completed = models.BooleanField(null=True)  # [Analysis Status] COMPLETED_ALL_STEPS for the library, the same on each of its rows
    metrics = models.JSONField(default=dict)    # {metric: {"value", "unit", "lsl", "usl", "passed"}} - the section, so 'which number' stays answerable
    method = models.TextField(blank=True)       # 'DRAGEN TSO500 MetricsOutput 2.1.1.4' off [Header]
    measured_date = models.DateTimeField(null=True, blank=True)
    file_upload = models.ForeignKey(FileUpload, null=True, on_delete=SET_NULL)  # provenance, as SpecimenMeasure keeps source_payload
    user = models.ForeignKey(User, null=True, on_delete=SET_NULL)

    class Meta:
        unique_together = ("sample_name", "category")   # a re-analysis of the same library replaces its rows
```

**Superseded by the revision below**, which is what landed: a file column turned out to be a pair carrying both arms,
so the row is keyed on (run, pair, category) and claims a Specimen rather than an Extraction.

`ExtractionMatchMixin` (`patients/models.py:ExtractionMatchMixin`) gives it `extraction` (nullable), `extraction_reference`,
`extraction_match_status`, `extraction_match_error`, `extraction_match_date`.
`patients/tasks/extraction_matching_tasks.py:reconcile_pending_extractions` adds `LibraryQC` to the models it
re-resolves (its `user` is the row's own `user`), so a file that lands before the CVO or Mocha has made the
extraction attaches itself on the next reconcile, as a Sample does. No Guardian mixin: rows are only ever read through
an extraction's specimen (the specimen page, the case report), which has already checked the patient's permission,
and the extraction FK can be null until matched.

Migration: seqauto `0048_libraryqc`.

`patients/models_enums.py:MatchStatus` is unchanged. `Extraction` and `SpecimenMeasure` are unchanged.

### The upload record

`upload/models/models.py`: `UploadedDragenTSO500MetricsOutput(UploadData)` with `file_upload` OneToOne, the shape of
the other single-shot satellites - it is what `upload/uploaded_file_type.py:get_upload_data_for_uploaded_file` needs to
call the upload processed, which is what the API's sha256 de-duplication keys on
(`upload/views/views_rest.py:APIFileUploadView`). `get_data_url` is the specimen page of the first matched row, or None.

`upload/models/models_enums.py:UploadedFileTypes` gains `DRAGEN_TSO500_METRICS_OUTPUT = 'Q', 'DRAGEN TSO500 Metrics Output'`,
and `upload/file_type_icons.py` an entry for it (`test_file_type_icons` checks).

### The `case_field` keys

A bool `case_field` names a library QC category with `qc` the way it names a measure with `measure`; `tick_when` takes
two more rules, and a field carries exactly one of `measure` / `qc`:

```python
{"key": "assay_success_amplifications", "label": "Amplifications", "type": "bool", "default": True,
 "group": "assay_success",
 "qc": "cnv",                              # a LIBRARY_QC_CONTEXT_KEYS value
 "tick_when": {"passed": True}}

# tick_when rules against a qc row (a list ticks when any holds, as for measures):
{"passed": True} | {"passed": False}       # LibraryQC.passed
{"completed": False}                       # LibraryQC.completed - the run did not finish for the library
```

## Implementation

1. **Parser.** *upload/tso500/dragen_metrics_output_parser.py* (new): the section reader is the CVO's
   (`read_combined_variant_output` - a `read_sections` alias is fine, the format is the module's), `can_process_file`
   off a banner pattern ending `Metrics Output` beside the CVO's `FIRST_LINE_PATTERN`. It maps the section names
   above to `LibraryQCCategory`, reads `Metric (UOM)` into name and unit, `LSL Guideline` / `USL Guideline` as floats
   or None for `NA`, and yields per sample column `{category: {metric: {value, unit, lsl, usl, passed}}}` plus
   `completed` from `[Analysis Status]` (whose header's first cell is empty - the sample names start at column two).
   A section or metric it has no name for is skipped, since 2.6 adds and drops metrics (README). `passed` for a
   category is all-of its metrics' `passed`, and None when every value in the section is `NA` for that sample.

2. **Factory and task.** `upload/import_task_factories/import_task_factories.py`:
   `DragenTSO500MetricsOutputImportTaskFactory(ImportTaskFactory)` - single shot, not a VCF pipeline, so it subclasses
   `upload/import_task_factories/import_task_factory.py:ImportTaskFactory` like `PatientRecordsImportTaskFactory`;
   extensions `['tsv']`, ability 1000 off `can_process_file` (a tsv would otherwise be a gene list),
   `get_metadata_keys` the default empty set (the file has no coordinates, so a `genome_build` is a 400 - the
   client note in variantgrid_api#24), enabled under `settings.VARIANT_GENE_LEVEL_ENABLED` like its siblings so
   the capabilities endpoint (`variantgrid/views_rest.py:CapabilitiesView`) lists `dragen_tso500_metrics_output`
   on the deployments that take TSO 500 files. *upload/tasks/import_dragen_tso500_metrics_output_task.py* (new):
   `ImportDragenTSO500MetricsOutputTask(ImportTask)` returning the row count; registered at module bottom, listed in
   `variantgrid/settings/components/celery_settings.py` `CELERY_IMPORTS` and routed to `WEB_WORKERS` beside
   `ImportPatientRecords`.

3. **Records.** *upload/tso500/dragen_metrics_output_records.py* (new): per sample column, the extraction reference off
   `SAMPLE_ID_ACCESSION_PATTERN` (a column carrying no accession is a `SimpleVCFImportInfo`-style message on the
   pipeline and no rows), `resolve_reference(Extraction, ...)` scoped to the uploading user, and one `LibraryQC` per
   category `update_or_create`d on `(sample_name, category)` with `apply_extraction_match` - matched, pending or
   needs attention, the rows are written either way. Both arms come out of one file; the DNA arm gets rows for the
   five DNA categories and the RNA arm for `RNA`, the other arm's rows being None-passed and not written.
   After the writes, `reconcile_pending_extractions.delay()` as `patients/views_rest.py` does after an extraction
   lands, so a MetricsOutput that arrived first is matched as soon as the CVO makes the extraction.

4. **Case report.** `classification/report/case_report_context.py`: `specimen_library_qc(specimen) -> dict[str, LibraryQC]`
   beside `specimen_measures` - the specimen's extractions' rows, newest `measured_date` per category (a repeat
   library supersedes), keyed by `LIBRARY_QC_CONTEXT_KEYS`. `analysis/views/views_classify_report.py:case_report_build_dialog`
   passes it to `_case_values_for_form` and `_measure_notes`.
   `classification/models/classification_report_models.py`: `measure_tick` grows a sibling `library_qc_tick(tick_when, qc)`
   (or one `tick_for(field, measures, library_qc)` dispatching on which key the field carries), `describe_tick_when`
   words the two new rules ("ticked when the library passed QC", "ticked when the run did not complete for the
   library"), and `validate_case_fields` requires exactly one of `measure` / `qc`, a known category key, and only
   `passed` / `completed` rules against a `qc`.
   `analysis/templates/analysis/case_report_build_dialog.html`: a `qc` field shows `passed` / `failed` / `no QC`
   after its label, and the note under the group lists the category's metrics against their guideline -
   `CNV: GENE_SCALED_MAD 0.059 (<= 0.134), MEDIAN_BIN_COUNT_CNV_TARGET 6.4 (>= 1) · passed` - so a scientist sees
   which number failed, the way phase 1 shows a measure's threshold.

5. **Specimen page.** `patients/templates/patients/view_specimen.html` gets a Library QC table under Measures: one row
   per (extraction, category) with passed / completed / measured / method, and the metrics expanded in a title or
   collapsed row. `seqauto/admin.py` registers the model.

6. **sapath.** *sapath/tso500_case_template.py*: Amplifications `qc: cnv, passed`; Variants
   `qc: small_variant_tmb, passed`; Fusions `qc: rna, passed`; Quality caveat `qc: rna, [passed False]`; Fail caveat
   `qc: dna, [passed False, completed False]`; TMB, MSI, Purity and Deletions as they are. Migration `0018` updates
   the "TSO 500" row from `CASE_FIELDS` as `0017` does. No env setting changes.

7. **Tests.** *upload/tests/test_import_dragen_tso500_metrics_output.py* (new), on the fixture: the parser reads every
   category for the right arm with the other arm None; a metric outside its guideline fails the category; `NA`
   guideline is no bound; `COMPLETED_ALL_STEPS` FALSE lands as `completed=False`; an existing extraction matches, a
   missing one parks Pending and `reconcile_pending_extractions` attaches it once the extraction exists; a
   re-analysis replaces the rows; the factory claims the file over the gene list and the CVO factories.
   `analysis/tests/test_case_report.py`: the `passed` and `completed` rules, a missing row falling back to
   `default`, newest library winning, and `validate_case_fields` rejecting a field with both `measure` and `qc`.
   `variantgrid/tests/test_capabilities.py`: the new type is listed.

8. **Docs.** `classification/AGENTS.md` case report bullet: the `qc` key and its two rules. `upload/AGENTS.md`: one
   gotcha line - the MetricsOutput import is single shot, claims extractions by accession and parks the rest.
   `upload/test_data/tso500/README.md` upload table row for *MetricsOutput.tsv*: no metadata at all. `scripts/vg map`
   after the model and task land.

## Revision, 2026-09-22: one run-level file, not one per pair

Found while doing the pipeline side ([the issue comment](https://github.com/SACGF/variantgrid_sapath/issues/455#issuecomment-5771815899)):
the lab's *tso500_run_wrapper.py* rewrites every *MetricsOutput.tsv* in place - the run-level one and each pair's -
before the results reach TAU (a sex metrics block is inserted, `[Analysis Status]` rows are padded, the research-use
text is stripped). The one untouched DRAGEN copy is the run-level *Results/MetricsOutput_orig.tsv* the wrapper keeps,
which has one column per sample for **every pair in the run**. So the pipeline now sends that, once per sequencing
run, with the run named in upload metadata:

```python
vg_api.upload_file("Results/MetricsOutput_orig.tsv", path=None,
                   metadata={"sequencing_run": "<SequencingRun.name>"},
                   file_type=UploadFileType.DRAGEN_TSO500_METRICS_OUTPUT)
```

An unknown metadata key is a 400 the client does not skip, so the importer has to accept `sequencing_run` before any
server lists `dragen_tso500_metrics_output` in its capabilities. What lands so far is untouched by this except where
the steps below say.

### Data change

A real run-level 2.6.2 file settled the column model: a column is the CVO's **Pair ID**
(`5_C0000001_FCUP_2600000001` - sequencing number, patient C-number, initials, ten-digit accession, no container
suffix) and carries **both arms** - the DNA sections and the RNA sections have values in the same column. So a
column names no extraction. It claims the Specimen whose accession it ends in, resolve-only (accessioning a case is
the CombinedVariantOutput's job), parked when the specimen is not there yet.

The lab writes the Pair ID inconsistently - the long form above, or the patient code alone (`C0000001`), which the
NGS-pipelines fixture uses. So `TSO500_PAIR_ID_PATIENT_CODE_REGEX` reads both
(`^(?:\d+_)?(?P<patient_code>[^_]+)(?:_|$)`), and a bare column, having no accession, takes one off the run: the
current sample sheet's `SequencingSample` names carry the code and the accession together
(`1_TSO_DNAHRD_C17817_2517114977C_B4`). A column nothing on the run is named for keeps its rows with the claim
parked saying so. Nothing else changes in the CVO - its arms' sample IDs still carry accession and container.

The row's identity is the pair on its run, which the column name alone does not give: the lab's pair IDs carry a
leading sequencing number, but nothing guarantees a re-sequenced pair is renamed. The run name is always present
(the key is required for this file type), so it is part of the key; the SequencingRun is linked when registered and
left null when it has not been, the way the CVO's `link_to_sequencing_run` tolerates a missing `SequencingSample`.
A nullable FK cannot be the unique key (Postgres treats nulls as distinct), which is why the name is stored as well.
The arms' samples are linked by the CVO, so no `SequencingSample` is kept here.

```python
class LibraryQC(TimeStampedModel):
    sequencing_run_name = models.TextField()    # the upload's 'sequencing_run' metadata - a SequencingRun.name
    pair_id = models.TextField()                # the file's column - the CVO's 'Pair ID'
    specimen_reference = models.TextField()     # the ten-digit accession inside the pair ID
    specimen = models.ForeignKey(Specimen, null=True, blank=True, on_delete=SET_NULL)
    specimen_match_status = models.CharField(max_length=1, choices=MatchStatus.choices, null=True, blank=True)
    specimen_match_error = models.TextField(null=True, blank=True)
    specimen_match_date = models.DateTimeField(null=True, blank=True)
    sequencing_run = models.ForeignKey(SequencingRun, null=True, blank=True, on_delete=SET_NULL)
    sequencing_sample = models.ForeignKey(SequencingSample, null=True, blank=True, on_delete=SET_NULL)  # the arm's sheet row
    category = models.CharField(max_length=1, choices=LibraryQCCategory.choices)
    nucleic_acid = models.CharField(max_length=1, choices=NucleicAcid.choices)  # the arm the category is about
    passed = models.BooleanField(null=True)
    completed = models.BooleanField(null=True)
    metrics = models.JSONField(default=dict)    # {metric: {"value", "unit", "lsl", "usl", "passed", "guideline_source"}}
    method = models.TextField(blank=True)
    measured_date = models.DateTimeField(null=True, blank=True)
    file_upload = models.ForeignKey("upload.FileUpload", null=True, on_delete=SET_NULL)
    user = models.ForeignKey(User, null=True, on_delete=SET_NULL)

    class Meta:
        unique_together = ("sequencing_run_name", "pair_id", "category")   # a re-analysis of the run replaces its rows
```

Each row is one arm, and links that arm's `SequencingSample` by the sheet's own `Pair_ID` / `Sample_Type` columns
(`seqauto/models/models_seqauto.py:sequencing_sample_for_pair`), posted per sample as `SequencingSampleData` - which needs
sapath's `SEQAUTO_SAMPLE_SHEET_EXTRA_COLUMNS` to list them and the pipeline to send them (../NGS-pipelines
*claude/plans/tso500_sample_sheet_pair_columns_plan.md*). Both arms of a pair are linked across its rows: the five DNA
categories to the DNA sample, RNA to the RNA sample. The link is left null and
`patients/tasks/extraction_matching_tasks.py:link_library_qc_to_sequencing_samples` fills it once the sheet is registered
or re-sent. The linked arm's name also supplies a bare column's accession before the patient-code scan does.

`ExtractionMatchMixin` is not used - the claim is on the Specimen, so `LibraryQC.apply_specimen_match` mirrors it over
the `specimen_*` columns, and `patients/tasks/extraction_matching_tasks.py:reconcile_pending_extractions` re-resolves
those rows beside the extraction ones. `seqauto/migrations/0048_libraryqc.py` is unpushed, so it is regenerated
rather than added to.

Two more things the real file settled:

- Every metric of a known section counts towards its category, not just the ones the table above lists - 2.6.2 adds
  `PCT_CHIMERIC_READS` (USL 8) to the small-variant section and `EXCESSIVE_TF` (USL 0) to GIS. Only the section names
  are keyed on.
- The file's own guideline is still the policy, with one escape hatch: `TSO500_LIBRARY_QC_GUIDELINES`
  (`variantgrid/settings/components/default_settings.py`), `{(section, metric): (lsl, usl)}` - keyed on the section
  too, since `MEDIAN_INSERT_SIZE` is in two of them. The lab's methods paragraph quotes 9M `TOTAL_ON_TARGET_READS`
  where the 2.6.2 file's LSL is 2,500,000. Each metric records the bounds applied and `guideline_source`.

### Steps

1. **Metadata.** `upload/upload_metadata.py`: `SEQUENCING_RUN = "sequencing_run"`, a validator (non-empty string;
   existence is not checked - the run is normally registered first by the pipeline's `sequencing_run` step, and an
   unregistered one is linked by nobody rather than rejected), added to `_VALIDATORS`. The factory's
   `get_metadata_keys` returns `{SEQUENCING_RUN}`, and the import fails with a clear message when the upload carries
   none - a run-level file with no run named cannot be keyed.
2. **Parser.** `upload/tso500/dragen_metrics_output_parser.py`: a column is a `pair_id`, every metric of a known
   section is read, and `guideline_for(section, metric, lsl, usl)` applies `TSO500_LIBRARY_QC_GUIDELINES`.
   `CATEGORY_NUCLEIC_ACID` says which arm each category is about.
3. **Records.** `upload/tso500/dragen_metrics_output_records.py:write_library_qc` takes the run name, resolves
   `SequencingRun` by name, reads each column's specimen accession off `PAIR_ID_ACCESSION_PATTERN` (the trailing ten
   digits) or, for a bare-patient-code column, off the run's current sample sheet (`specimen_reference_from_run`),
   `resolve_reference(Specimen, ...)` scoped to the uploading user, and `update_or_create`s on the three-column key.
   A column nothing names a specimen for keeps its rows with `LibraryQC.park_specimen_claim` saying so - a run holds
   many pairs, and the controls and other assays sharing its flowcell name none.
   `UploadedDragenTSO500MetricsOutput.get_data_url` points at the sequencing run page when the run resolved, else the
   first matched specimen.
4. **Case report.** `classification/report/case_report_context.py:specimen_library_qc` reads the specimen's rows
   directly (`filter(specimen=specimen)`), newest `measured_date` per category - a re-sequenced pair now keeps its
   earlier run's rows, which is the point of the key.
5. **Specimen page.** The Library QC table's columns are run (linked when `sequencing_run` is set), pair, arm,
   category, QC, metrics, method, measured.
6. **Fixture.** A run-level `upload/test_data/tso500/MetricsOutput_orig.tsv` (new) written from a real 2.6.2 file's
   layout with synthetic identifiers: `5_C0000001_FCUP_2600000001` (the ExampleSample pair, within guideline),
   `7_C0000002_ABCD_2600000002` (one DNA metric outside guideline, `COMPLETED_ALL_STEPS` FALSE with a `FAILED_STEPS`
   entry), `9_0PRI_2600000003` (no C-number, as controls are named), `C0000004` (the bare patient-code
   form, whose accession comes off the run's sample sheet) and `11_NTC` (naming nothing at all). The per-pair
   file stays for the parser tests and for variantgrid_api's copy. README: describe the run-level file, and the upload
   table row for `MetricsOutput` becomes metadata `sequencing_run` **required**.
7. **Tests.** In `upload/tests/test_import_dragen_tso500_metrics_output.py`: rows keyed on run and pair and resolved
   to the specimen; the second pair's failing metric and unfinished run; the control pair parked; the bare-code column resolved
   through the run's sample sheet and parked without one; a registered run linked and an unregistered one null; the same file re-imported replaces rather
   than duplicates; the same pair on another run is a second row; a guideline override changes `passed`; an upload
   without `sequencing_run` fails with the message. `upload/tests/test_api.py`: `sequencing_run` accepted for this
   type and a `genome_build` on it still a 400.
8. **Docs.** `upload/AGENTS.md` gotcha line: keyed on the run and the pair, `sequencing_run` metadata required.
   variantgrid_api: the 1.8.0 CHANGELOG says "no metadata" - a docs-only fix there (README `upload_file` metadata keys
   and the changelog line), no release needed.

## Order across the repos

1. variantgrid_api#24 - done, released as 1.8.0 on 2026-09-22 (a docs-only follow-up for the metadata key).
2. This plan (server) including the revision above, so `api/v1/capabilities` lists `dragen_tso500_metrics_output`
   with `sequencing_run` already accepted - the name the client gates on.
3. ../NGS-pipelines plan - written on branch `sapath455_metrics_output`, sending the run-level file; inert against a
   server that does not list the type.
4. sapath step 6 can land with step 2 or after it - the template row only ticks once rows exist.
