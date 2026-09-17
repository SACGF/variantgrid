# #647 — VCF page: bulk set sample BAM/CRAM files from a pattern

Written by Claude Opus 5 (claude-opus-5), 2026-09-16
Status: draft

[#647](https://github.com/SACGF/variantgrid/issues/647) (milestone SA Path VG 4): "Make a new tab on VCF page where you
can use a pattern to create lots of BAMs for a VCF (like old functionality)".

## Background

The VCF page used to have a "Bulk Set Bam files" link: a text box defaulting to `/data/{name}.bam` and a JS
`setBamPath()` that filled a BAM path column in the VCF's sample formset from each sample's name. Commit `7c5ee7c03`
(2022, "#2 - multiple bams per sample") removed it when the single `Sample.bam_file_path` column became its own model
(`git show 7c5ee7c03 -- snpdb/templates/snpdb/data/view_vcf.html`); the page has since become
`snpdb/templates/snpdb/data/view_vcf_cohort.html`, shared with the cohort page.

Today files are only added one sample at a time, on the sample page's "Files" tab
(`snpdb/views/views_data.py:sample_files_tab`, `snpdb/templates/snpdb/data/sample_files_tab.html`, an
inline formset `SampleFilesFormSet` in `snpdb/forms.py`), or automatically by seqauto on VCF import
(`upload/vcf/vcf_import.py`, a `get_or_create` of the sequencing sample's single BAM). A 96-sample VCF from outside
seqauto means 96 visits to the Files tab.

Rows are read by `snpdb/models/models_vcf.py:Sample.get_bam_files` (the sample page's IGV links),
`analysis/models/nodes/analysis_node.py:AnalysisNode.get_bams_dict` (IGV from analysis grids), and the display-side
prefix mapping in `snpdb/sample_file_path.py:get_bam_paths_and_user_data_paths`. #1805 (IGV link per sample on
variant details) will read them too.

## Data

No model change and no migration. The rows this feature writes are the existing model, unchanged:

```python
class SampleFilePath(models.Model):
    sample = models.ForeignKey(Sample, on_delete=CASCADE)
    file_type = models.CharField(max_length=1, null=True, blank=True, choices=SampleFileType.choices)
    label = models.TextField(null=True, blank=True)
    file_path = models.TextField()
```

(`snpdb/models/models_vcf.py:SampleFilePath`, `file_type` from `snpdb/models/models_enums.py:SampleFileType`.)

There is no unique constraint on `(sample, file_type, file_path)`, and existing databases may already hold duplicates
(manual Files-tab entries plus seqauto's `get_or_create`), so dedupe is done by the bulk writer reading the existing
set, rather than by adding a constraint that would need a data cleanup migration first.

Per-sample preview/save result, a dataclass in `snpdb/sample_file_path.py` (members only):

```python
@dataclass
class SampleFilePathResolution:
    sample: Sample
    file_path: Optional[str]     # None when the pattern could not be resolved for this sample
    error: Optional[str]         # eg "No patient_code" - shown in the preview, sample skipped on save
    already_exists: bool         # a SampleFilePath with this sample/file_type/file_path is already stored
```

## Behaviour (defaults chosen)

- **Where:** a new "Sample Files" tab on the VCF page, VCF only (not the cohort page that shares the template),
  AJAX-loaded on first click via `data-href` like the sample page's Files tab.
- **Pattern syntax:** Python %-style with named keys, the same syntax and placeholders as the user's grid sample label
  template (`snpdb/models/models_vcf.py:Sample._get_sample_formatter_params`): `sample_id`, `sample`, `patient_id`,
  `patient_code`, `patient`, `specimen_id`, `specimen`, plus a new `vcf_sample_name` (the column name in the VCF file,
  which is usually what a pipeline named the BAM, and does not change when a user renames the sample). A literal `%`
  is written `%%`. Initial value `/data/%(vcf_sample_name)s.bam`.
- **Form fields:** pattern (required), file type (BAM or CRAM only; BED/VCF stay on the per-sample tab), label
  (optional, applied to rows this save creates).
- **Form validation** (whole pattern, before any sample is looked at): at least one `%(key)s` placeholder (a pattern
  with none would give every sample the same file); every key is one of the known placeholders; trial-formatting
  against a dict holding every key catches bad conversions (`%(sample)d`, a dangling `%`).
- **Per-sample resolution:** keys whose value is missing or blank for this sample (no patient, no specimen, patient
  with no `patient_code`) are removed before formatting, so a pattern using them raises `KeyError` and the sample gets
  `error = "No <key>"`. Blank values are an error rather than a path with an empty segment.
- **Flow:** the form has two submit buttons, Preview and Save.
  - *Preview* renders a table, one row per sample (VCF column order): sample, existing files (type, label, path),
    resolved path, and a status of New / Already present / error message. Nothing is written.
  - *Save* re-resolves on the server (the preview is not trusted), writes the New rows, and skips errored and
    already-present samples. The message says "Added N, M already present, K skipped" with the skipped sample names.
    The tab reloads in place showing the updated existing-files column.
- **Additive only:** never deletes or edits an existing row; a present path with a different label is left alone.
  Removing or relabelling stays on the per-sample Files tab.
- **No filesystem check:** the web server often cannot see the sequencing storage (hence `UserDataPrefix`), so paths
  are stored as typed.
- **Permissions:** the tab and its existing-files table are visible to anyone who can view the VCF
  (`VCF.get_for_user`). The form is only rendered, and any POST (preview or save) only accepted, when
  `vcf.can_write(request.user)` (`library/django_utils/guardian_permissions_mixin.py:GuardianPermissionsMixin.check_can_write`
  raises `PermissionDenied`). Write on the VCF implies write on each of its samples
  (`snpdb/models/models_vcf.py:Sample.can_write` defers to the VCF), so no per-sample check is needed.

## Implementation

### 1. Placeholder

`snpdb/models/models_vcf.py:Sample._get_sample_formatter_params` gains `"vcf_sample_name": self.vcf_sample_name`. The
grid sample label template picks it up for free; `snpdb/forms.py:SettingsOverrideForm._validate_sample_formatter_func`
builds its dummy `Sample` without it, so give that dummy a `vcf_sample_name` too.

### 2. Resolution and bulk write — `snpdb/sample_file_path.py`

Add a module docstring (it has none): owns sample file path display mapping and bulk creation from a pattern; entry
points `get_bam_paths_and_user_data_paths`, the resolver and the writer below.

- The `SampleFilePathResolution` dataclass above.
- A `PATTERN_KEYS` tuple (the placeholder names) and a `validate_sample_file_path_pattern(pattern)` that raises
  `ValueError` for the form-level rules, so the form's `clean_pattern` is a thin wrapper.
- `resolve_sample_file_paths(vcf, pattern, file_type) -> list[SampleFilePathResolution]`: one queryset
  `vcf.sample_set.select_related("patient", "extraction__specimen").prefetch_related("samplefilepath_set")`, ordered by
  pk (VCF column order). `already_exists` is computed from the prefetched rows, so preview is a fixed number of
  queries whatever the sample count.
- `create_sample_file_paths(resolutions, file_type, label) -> list[SampleFilePath]`:
  `SampleFilePath.objects.bulk_create` of the resolutions with a `file_path`, no `error` and not `already_exists`.
  Called inside `transaction.atomic` from the view.

### 3. Form — `snpdb/forms.py`

`VCFSampleFilesPatternForm(forms.Form)`: `pattern` (`CharField`, `TextInput`, help text listing the placeholders and
`%%`), `file_type` (`ChoiceField` restricted to `SampleFileType.BAM` / `SampleFileType.CRAM`), `label`
(`CharField`, `required=False`). `clean_pattern` calls `validate_sample_file_path_pattern`.

### 4. View and URL

`snpdb/views/views_data.py`, `vcf_sample_files_tab(request, vcf_id)`, next to `sample_files_tab`:

- `vcf = VCF.get_for_user(request.user, vcf_id)`; `has_write_permission = vcf.can_write(request.user)`.
- GET: unbound form, `resolutions = None`, table of samples with existing files (reuse the resolver's queryset shape
  so this is not N+1).
- POST: `vcf.check_can_write(request.user)`, bind the form. `request.POST["action"]` is `preview` or `save` (the two
  submit buttons share `name="action"`; jquery.form's `ajaxForm` posts the clicked button). Valid form → resolve;
  on `save`, `create_sample_file_paths` then re-resolve so the table shows the new rows as present, and add the
  counts message (`messages.add_message`, surfaced by `{% update_django_messages %}`). Invalid form →
  `add_save_message(request, False, ...)` and re-render with errors.
- Renders *snpdb/templates/snpdb/data/vcf_sample_files_tab.html*.

`snpdb/urls.py`: `path('vcf_sample_files_tab/<int:vcf_id>', views_data.vcf_sample_files_tab,
name='vcf_sample_files_tab')`, beside `sample_files_tab`.

### 5. Templates

- `snpdb/templates/snpdb/data/view_vcf_cohort.html`: inside the existing `{% if vcf %}` block of the `nav-tabs`, add
  `<li class="nav-item"><a class="nav-link" data-toggle="tab" role="tab" href="#vcf-sample-files" data-href="{% url 'vcf_sample_files_tab' vcf.pk %}">Sample Files</a></li>`
  after "VCF Info", and an empty `<div class="tab-pane" role="tabpanel" id="vcf-sample-files"></div>` in the
  `tab-content` `{% if vcf %}` section. The generic `.nav-tabs a[data-href][data-toggle="tab"]` processor in
  `variantgrid/static_files/default_static/js/global.js` loads it.
- New *vcf_sample_files_tab.html*: copy the shape of `snpdb/templates/snpdb/data/sample_files_tab.html` - a wrapper
  div whose script does `{% update_django_messages %}` and `$('form#vcf-sample-files-form').ajaxForm({target: parent})`
  so preview and save reload the tab in place. Form rendered with crispy (`form_helper.horizontal`) and the Preview
  (`btn-secondary`) / Save (`btn-primary`) buttons only when `has_write_permission`; otherwise "You can view but not
  modify this VCF." Below it a plain Bootstrap `table` of samples: sample (`{% preview sample %}` or a link to
  `view_sample`), existing files, and when `resolutions` is set the resolved path and a status badge
  (success New / secondary Already present / danger error). Save is shown after a preview too, so the flow is
  type pattern → Preview → check → Save; Save is not disabled before preview.

### 6. Close the write gap on the per-sample Files tab

`snpdb/views/views_data.py:sample_files_tab` saves the POSTed formset after only `Sample.get_for_user` (a read check).
Add `sample.check_can_write(request.user)` at the top of its POST branch, the same rule the new tab enforces.

## Tests

New *snpdb/tests/test_vcf_sample_file_paths.py*, fixtures from `snpdb/tests/utils/fake_cohort_data.py:create_fake_cohort`
(one VCF, samples `proband`/`mother`/`father`, permissions for the creating user). Give `proband` a `Patient` with a
`patient_code` and `mother` a patient without one; leave `father` with none.

1. `test_resolve_pattern` - `/data/%(vcf_sample_name)s.bam` resolves for all three; `/data/%(patient_code)s.bam`
   resolves for proband and errors ("No patient_code") for mother and father; `validate_sample_file_path_pattern`
   rejects a pattern with no placeholder and one with an unknown key.
2. `test_create_is_additive` - pre-create a BAM row for proband at the resolved path with label "old" and a BED row
   for mother; create twice with the same pattern. Proband's row is untouched (one row, label "old"), mother gets a
   BAM row alongside the BED, and the second create adds nothing.
3. `test_write_permission_required` - a second user given only the VCF read perm (`assign_perm(VCF.get_read_perm(), user, vcf)`)
   POSTs `action=save` to `vcf_sample_files_tab` and to `sample_files_tab`: both 403, no rows created.

`snpdb/tests/test_urls.py`: add `('vcf_sample_files_tab', {"vcf_id": cls.vcf.pk}, 200)` to
`PRIVATE_OBJECT_URL_NAMES_AND_KWARGS` (render for owner, 404/403 for non-owner).

## Verification

- `scripts/vg tests --explain --run`, plus `python3 manage.py test --keepdb snpdb.tests.test_vcf_sample_file_paths snpdb.tests.test_urls snpdb.tests.test_sample_formatter`.
- `python3 manage.py vg page /snpdb/view_vcf/<id> --queries` before and after (the page itself should not change
  query count) and `python3 manage.py vg page /snpdb/vcf_sample_files_tab/<id> --queries` on a many-sample VCF: the
  count must not grow with the number of samples.
- In a browser (`run` skill): Preview then Save on a multi-sample VCF; the sample page Files tab shows the new rows and
  its IGV links use them.

## Out of scope

- Editing or deleting files in bulk, and BED/VCF types (per-sample Files tab).
- Checking the files exist, or finding them by globbing a directory.
- Bulk-setting across several VCFs or a cohort.
