# #29 — Clone node count settings, show sample owner

[#29](https://github.com/SACGF/variantgrid/issues/29) is a 2020 grab-bag of analysis fixes. Most of it
has landed or gone stale (see the audit comment on the issue). Two items are still live and are small
enough to do together:

1. **Clone doesn't copy node count settings** — a real bug in `Analysis.clone()`.
2. **Show the owner on the sample page** — the analysis settings and VCF pages already show theirs.

The two are independent; either can ship alone.

---

## 1. `Analysis.clone()` loses node count settings

### How node counts are configured

A per-analysis override lives in two models (`analysis/models/models_analysis.py:370`):

```python
class AnalysisNodeCountConfiguration(models.Model):
    analysis = models.OneToOneField(Analysis, on_delete=CASCADE)

class AnalysisNodeCountConfigRecord(AbstractNodeCountSettings):
    node_count_config = models.ForeignKey(AnalysisNodeCountConfiguration, on_delete=CASCADE)
```

`Analysis.get_node_count_types()` (`models_analysis.py:245`) reads the records in `sort_order`, and
falls back to `BuiltInFilters.DEFAULT_NODE_COUNT_FILTERS` when there's no configuration at all.
`set_node_count_types()` (`models_analysis.py:258`) `get_or_create`s the configuration and hands the
record set to `AbstractNodeCountSettings.save_count_configs_from_array()`
(`snpdb/models/models_user_settings.py:724`), which deletes and recreates.

Three states are therefore distinct, and all three need to survive a clone:

| Source state | `get_node_count_types()` returns |
|---|---|
| No `AnalysisNodeCountConfiguration` | `DEFAULT_NODE_COUNT_FILTERS` |
| Configuration with records | those records, in `sort_order` |
| Configuration with zero records | `[]` — the user explicitly turned every count off |

The third is reachable: the node counts tab posts a `node_counts` string and
`_analysis_settings_node_counts_tab` (`analysis/views/views.py:1022`) passes whatever it splits
straight to `set_node_count_types`. So copying only the records would silently promote
"all counts off" back to the defaults — the configuration's *existence* is the state to carry.

### The gap

`Analysis.clone()` (`models_analysis.py:277`) copies nodes (`save_clone`), edges (`AnalysisEdge`) and
`AnalysisVariable` rows, then returns. Nothing copies the node count configuration, so every clone
falls into the first row of the table above.

Note the shape of the method: it mutates `self` into the copy —

```python
analysis_id = self.pk
analysis_copy = self
analysis_copy.pk = None  # to save as new
```

so after line 290 `self` is the *new* analysis and `analysis_id` is the only handle on the original.
The copy must read the source config through `analysis_id`, the way the `AnalysisEdge` and
`AnalysisVariable` loops already do.

### The change

In `Analysis.clone()`, after the `AnalysisVariable` loop and inside the `disable_auditlog()` block:

```python
if source_config := AnalysisNodeCountConfiguration.objects.filter(analysis_id=analysis_id).first():
    # An empty record set is a real setting (every count switched off), so copy the
    # configuration itself - falling back to the defaults would be a different analysis
    config_copy = AnalysisNodeCountConfiguration.objects.create(analysis=analysis_copy)
    for record in source_config.analysisnodecountconfigrecord_set.all().order_by("sort_order"):
        record.pk = None
        record.node_count_config = config_copy
        record.save()
```

`AnalysisNodeCountConfiguration` is already imported at module scope — it's defined in this file.

### What else this fixes

`Analysis.clone()` is the copy primitive behind every analysis-duplicating path, so one change covers
them all:

- `analysis_clone` (`analysis/views/views_json.py:56`) — the user-facing "clone" the issue reported.
- `AnalysisTemplate.save_snapshot` (`models_analysis.py:544`) — a template version's snapshot now
  captures the node counts the template author configured.
- `AnalysisTemplate.clone` (`models_analysis.py:564`).
- `AnalysisTemplateRun.create` (`models_analysis.py:625`) — analyses generated from a template now
  inherit the template's node counts rather than the build-wide defaults.
- `analysis_create_default_templates` (`analysis/management/commands/`), and the form path at
  `analysis/forms/forms.py:190`.

The template paths are the bigger win: a lab that sets up "show me ClinVar and OMIM counts" on a
template currently loses it on every run.

Existing analyses are unaffected — this only changes what future clones get, so no migration.

### Test

`analysis/tests/test_clone_nodes.py` has no node count coverage. Add a small `TestCase` there using
`AnalysisSetupMixin` (`analysis/tests/utils.py`), which gives `cls.analysis` and `cls.grch37` without the
heavy node fixtures `TestCloneAnalysisNodes` builds — node counts don't need nodes.

Set the state explicitly in each test rather than relying on the fixture: `set_defaults_and_save()`
(`models_analysis.py:232`) calls `set_node_count_types()` when the user has a
`NodeCountSettingsCollection`, so whether the fixture analysis starts with a configuration depends on
the test user's settings. Cover the three states, since the empty-record case is the one that's easy to
regress:

- source with records → clone returns the same labels in the same order.
- source with a configuration and zero records → clone returns `[]`, not the defaults.
- source with no configuration → clone returns the defaults.

Compare on `get_node_count_types()`, which is what the settings tab and node loading both read.

---

## 2. Owner on the sample page

### Where owners already show

- **Analysis settings** — `AnalysisForm.Meta` (`analysis/forms/forms.py:227`) lists `user` in `fields`
  and in `read_only`, so `ROFormMixin` (`library/forms.py:58`) disables it and
  `{% crispy form form_helper.horizontal %}` renders it. Already done.
- **VCF** — `VCFForm.Meta.fields` (`snpdb/forms.py:269`) includes `user`, rendered by
  `{% crispy vcf_form form_helper.horizontal %}` (`view_vcf.html:356`). Already visible.
- **Grids** — the VCF, sample and analysis list grids all carry an owner column
  (`snpdb/grids.py:80,144`, `analysis/grids.py:466`).

That leaves the sample detail page.

### Why the sample page has none

`Sample` (`snpdb/models/models_vcf.py:323`) has no `user` field — a sample's owner is its VCF's
uploader, `sample.vcf.user`. `SampleForm` (`snpdb/forms.py:362`) uses `exclude`, so there's no model
field to surface.

### The change

`SampleForm` already has the pattern for exactly this — `genome_build` is a plain `CharField` declared
on the form, initialised from `self.instance.vcf.genome_build`, and listed in `read_only`:

```python
genome_build = forms.CharField()  # From VCF.genome_build
```

Add an `uploaded_by` alongside it:

- Declare `uploaded_by = forms.CharField()  # From VCF.user` next to `genome_build`
  (`snpdb/forms.py:363`).
- Add `'uploaded_by'` to `Meta.read_only` (`snpdb/forms.py:369`).
- Set `self.fields['uploaded_by'].initial = self.instance.vcf.user` in `__init__`, beside the existing
  `genome_build` initial (`snpdb/forms.py:385`).

Label it "Uploaded by" to match the samples grid column (`snpdb/grids.py:144`) — the same fact, so it
should read the same in both places. `CharField` renders the username via `User.__str__`.

`_field_order()` (`snpdb/forms.py:388`) puts `genome_build` first, then everything that isn't an
extraction-match field, then the extraction-match fields. `uploaded_by` lands in the middle group with
no change needed; if it should sit next to `genome_build`, extend the leading list to
`["genome_build", "uploaded_by"]` and the `editable` comprehension's exclusion to match.

No template change — `view_sample.html:164` renders the whole form through crispy.

### Test

`SampleForm` is worth one test: that `uploaded_by` initialises from `sample.vcf.user`, since it's our
derivation across a relation rather than a field declaration. The `disabled=True` behaviour comes from
`ROFormMixin` and Django, which are already covered.

---

## Related, tracked separately

[#1809](https://github.com/SACGF/variantgrid/issues/1809) covers the `User` foreign keys on forms across
the project: six render a plain `Select` listing every user account, and `user_autocomplete` already
exists as the fix.

That audit found something that constrains this plan: `ROFormMixin` sets `disabled = True`, and a
disabled Django `Select` still renders every `<option>`. So `AnalysisForm.user` being `read_only` shows
the owner but also ships the full user list — reading the owner off a disabled picker is not the pattern
to copy.

This is why item 2 above adds `uploaded_by` to `SampleForm` as a read-only `CharField` (following
`genome_build`) rather than a `ModelChoiceField`: a `CharField` renders the one username and nothing
else. #1809 can revisit the existing `read_only` `ModelChoiceField`s on their own.
