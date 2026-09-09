# Sample-bound filter nodes that apply to a patient

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-09
Status: in progress

Issue: [#1855](https://github.com/SACGF/variantgrid/issues/1855). Follows #1854, which gave every node a
patient as well as a sample (`analysis/models/nodes/analysis_node.py:AnalysisNode.get_proband`).

## Problem

Four filter nodes hang their genotype filter off one `sample` FK, auto-set from the ancestors' proband
(`analysis/models/nodes/cohort_mixin.py:AncestorSampleMixin.handle_ancestor_input_samples_changed`):
Zygosity, Allele Frequency, Mode of Inheritance and Gene List (its sample QC gene list panel).

A group level SampleNode (extraction / specimen / patient, `analysis/models/nodes/sources/sample_node.py:SampleNode`)
feeds several samples in - a TSO 500 DNA extraction is a small variant caller and a CNV caller, each its own
VCF and sample. `get_proband_sample()` then has no single answer, so the filter node stays unset and the user
has to pick one caller's sample. Picking one is wrong, not just inconvenient: a zygosity filter keyed on the
small variant sample's genotype column is NULL on every CNV row, so the node silently drops the other caller's
variants.

## Decision

Each of the four nodes applies to **either one sample or one patient**. In patient mode the filter applies to
every ancestor sample that belongs to that patient, and the node's filter is the OR of the per-sample filters -
the shape the group level SampleNode already produces. The proband patient (`AnalysisNode.get_proband_patient`)
is what auto-sets it where the proband sample is ambiguous.

"Every ancestor sample" is the scope, rather than every sample of the patient, so an extraction level node
above gives the DNA arm only and a patient level node gives everything - the source node decides the reach,
the filter follows it. The analysis genome build is already enforced by the ancestors
(`patients/sample_grouping.py:get_sample_group`).

## Data

`patient` joins `sample` on each of the four models. Exactly one of the two is set; both null means unset.
The form, `_set_sample` and `_set_patient` keep that invariant, so no mode column is needed.

```python
# analysis/models/nodes/filters/zygosity_node.py
class ZygosityNode(AncestorSampleMixin, AnalysisNode):
    sample = models.ForeignKey(Sample, null=True, on_delete=SET_NULL)
    patient = models.ForeignKey(Patient, null=True, blank=True, on_delete=SET_NULL)
    zygosity = models.CharField(max_length=1, choices=ZygosityNodeZygosity.CHOICES, null=True)
    exclude = models.BooleanField(default=False)

# analysis/models/nodes/filters/allele_frequency_node.py
class AlleleFrequencyNode(AncestorSampleMixin, AnalysisNode):
    sample = models.ForeignKey(Sample, null=True, on_delete=SET_NULL)
    patient = models.ForeignKey(Patient, null=True, blank=True, on_delete=SET_NULL)

# analysis/models/nodes/filters/moi_node.py
class MOINode(AncestorSampleMixin, AnalysisNode):
    sample = models.ForeignKey(Sample, null=True, blank=True, on_delete=SET_NULL)
    patient = models.ForeignKey(Patient, null=True, blank=True, on_delete=SET_NULL)
    ...  # unchanged

# analysis/models/nodes/filters/gene_list_node.py
class GeneListNode(AncestorSampleMixin, GeneCoverageMixin, AnalysisNode):
    sample = models.ForeignKey(Sample, null=True, blank=True, on_delete=SET_NULL)
    patient = models.ForeignKey(Patient, null=True, blank=True, on_delete=SET_NULL)
    sample_gene_list = models.ForeignKey(SampleGeneList, null=True, blank=True, on_delete=SET_NULL)
    ...  # unchanged
```

One migration in `analysis/migrations/` (0144_filter_node_patient), four `AddField`s. Existing nodes keep their
sample; nothing is backfilled - a node that is set stays as it is, and one that is unset picks the patient up
the next time its ancestors change.

## The mixin

`analysis/models/nodes/cohort_mixin.py:AncestorSampleMixin` owns the two-way choice. Its docstring says the
model needs `sample` and `patient` fields and that one at most is set.

- `_set_sample(sample)` clears `patient`; new `_set_patient(patient)` clears `sample`. GeneListNode's
  `_set_sample` override keeps doing its `sample_gene_list` work and `_set_patient` clears `sample_gene_list`.
- `get_filter_patient() -> Optional[Patient]`: `patient` in patient mode, otherwise the sample's patient by
  `patients/sample_grouping.py:get_patient_for_source` at `SampleSourceLevel.SAMPLE` (the sample may be linked
  directly or through its extraction). MOINode's patient panel and the MOI editor read this.
- `get_filter_samples() -> list[Sample]`: `[sample]` in sample mode; in patient mode the ancestor samples
  (`_get_ancestor_samples`) whose patient is `self.patient`, sorted by pk. Every path that used
  `self.sample` for a genotype join goes through this.
- `_get_sample()` (the `SampleMixin` hook) returns None in patient mode, as the group SampleNode does, so the
  single-cohort machinery in `CohortMixin` stays out of the way; `_get_cohorts_and_sample_visibility_for_node`
  and `_get_annotation_kwargs_for_node` cover every filter sample instead, copying
  `SampleNode._get_cohorts_and_sample_visibility_for_node` / `SampleNode._get_annotation_kwargs_for_node`.
- `_get_filter_samples_arg_q_dict(per_sample) -> dict[Optional[str], dict[str, Q]]` is the one place the OR
  is built. `per_sample(sample)` returns that sample's arg_q_dict keyed on its own alias, exactly what each
  node builds today. Sample mode returns it as is, so a single-sample node's query is byte for byte what it is
  now. Patient mode wraps each sample's dict in a `pk__in` subquery annotated with only that sample's
  genotype join and ORs them under the `None` key, hashed with `_get_node_q_hash()`. That is
  `SampleNode._get_sample_pk_q`, so that method's body moves to a module function in `cohort_mixin.py`,
  `get_sample_pk_in_q(node, sample, arg_q_dict) -> Q`, and SampleNode calls it. The reason it is subqueries and
  not one Q over several aliases is `analysis/models/nodes/analysis_node.py:annotate_and_filter_queryset`: a
  Q keyed on an alias runs as soon as that alias is annotated, and an OR across two aliases has nowhere to
  hang until both are.
- `handle_ancestor_input_samples_changed`: a set sample or patient that is no longer reachable from the
  ancestors is cleared (the existing rule, extended to the patient). When nothing is set: the proband sample if
  there is one, else the proband patient, else the single ancestor sample. Single-sample analyses therefore
  keep setting the sample and stay unchanged.
- `_get_configuration_errors`: a patient that is not the patient of any ancestor sample, or one with no
  ancestor samples, is an error, worded like the existing sample one. Zygosity and Allele Frequency require a
  sample or patient ("No sample or patient selected."); MOI and Gene List stay optional as they are.
- `analysis/signals/source_data_invalidation.py:handle_sample_pre_delete` bumps patient mode nodes as well:
  the four `<node>__patient=` terms for the deleted sample's patient join the existing `__sample=` ones, and
  the comment above them says why.

## Per node

**Zygosity** (`analysis/models/nodes/filters/zygosity_node.py:ZygosityNode`). The per-sample filter is the
Q it builds now on `get_cohort_genotype_alias_and_field("zygosity")`, `exclude` applied inside it, so in
patient mode "exclude HET" is "this sample's rows that are not HET", unioned. A sample whose VCF has no GT
(`Sample.has_genotype`, a fusion caller) has nothing to filter on and contributes its rows unfiltered -
zygosity IN all four codes, as `SampleNode._get_sample_arg_q_dict` does. The multiple-hit branch works off
the parent queryset and is untouched. Method summary lists which samples were filtered and which passed
through.

**Allele Frequency** (`analysis/models/nodes/filters/allele_frequency_node.py:AlleleFrequencyNode`). Per
sample: `NodeAlleleFrequencyFilter.get_sample_arg_q_dict(self, sample)` where `sample.has_allele_frequency`,
else the sample's rows unfiltered (the rule `SampleNode._get_sample_arg_q_dict` applies). That rule now
applies in sample mode too - a single sample with no AF column passes through with the method summary saying
so, where the node used to empty itself. `modifies_parents` is true when the filter has a restricted range.

**Mode of Inheritance** (`analysis/models/nodes/filters/moi_node.py:MOINode`). The patient panel's terms come
from `get_filter_patient()`. The zygosity-and-genes OR is built per sample by the helper; with no filter
samples the gene-only branch runs as now. `analysis/views/nodes/node_views.py:MOINodeView.get_context_data`
passes the patient gene/disease data for `get_filter_patient()`, and the editor's "From Patient" panel keys
on the picker rather than on `#id_sample`.

**Gene List** (`analysis/models/nodes/filters/gene_list_node.py:GeneListNode`). In patient mode the sample QC
panel's gene lists are the active sample gene list of every filter sample, resolved when `get_gene_lists()`
is called (deduplicated, samples without one skipped) and `sample_gene_list` stays null - `_get_node_q_hash`
already folds in the gene list pks, so nothing stored goes stale. Sample mode keeps `sample_gene_list` as the
stored choice. `GeneCoverageMixin` iterates `get_samples()` and needs no change.

## Editor

One picker replaces the sample `<select>` in the four editors, following
`analysis/forms/forms_nodes.py:SampleNodeForm` (`source` carries `"<level>:<pk>"`): a form field `applies_to`
with choices `sample:<pk>` for each ancestor sample and `patient:<pk>` for each patient those samples resolve
to, labelled "<patient> (all N samples)". The choices are a plain `Select` - the set is the ancestors', small
and already permission-checked. `save()` unpacks it through `_set_sample` / `_set_patient`, and
`get_analysis_variable_field("applies_to")` answers `"patient"` or `"sample"` from what the instance holds, so
an analysis template variable binds to the FK that is set. Node name and method summary say
"Patient X (3 samples)" in patient mode.

## Tests

In `analysis/tests/test_sample_node_levels.py`, on `SampleNodeLevelsTestCase` (the TSO 500 patient: a DNA
extraction with two callers, an RNA arm, a blood draw, a sample linked straight to the patient), a new class
for the filter nodes:

1. A ZygosityNode under the extraction node auto-sets the patient, sample stays null; under a single-sample
   node it sets the sample and its arg_q_dict is keyed on the sample's alias, as today.
2. HET under the extraction node gives the small variant and CNV HET calls; exclude HET gives the rest of both
   callers' rows; the RNA, blood and unlinked samples' variants are out of reach at extraction level and in
   reach under a patient level node.
3. An allele frequency range applies per sample; a sample whose VCF has no AF passes through.
4. Reconnecting the node under another patient's source clears the patient.
5. The picker round-trips `patient:<pk>` and `sample:<pk>`, and an AnalysisVariable on it binds to the right FK.
6. MOI patient panel terms come from the node's patient; Gene List sample QC panel unions the samples' lists.

`analysis/tests/test_urls.py` covers the four editors. `scripts/vg tests --explain` names the rest
(`test_zygosity_nodes`, `test_allele_frequency_node`, `test_moi_node`, `test_gene_list_node`,
`test_clone_nodes`, `test_serializers`).

## Docs

`analysis/CLAUDE.md`: the `AncestorSampleMixin` line under Patterns says a filter node applies to a sample or a
patient and how the OR is built; the invalidation comment in `source_data_invalidation.py` is updated in place.
