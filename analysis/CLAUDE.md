# analysis — agent notes
Owns: Analysis (a DAG of AnalysisNode subclasses filtering Variants), AnalysisEdge, NodeVersion/NodeCache/NodeTask,
the lease-based node scheduler, templates (AnalysisTemplate/Version/Run, AnalysisVariable), VariantTag, node editors, node grid + export.
Start with:
- models/nodes/analysis_node.py — AnalysisNode base (Q composition, versioning, load/counts, save cascade) + NodeVersion/Cache/Task
- models/nodes/sources/*.py, models/nodes/filters/*.py — one module per concrete node
- models/models_analysis.py — Analysis (locking, can_write), AnalysisTemplate/Version/Run
- tasks/analysis_update_tasks.py + tasks/node_update_tasks.py — dispatcher and the per-node celery tasks
- models/nodes/node_utils.py — update_analysis (entry point after any edit), reload_analysis_nodes, tag recounts
- views/nodes/node_view.py + forms/forms_nodes.py — node editors; grids.py + grid_export.py — node grid/export
Patterns here:
- A node is a filter, not a result set. Override `_get_node_q` (or `_get_node_arg_q_dict` when the Q needs an annotation
  alias such as a cohort genotype join); leave `get_queryset` alone. `analysis/models/nodes/analysis_node.py:AnalysisNode.get_arg_q_dict`
  merges parent dicts with the node's own via `AnalysisNode.merge_arg_q_dicts`, and
  `analysis/models/nodes/analysis_node.py:annotate_and_filter_queryset` annotates once per alias then filters — a chain
  of nodes is one SQL query.
- Declare `min_inputs`/`max_inputs` on the subclass (0/0 = source, `AnalysisNode.PARENT_CAP_NOT_SET` = unbounded,
  `analysis/models/nodes/filters/merge_node.py:MergeNode`). Set `uses_parent_queryset = False` and override
  `_get_arg_q_dict_from_parents_and_node` only when the node combines parents itself
  (`analysis/models/nodes/filters/venn_node.py:VennNode`). A pass-through node returns False from `modifies_parents` to
  reuse the parent's counts and cache (`analysis/models/nodes/filters/gene_list_node.py:GeneListNode.modifies_parents`).
- Editing a node: set `queryset_dirty = True`, `save()`, then `analysis/models/nodes/node_utils.py:update_analysis`.
  `analysis/models/nodes/analysis_node.py:AnalysisNode._save` bumps `version`, sets DIRTY, cascades the bump to every
  descendant and creates the NodeVersion row. Everything version-scoped (Redis Q cache keyed on NodeVersion pk, NodeCache,
  NodeTask) cascades from NodeVersion; `analysis/tasks/node_update_tasks.py:delete_analysis_old_node_versions`
  drops the old ones. Move/connect edits go through `analysis/views/views_json.py:NodeUpdate`.
- Scheduling is state-driven, not a prebuilt chain: `analysis/tasks/analysis_update_tasks.py:create_and_launch_analysis_tasks`
  (on `scheduling_single_worker`) calls `analysis/tasks/analysis_update_tasks.py:lease_ready_nodes`, which claims DIRTY
  nodes whose parents are all in `analysis/models/enums.py:NodeStatus.READY_STATUSES` (error parents count — the child
  fails fast with ERROR_WITH_PARENT); each `analysis/tasks/node_update_tasks.py:update_node_task` re-triggers the
  dispatcher when it finishes. `analysis/tasks/analysis_update_tasks.py:dispatch_analysis_backlog` is the periodic sweep.
- Node-specific load work goes in `_load`, returning a dict of fields to persist;
  `analysis/models/nodes/analysis_node.py:AnalysisNode.load` then runs `node_counts` (one aggregate per configured count
  label, plus the exact pk list for nodes under `ANALYSIS_NODE_STORE_ID_SIZE_MAX`) and writes via `AnalysisNode.update`.
  The load's products live on the NodeVersion row: `variant_ids` (the pk list) and `load_data` — `{"counts": {label: count}}`
  plus `{"timings": {phase: seconds}}` and whatever `_get_load_data()` contributes. Every write of those goes through
  `.update()`/raw SQL, so it sets `modified` explicitly.
- Expensive set operations materialise instead of composing: override `use_cache`/`write_cache` to fill a VariantCollection
  (`analysis/models/nodes/filters/intersection_node.py:IntersectionNode.use_cache`); VennNode keeps its own
  `analysis/models/nodes/filters/venn_node.py:VennNodeCache` keyed on the two parent NodeVersions.
- The Duo/Trio/Quad wizards are one view and one template: `analysis/views/views_wizard.py:FamilyWizardView` +
  `analysis/templates/analysis/family_wizard.html`, driven by the subclass's `family_*`/`role_*` attributes. A sample's
  role says which family member it is (Mother/Father/Proband/Sibling), and the sample's sex narrows the roles on offer
  (`analysis/forms/forms.py:FamilyWizardForm`) - a Duo stores its parent's role as the Duo's relationship.
- Node editor = ModelForm subclass of `analysis/forms/forms_nodes.py:BaseNodeForm` + a `NodeView` subclass in views/nodes/
  with `model` set + template `analysis/node_editors/<classname>_editor.html`.
  `analysis/views/views_node.py:get_node_views_by_class` finds the view by `model`, so defining the class registers it;
  `analysis/views/nodes/node_view.py:NodeView.form_valid` does the dirty/save/update_analysis dance for you.
- The sample / patient page's Classify & Report tab is `analysis/classify_report.py` +
  `analysis/views/views_classify_report.py` (it lives here because it is built on VariantTag - analysis may import
  classification, never the other way round). A tagging is in a case's queue when `Tag.requires_classification` and it
  belongs to one of the case's samples: its own `VariantTag.sample` (the study's proband, from
  `analysis/models/nodes/analysis_node.py:AnalysisNode.get_proband_sample` at tag time), else its `VariantTag.patient`
  (`AnalysisNode.get_proband_patient` - a node above sample level names the person while leaving which of their VCFs
  open, so the tagging is on the patient's tab and on each of that patient's sample tabs), else the analysis it was made
  in contains the sample. Carrying the variant is only a display filter for a tagging with no sample - a relative who is
  HET for the proband's variant does not need their own classification, so it never assigns ownership.
  `analysis/variant_tag_operations.py:sample_carries_variant` treats a genotype row as the call when the sample's VCF
  has no GT field (`Sample.has_genotype`) - a fusion caller reports read support, so every zygosity is unknown.
- The queue row's button only opens the dialog when the allele has been curated before - there is something to choose
  between. With nothing previous it links straight to the full create page in a new tab
  (`analysis/classify_report.py:ClassifyQueueRow.full_form_url`), which offers transcript, bucket and sample properly.
  The dialog's own "full form" button goes to the same URL, so starting from scratch is always the full page and the
  dialog is only ever the copy-from-previous shortcut.
- Every route into the full create page goes through the tagging, analysis or not
  (`create_classification_for_variant_tag`, posting to `analysis/views/views.py:create_classification_from_variant_tag`) -
  the analysis alone cannot say which of its taggings a record is for. The page starts on the tagging's own sample, since
  a record without it never reaches the case's queue row or its report. Its sample dropdown is a `ModelSelect2` fed by
  `snpdb/views/views_autocomplete.py:SampleAutocompleteView`, so narrowing `fields['sample'].queryset` only validates the
  POST - what the user sees is narrowed by forwarding the tagging's patient to the autocomplete
  (`analysis/views/views.py:CreateClassificationForVariantTagView`, the same `forward.Const` the specimen autocomplete
  takes).
- The Classify & Report tab's label carries the counts (`analysis/views/views_classify_report.py:classify_report_summary`,
  drawn by `analysis/templates/analysis/classify_report_tab_counts.html`), so a page says whether there is anything to do
  before the tab is opened. It is fetched after render - deciding which taggings are the case's walks every analysis its
  samples are in, which is too much for page load.
- "Classify all" walks `analysis/classify_report.py:ClassifyQueueRow.needs_classification` rows, not unresolved ones -
  a row whose allele has been classified is waiting on "Clear tag", and offering it again made a second record every
  time the wizard was run.
- The queue dialog's previous classifications are filtered to the tag's allele origin bucket
  (`analysis/classify_report.py:tag_allele_origin_bucket`; "Both" means no filter), and an external lab's record is listed
  without "Apply to this sample". Where a tagged allele has nothing of the lab's own, the dialog offers the gene level
  candidates instead - the same `ClassificationConsensus.gene_consensus_groups` rows the create page shows.
- A to-do tagging is resolved against a classification rather than deleted
  (`analysis/variant_tag_operations.py:resolve_variant_tag`), so it stays as the record of what was flagged. That happens
  by itself when the classification is of the tagging's own sample, when the create form was launched from the tagging
  (`analysis/variant_tag_operations.py:resolve_launching_variant_tag` - clicking its "New classification" is the scientist
  saying whose it is), and via the queue row's "Clear tag" button otherwise. A withdrawn `resolved_classification` puts
  the to-do back (`VariantTag.is_resolved`).
- Every "New Classification" button scopes with `snpdb/models/models.py:Tag.classify_queue_qs` rather than naming a tag,
  so a lab's own queue tag is offered and resolved the same way: the tag node editor's Classifications tab
  (`analysis/views/nodes/node_views.py:TagNodeView`), the variant tags grid (`variantopedia/grids.py:VariantTagsColumns`)
  and the Classify & Report tab. Retiring a tag takes it out of the vocabulary, so its taggings stop being to-dos.
- A tagging's identity in an analysis is (variant, tag, analysis, user, sample, patient) - both being the tagged
  node's proband (`AnalysisNode.get_proband`, one walk for the two), worked out before the `get_or_create` in
  `analysis/views/views_json.py:set_variant_tag` and enforced by `varianttag_one_per_person_in_analysis`
  (`nulls_distinct=False`, so an analysis has at most one tagging that names nobody). Sample and patient resolve
  independently, so an extraction level node whose DNA arm has two callers still tags for the person. A tagging never
  changes who it is about: tagging for this proband adds a row rather than taking the tag off a sibling, and the X on a
  pill deletes that one tagging by pk. Older taggings were given their patient by
  `analysis/management/commands/one_off_backfill_variant_tag_patient.py` (#1854).
- The analysis grid draws one pill per tagging, marked with whose it is: a tagging that names someone - a sample, or
  just the patient - gets a solid person marker naming them, boxed as well where that isn't the proband of the node the
  grid is showing (`nodeProbandSampleId` / `nodeProbandPatientId`), and one that names nobody gets a hollow person.
  `variantTags` is `{variant_id: [{id, tag, sample, patient, patient_name, resolved}]}` - one entry per tagging,
  resolution included. The patient's name rides on the tagging (a de-identified patient shows as their code), while a
  sample's comes from `render_analysis_samples_dict` (@see `render_variant_tags_dict`, `VariantGridFormat.tags`,
  `variantTaggingPillOptions` in `grid.js`).
- A resolved tagging is hidden from the work lists: the tags node (`TagNode.include_resolved`, off by default), the
  variant page's tag list and the variant tags page (both on `UserGridConfig.show_hidden_data` under grid name
  `Variant Tags`, shown as a "Show resolved" checkbox). They all filter with
  `analysis/models/models_variant_tag.py:VariantTag.unresolved_q` - the SQL twin of `is_resolved` - never a bare
  `resolved__isnull=True`. The analysis grid keeps the pill, since clicking it is how a tag is removed, and draws it
  as done. Tag stats and the analyses list pills are history and count everything.
- A sample-bound filter node (Zygosity, Allele Frequency, MOI, Gene List) applies to one sample **or** one patient
  (`analysis/models/nodes/cohort_mixin.py:AncestorSampleMixin` - exactly one of `sample`/`patient` is set, both null is
  unset, and the editors pick either through the one `sample_source` control,
  `analysis/forms/forms_nodes.py:AncestorSampleSourceMixin`). `handle_ancestor_input_samples_changed` sets
  the proband sample where there is one, else the proband patient - which is what a group level SampleNode gives it
  (#1855). In patient mode `get_filter_samples()` is every ancestor sample of that patient, so the source node decides
  the reach, and `AncestorSampleMixin._get_filter_samples_arg_q_dict` ORs one `pk IN (subquery)` per sample
  (`analysis/models/nodes/cohort_mixin.py:get_sample_pk_in_q`, the same helper SampleNode's group levels use) - a Q keyed
  on an alias runs as soon as that alias is annotated, so an OR across two samples' aliases has nowhere to hang. A sample
  with nothing to filter on (no GT, no AF) contributes its rows unfiltered rather than emptying the node.
- Load nodes with `AnalysisNode.objects.get_subclass(pk=...)` / `.select_subclasses()` (`analysis/models/nodes/analysis_node.py:NodeInheritanceManager`);
  in views use `analysis/views/analysis_permissions.py:get_node_subclass_or_404`, which enforces
  `analysis/models/models_analysis.py:Analysis.can_write` (locked analyses and template snapshots are read-only).
Gotchas:
- Never call `node.save()` inside a celery task — use `analysis/models/nodes/analysis_node.py:AnalysisNode.update`, a
  conditional UPDATE on (pk, version) that raises `analysis/exceptions.py:NodeOutOfDateException` if the user bumped the
  node meanwhile; `update_node_task` treats that as "exit quietly, a newer version is coming".
- `AnalysisNode.save` locks the Analysis row first (`select_for_update`) so concurrent subtree cascades take NodeVersion
  locks in one order; copy `analysis/signals/source_data_invalidation.py:_bump_nodes` or
  `analysis/models/nodes/node_utils.py:reload_analysis_nodes` when writing a new bulk bump.
- The Q-object cache is on by default (`ANALYSIS_NODE_CACHE_Q`) and keyed only on NodeVersion pk. A test that changes a
  node and expects a different queryset without saving needs `@override_settings(ANALYSIS_NODE_CACHE_Q=False)`.
- A small node answers `get_arg_q_dict` with `Q(pk__in=[...])` from `NodeVersion.variant_ids` (issue #546) - for its
  own grid and export as much as for the children composing it, and ahead of the Redis Q cache, which holds the real
  filter the load ran. A node loaded before those pks were stored falls back to its filter chain - a count mismatch
  between the two paths is a real bug.
- Anything that evaluates a node's Variant queryset runs under
  `analysis/models/nodes/analysis_node.py:node_query_planner_settings` (grid handler, export task, `update_node_task`,
  tag recounts, Venn cache). With the grid's columns selected the query joins 50+ relations, past Postgres's
  `join_collapse_limit` of 8, and the planner keeps the SQL's join order - whole variant table first, the node's
  selective filter last. It is per operation rather than server-wide because the extra planning time regressed
  unrelated queries.
- `load_seconds` is one number - `NodeVersion.load_data["timings"]` is where to look when it's big: seconds per load
  phase (`load`, `live_data_sources`, `counts`, `variant_ids`, `load_data`), also logged as a warning past
  `settings.ANALYSIS_NODE_SLOW_LOAD_SECONDS` and printed by `manage.py profile_analysis_nodes` (#1838).
- The TagNode editor's tag pills are counted at render from `analysis_varianttag` alone
  (`analysis/models/nodes/filters/tag_node.py:TagNode.get_tag_counts`) - taggings in scope, not variants in the node's
  input. Counting them over the node's queryset put a per-tag aggregate over `snpdb_variant` in every TagNode load.
- Tagging does not bump versions: `analysis/signals/signal_handlers.py:variant_tag_create` marks tag nodes dirty and
  `analysis/models/nodes/node_utils.py:update_analysis_tag_node_counts` recounts tag labels in place on the existing
  NodeVersion. It runs concurrently with loads and only computes the tag labels, so it merges into `load_data["counts"]`
  DB-side (`jsonb_set` + `||`) rather than writing the whole dict.
- Counts are sanity-checked at load (`analysis/models/nodes/analysis_node.py:AnalysisNode._raise_or_warn_count_mismatch`):
  any label count > total, or a single-parent node with more variants than its parent. A filter that fans out over a
  multi-valued join (transcript annotation, gene lists) must set `queryset_requires_distinct` or use a subquery.
- A SampleNode threshold only filters when the sample's VCF carries that column
  (`analysis/models/nodes/sources/sample_node.py:SampleNode.get_applied_thresholds`) - PL<=0 against a fusion caller's
  VCF (AD, no PL) emptied the node while the cached stats still answered with its whole count, and the two disagreeing
  is what `_raise_or_warn_count_mismatch` raises on. Anything reading a threshold - query, cache check, method summary -
  goes through there.
- The Duo/Trio/Quad editors' zygosity table is built once per page from a stub node
  (`analysis/models/nodes/sources/trio_node.py:TrioNode.get_zygosity_table_data`), so it cannot see a field's live
  value - anything that depends on one ships as a `{placeholder}` the editor's `updateZygosityTable()` substitutes
  (the mosaic thresholds, #1830). A member with no `other_filters_<member>` key renders a blank cell.
- Changing `Analysis.VERSION_BUMP_FIELDS` (custom columns, default sort) must bump `Analysis.version`
  (`analysis/forms/forms.py:AnalysisForm`): node grids and editors are `cache_page`d under the analysis version in the
  URL (`analysis/views/views_node.py:node_view`).
- A worker that dies mid-LOADING perma-fails the node on the next sweep (`lease_ready_nodes`) on the assumption it OOM'd
  the box; only QUEUED nodes are re-leased, up to `analysis/tasks/node_update_tasks.py:MAX_NODE_ATTEMPTS`.
- Templates run by cloning the snapshot then setting AnalysisVariable-bound fields in toposort order
  (`analysis/models/models_analysis.py:AnalysisTemplateRun.populate_arguments`); a source node's editor gets the variable
  widget only when `analysis.template_type == TEMPLATE` (`analysis/views/nodes/node_view.py:NodeView.get_form`).
- Renaming, moving or removing a class whose instances are pickled into Redis - enum choices inside a node's
  `arg_q_dict` Q objects, dataclasses, model subclasses - needs a `CACHE_VERSION` bump in
  `variantgrid/settings/components/default_settings.py`. Cache keys are pks (`analysis/models/nodes/analysis_node.py`
  keys node Q dicts on `node_version.pk`), so old entries survive a deploy and fail on unpickle with
  `AttributeError: Can't get attribute 'OldName'`; the bump flushes every deployment's cache at once.
Tests:
- `analysis/tests/utils.py:AnalysisSetupMixin` gives `cls.analysis` + `cls.grch37` with a fake annotation version
  (`annotation/fake_annotation.py:get_fake_annotation_version`); samples/cohorts/trios/quads/pedigrees from
  `snpdb/tests/utils/fake_cohort_data.py:create_fake_cohort` and friends;
  `snpdb/tests/utils/vcf_testing_utils.py:slowly_create_test_variant` for a handful of real Variants.
- `analysis/tests/test_urls.py:Test` is the URLTestCase (every analysis/node URL incl. editors and grid exports, owner vs
  non-owner) — run it after touching urls.py, a view signature or an editor template. Celery is eager under URLTestCase,
  so `update_analysis` executes inline. `analysis/tests/test_scheduler.py` covers lease/backoff/reclaim.
- A test that creates a classification needs `@override_settings(CLINGEN_ALLELE_REGISTRY_LOGIN=None)` — autopopulate
  asks ClinGen for the variant and `snpdb/tests/utils/mock_clingen_api.py` raises on an HGVS it has no recorded response
  for, which surfaces as a 500 from the create POST rather than as a ClinGen error
  (`analysis/tests/test_classify_report.py`).
- `manage.py profile_analysis_nodes --analysis <id> --rerun --explain` times every node's queryset and dumps EXPLAIN plans
  (also `--sample/--trio/--cohort` synthetic runs) — run it before and after changing a Q.
Deep reference: __analysis_readme.md · claude/research/analysis.md · claude/maps/models.md#analysis
