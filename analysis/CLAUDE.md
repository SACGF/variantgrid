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
  plus whatever `_get_load_data()` contributes (TagNode snapshots its editor's `tag_counts` there, as counting them in
  global mode is slow). Every write of those goes through `.update()`/raw SQL, so it sets `modified` explicitly.
- Expensive set operations materialise instead of composing: override `use_cache`/`write_cache` to fill a VariantCollection
  (`analysis/models/nodes/filters/intersection_node.py:IntersectionNode.use_cache`); VennNode keeps its own
  `analysis/models/nodes/filters/venn_node.py:VennNodeCache` keyed on the two parent NodeVersions.
- Node editor = ModelForm subclass of `analysis/forms/forms_nodes.py:BaseNodeForm` + a `NodeView` subclass in views/nodes/
  with `model` set + template `analysis/node_editors/<classname>_editor.html`.
  `analysis/views/views_node.py:get_node_views_by_class` finds the view by `model`, so defining the class registers it;
  `analysis/views/nodes/node_view.py:NodeView.form_valid` does the dirty/save/update_analysis dance for you.
- The sample / patient page's Classify & Report tab is `analysis/classify_report.py` +
  `analysis/views/views_classify_report.py` (it lives here because it is built on VariantTag - analysis may import
  classification, never the other way round). A tagging is in a case's queue when `Tag.requires_classification` and it
  belongs to one of the case's samples: its own `VariantTag.sample` (the study's proband, from
  `analysis/models/nodes/analysis_node.py:AnalysisNode.get_proband_sample` at tag time), else the analysis it was made
  in contains the sample. Carrying the variant is only a display filter for a tagging with no sample - a relative who is
  HET for the proband's variant does not need their own classification, so it never assigns ownership.
- The queue dialog's previous classifications are filtered to the tag's allele origin bucket
  (`analysis/classify_report.py:tag_allele_origin_bucket`; "Both" means no filter), and an external lab's record is listed
  without "Apply to this sample". Where a tagged allele has nothing of the lab's own, the dialog offers the gene level
  candidates instead - the same `ClassificationConsensus.gene_consensus_groups` rows the create page shows.
- A to-do tagging is resolved against a classification rather than deleted
  (`analysis/variant_tag_operations.py:resolve_variant_tag`), so it stays as the record of what was flagged. That happens
  by itself when the classification is of the tagging's own sample and via the queue row's "Clear tag" button otherwise.
  A withdrawn `resolved_classification` puts the to-do back (`VariantTag.is_resolved`).
- Every "New Classification" button scopes with `snpdb/models/models.py:Tag.classify_queue_qs` rather than naming a tag,
  so a lab's own queue tag is offered and resolved the same way: the tag node editor's Classifications tab
  (`analysis/views/nodes/node_views.py:TagNodeView`), the variant tags grid (`variantopedia/grids.py:VariantTagsColumns`)
  and the Classify & Report tab. Retiring a tag takes it out of the vocabulary, so its taggings stop being to-dos.
- A resolved tagging is hidden from the work lists: the tags node (`TagNode.include_resolved`, off by default), the
  variant page's tag list and the variant tags page (both on `UserGridConfig.show_hidden_data` under grid name
  `Variant Tags`, shown as a "Show resolved" checkbox). They all filter with
  `analysis/models/models_variant_tag.py:VariantTag.unresolved_q` - the SQL twin of `is_resolved` - never a bare
  `resolved__isnull=True`. The analysis grid keeps the pill, since clicking it is how a tag is removed, and draws it
  as done. Tag stats and the analyses list pills are history and count everything.
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
- Small parents are inlined as `Q(pk__in=[...])` from `NodeVersion.variant_ids`
  (`analysis/models/nodes/analysis_node.py:AnalysisNode.get_small_parent_arg_q_dict`, issue #546); a parent loaded before
  those pks were stored falls back to the subquery — a count mismatch between the two paths is a real bug.
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
- Changing `Analysis.VERSION_BUMP_FIELDS` (custom columns, default sort) must bump `Analysis.version`
  (`analysis/forms/forms.py:AnalysisForm`): node grids and editors are `cache_page`d under the analysis version in the
  URL (`analysis/views/views_node.py:node_view`).
- A worker that dies mid-LOADING perma-fails the node on the next sweep (`lease_ready_nodes`) on the assumption it OOM'd
  the box; only QUEUED nodes are re-leased, up to `analysis/tasks/node_update_tasks.py:MAX_NODE_ATTEMPTS`.
- Templates run by cloning the snapshot then setting AnalysisVariable-bound fields in toposort order
  (`analysis/models/models_analysis.py:AnalysisTemplateRun.populate_arguments`); a source node's editor gets the variable
  widget only when `analysis.template_type == TEMPLATE` (`analysis/views/nodes/node_view.py:NodeView.get_form`).
Tests:
- `analysis/tests/utils.py:AnalysisSetupMixin` gives `cls.analysis` + `cls.grch37` with a fake annotation version
  (`annotation/fake_annotation.py:get_fake_annotation_version`); samples/cohorts/trios/quads/pedigrees from
  `snpdb/tests/utils/fake_cohort_data.py:create_fake_cohort` and friends;
  `snpdb/tests/utils/vcf_testing_utils.py:slowly_create_test_variant` for a handful of real Variants.
- `analysis/tests/test_urls.py:Test` is the URLTestCase (every analysis/node URL incl. editors and grid exports, owner vs
  non-owner) — run it after touching urls.py, a view signature or an editor template. Celery is eager under URLTestCase,
  so `update_analysis` executes inline. `analysis/tests/test_scheduler.py` covers lease/backoff/reclaim.
- `manage.py profile_analysis_nodes --analysis <id> --rerun --explain` times every node's queryset and dumps EXPLAIN plans
  (also `--sample/--trio/--cohort` synthetic runs) — run it before and after changing a Q.
Deep reference: __analysis_readme.md · claude/research/analysis.md · claude/maps/models.md#analysis
