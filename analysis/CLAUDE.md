# analysis — agent notes
Owns: Analysis (a DAG of AnalysisNode subclasses filtering Variants), AnalysisEdge, NodeVersion/NodeCount/NodeCache/NodeTask,
the lease-based node scheduler, templates (AnalysisTemplate/Version/Run, AnalysisVariable), VariantTag, node editors, node grid + export.
Start with:
- models/nodes/analysis_node.py — AnalysisNode base (Q composition, versioning, load/counts, save cascade) + NodeVersion/Cache/Count/Task
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
  descendant and creates the NodeVersion row. Everything version-scoped (Redis Q cache keyed on NodeVersion pk, NodeCount,
  NodeCache, NodeTask) cascades from NodeVersion; `analysis/tasks/node_update_tasks.py:delete_analysis_old_node_versions`
  drops the old ones. Move/connect edits go through `analysis/views/views_json.py:NodeUpdate`.
- Scheduling is state-driven, not a prebuilt chain: `analysis/tasks/analysis_update_tasks.py:create_and_launch_analysis_tasks`
  (on `scheduling_single_worker`) calls `analysis/tasks/analysis_update_tasks.py:lease_ready_nodes`, which claims DIRTY
  nodes whose parents are all in `analysis/models/enums.py:NodeStatus.READY_STATUSES` (error parents count — the child
  fails fast with ERROR_WITH_PARENT); each `analysis/tasks/node_update_tasks.py:update_node_task` re-triggers the
  dispatcher when it finishes. `analysis/tasks/analysis_update_tasks.py:dispatch_analysis_backlog` is the periodic sweep.
- Node-specific load work goes in `_load`, returning a dict of fields to persist;
  `analysis/models/nodes/analysis_node.py:AnalysisNode.load` then runs `node_counts` (one aggregate per configured count
  label, plus the exact pk list for nodes under `ANALYSIS_NODE_STORE_ID_SIZE_MAX`) and writes via `AnalysisNode.update`.
- Expensive set operations materialise instead of composing: override `use_cache`/`write_cache` to fill a VariantCollection
  (`analysis/models/nodes/filters/intersection_node.py:IntersectionNode.use_cache`); VennNode keeps its own
  `analysis/models/nodes/filters/venn_node.py:VennNodeCache` keyed on the two parent NodeVersions.
- Node editor = ModelForm subclass of `analysis/forms/forms_nodes.py:BaseNodeForm` + a `NodeView` subclass in views/nodes/
  with `model` set + template `analysis/node_editors/<classname>_editor.html`.
  `analysis/views/views_node.py:get_node_views_by_class` finds the view by `model`, so defining the class registers it;
  `analysis/views/nodes/node_view.py:NodeView.form_valid` does the dirty/save/update_analysis dance for you.
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
- Small parents are inlined as `Q(pk__in=[...])` from the pks stored in NodeCount
  (`analysis/models/nodes/analysis_node.py:AnalysisNode.get_small_parent_arg_q_dict`, issue #546); a parent loaded before
  those pks were stored falls back to the subquery — a count mismatch between the two paths is a real bug.
- Tagging does not bump versions: `analysis/signals/signal_handlers.py:variant_tag_create` marks tag nodes dirty and
  `analysis/models/nodes/node_utils.py:update_analysis_tag_node_counts` recounts tag labels in place on the existing NodeVersion.
- Counts are sanity-checked at load (`analysis/models/nodes/analysis_node.py:AnalysisNode._raise_or_warn_count_mismatch`):
  any label count > total, or a single-parent node with more variants than its parent. A filter that fans out over a
  multi-valued join (transcript annotation, gene lists) must set `queryset_requires_distinct` or use a subquery.
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
