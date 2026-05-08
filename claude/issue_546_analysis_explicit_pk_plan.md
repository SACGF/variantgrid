# Issue #546 — Speed up small jobs by using explicit variant IDs

GitHub: https://github.com/SACGF/variantgrid/issues/546

## What the issue is asking

When a parent node holds only a small number of variants (think: 200 from a diagnostics sample), the parent's full filter chain is needlessly re-executed every time a child or sibling node composes its own queryset. PostgreSQL's planner makes poor choices when these small results are wrapped in `pk IN (subquery)` against the 40M-row variant table.

The shortcut: when the parent's `count` is small, materialise its PKs once and substitute the parent's contribution to the child's `arg_q_dict` with a literal `Q(pk__in=[…])`. The child sees a parent that filters by PK only — no FilteredRelation joins, no annotation aliases — and Postgres plans a tight bitmap-or over `snpdb_variant_pkey`.

## What is already in place

`MergeNode` already does this for its own input parents (`analysis/models/nodes/filters/merge_node.py:132-145`):

```python
if settings.ANALYSIS_NODE_MERGE_STORE_ID_SIZE_MAX and \
        parent.count is not None and \
        parent.count <= settings.ANALYSIS_NODE_MERGE_STORE_ID_SIZE_MAX:
    variant_ids = AnalysisNode.get_parent_pks(parent)
    q = Q(pk__in=variant_ids)
    arg_q_dict = {None: {q: q}}
```

Backed by a `cache_memoize`-d helper keyed on `(parent.pk, parent.version)` (`analysis/models/nodes/analysis_node.py:530-539`), with a 15-minute TTL and a defensive ceiling at the same setting.

Setting: `ANALYSIS_NODE_MERGE_STORE_ID_SIZE_MAX = 1000` (`variantgrid/settings/components/default_settings.py:333`).

## What is left

Generalise the small-parent substitution from `MergeNode` to all descendant nodes that consume a parent through `get_parent_arg_q_dict`. The hook point is `AnalysisNode.get_single_parent_arg_q_dict` (`analysis/models/nodes/analysis_node.py:329-341`) for single-input nodes, plus the `MergeNode` multi-input override.

## Safety condition — annotation alias scope

Parents like `CohortNode`, `SampleNode`, `TrioNode`, `PedigreeNode` (via `CohortMixin._get_annotation_kwargs_for_node`, `analysis/models/nodes/cohort_mixin.py:52-56`) emit `arg_q_dict` keyed on a non-None alias such as `cohortgenotype_alias`. Children — and any descendant subtree under them — may produce their own Q objects keyed on that same alias (e.g. allele-frequency filters at `cohort_mixin.py:174`, VCF filter codes at `cohort_mixin.py:185-207`, zygosity filters at `cohort_mixin.py:137-155`).

The substitution from `parent.get_arg_q_dict()` to `{None: {q: Q(pk__in=…)}}` is safe when the resulting merged `arg_q_dict` (parent contribution + child contribution) carries no non-None outer keys — i.e. neither the parent nor any descendant references an alias the parent introduces. This matches the shape `MergeNode` path A already produces.

## Plan

### 1. Build a parent-substitution safety predicate

Add a method on `AnalysisNode` (e.g. `_can_substitute_parent_with_pks(parent)`) that returns `True` only when:

1. `parent.count` is not None, is `<= ANALYSIS_NODE_MERGE_STORE_ID_SIZE_MAX`, and is `> 0`.
2. After notional substitution, the merged `arg_q_dict` has no non-None outer keys — i.e. the child's own `_get_node_arg_q_dict` produces only `None`-keyed Qs and no descendant in the subtree introduces a Q against an alias the parent contributes.

To compute (2): collect `parent_aliases = set(parent.get_arg_q_dict().keys()) - {None}` plus any alias keys produced solely by the parent's ancestor chain in `get_annotation_kwargs`. Substitution applies when the child's own contribution and any further descendants only reference `None` keys.

### 2. Wire the substitution into `get_single_parent_arg_q_dict`

Change `AnalysisNode.get_single_parent_arg_q_dict` (`analysis/models/nodes/analysis_node.py:329-341`) so that, when the predicate from §1 holds, it returns:

```python
cached_pks = AnalysisNode.get_parent_pks(parent)
q = Q(pk__in=cached_pks)
return {None: {q: q}}
```

instead of `parent.get_arg_q_dict()`. The existing `get_parent_pks` cache (analysis_node.py:530-539) is reused as the source of truth — sibling and downstream substitutions for the same parent version share the cache hit. The `count == 0` early-exit at line 333-335 stays unchanged.

### 3. Keep `MergeNode` consistent

`MergeNode._get_arg_q_dict_from_parents_and_node` already implements the same idea per parent. Refactor it to call the same predicate / helper introduced in §1 so the alias-safety rule is uniform across single-input and merge paths. Behaviour stays identical for inputs where `count <= ANALYSIS_NODE_MERGE_STORE_ID_SIZE_MAX` already triggered path A.

### 4. Setting and ceiling

Keep the gate on `ANALYSIS_NODE_MERGE_STORE_ID_SIZE_MAX` (current default 1000) and the hard ceiling already enforced inside `get_parent_pks`.

### 5. Tests

- Unit test the safety predicate against representative parents:
  - Source nodes with `count <= threshold` and no descendant alias references → substitution applies.
  - `CohortNode` parent feeding a plain `FilterNode` (no alias references) → substitution applies.
  - `CohortNode` parent feeding a `ZygosityNode` / `AlleleFrequencyNode` child → predicate returns False (alias usage downstream).
- Integration test: build a small analysis (sample with ~200 variants → filter chain), assert `get_queryset()` SQL contains `pk IN (literal list)` rather than the parent's full filter chain. Compare row results between substitution-on and substitution-off to confirm equivalence.
- Regression test: verify the existing `MergeNode` path A produces identical results after refactor.

### 6. Profiling / measurement

Re-run the `/tmp/prof_random_page_cost` bundle (or equivalent) on an analysis where source parents have `count <= 1000`. Expect:

- Total `load_seconds` for descendant nodes drops in proportion to how much of their time was previously spent re-executing the parent's filter chain.
- EXPLAIN plans for affected children show `Bitmap Index Scan on snpdb_variant_pkey` with the literal IN list, in place of `Hash Semi Join` + FilteredRelation joins for the parent.

## Files touched (expected)

- `analysis/models/nodes/analysis_node.py` — add safety predicate; modify `get_single_parent_arg_q_dict`.
- `analysis/models/nodes/filters/merge_node.py` — refactor `_get_arg_q_dict_from_parents_and_node` to share the predicate.
- `analysis/tests/` — new tests for predicate and end-to-end substitution.
