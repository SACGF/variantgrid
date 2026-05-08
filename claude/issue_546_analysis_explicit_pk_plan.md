# Issue #546 — Speed up small jobs by using explicit variant IDs

GitHub: https://github.com/SACGF/variantgrid/issues/546

## What the issue is asking

When a parent node holds only a small number of variants (think: 200 from a diagnostics sample), the parent's full filter chain is needlessly re-executed every time a child or sibling node composes its own queryset. PostgreSQL's planner makes poor choices when these small results are wrapped in `pk IN (subquery)` against the 40M-row variant table.

The proposed shortcut: when the parent's `count` is small, materialise its PKs once and substitute the parent's contribution to the child's `arg_q_dict` with a literal `Q(pk__in=[…])`. The child sees a parent that filters by PK only — no FilteredRelation joins, no annotation aliases — and Postgres plans a tight bitmap-or over `snpdb_variant_pkey`.

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

Generalise the small-parent substitution from `MergeNode` to **all** descendant nodes that consume a parent through `get_parent_arg_q_dict`. The hook point is `AnalysisNode.get_single_parent_arg_q_dict` (`analysis/models/nodes/analysis_node.py:329-341`) for single-input nodes, plus any future multi-input override (only `MergeNode` exists today).

## The blocker — annotation alias scope

Parents like `CohortNode`, `SampleNode`, `TrioNode`, `PedigreeNode` (via `CohortMixin._get_annotation_kwargs_for_node`, `analysis/models/nodes/cohort_mixin.py:52-56`) emit `arg_q_dict` keyed on a non-None alias such as `cohortgenotype_alias`. Children — and any descendant subtree under them — may produce their own Q objects keyed on that same alias (e.g. allele-frequency filters at `cohort_mixin.py:174`, VCF filter codes at `cohort_mixin.py:185-207`, zygosity filters at `cohort_mixin.py:137-155`).

If we substitute the parent's `arg_q_dict` to a single `{None: {q: Q(pk__in=…)}}` form while a descendant still references the parent-introduced alias, `AnalysisNode.get_queryset` (`analysis/models/nodes/analysis_node.py:569-624`) will hit the `arg_q_dict filters … not applied` exception at line 602, or — if the alias still survives via `get_annotation_kwargs` — produce an SQL build error because the descendant's Q references a join that the substitution intended to elide.

The substitution is therefore only safe when **no descendant Q in the resulting `arg_q_dict` references a non-None outer key that originates in the parent being substituted**.

## Plan

### 1. Build a parent-substitution safety predicate

Add a method on `AnalysisNode` (e.g. `_can_substitute_parent_with_pks(parent)`) that returns `True` only when:

1. `parent.count` is not None, is `<= ANALYSIS_NODE_MERGE_STORE_ID_SIZE_MAX`, and is `> 0`.
2. None of the non-None outer keys in `self.get_arg_q_dict()` (computed *as if* the parent had been kept normal) trace back to aliases introduced by `parent.get_annotation_kwargs()`.

To compute (2): collect `parent_aliases = set(parent.get_arg_q_dict().keys()) - {None}` plus any alias keys produced solely by `parent`'s ancestor chain in `get_annotation_kwargs`. Compare against the set of non-None keys that the child's own `_get_node_arg_q_dict` and the parent's full `arg_q_dict` would contribute. If the only non-None keys all come from `parent`'s contribution, substitution drops *all* of them — including ones the child still needs — and is unsafe.

The simplest safe rule: substitution is allowed when, after notional substitution, the merged `arg_q_dict` would have no non-None keys at all (i.e. neither the parent nor the child needs an alias that the parent contributes). This is the same shape `MergeNode` path A already produces.

### 2. Wire the substitution into `get_single_parent_arg_q_dict`

Change `AnalysisNode.get_single_parent_arg_q_dict` (`analysis/models/nodes/analysis_node.py:329-341`) so that, when the predicate from §1 holds, it returns:

```python
cached_pks = AnalysisNode.get_parent_pks(parent)
q = Q(pk__in=cached_pks)
return {None: {q: q}}
```

instead of `parent.get_arg_q_dict()`. The existing `get_parent_pks` cache (analysis_node.py:530-539) is reused as the source of truth — sibling and downstream substitutions for the same parent version share the cache hit. The `count == 0` early-exit at line 333-335 stays unchanged.

### 3. Keep `MergeNode` consistent

`MergeNode._get_arg_q_dict_from_parents_and_node` already implements the same idea per parent. Refactor it to call the same predicate / helper introduced in §1 so the alias-safety rule is uniform across single-input and merge paths. Behaviour should be unchanged for inputs where `count <= ANALYSIS_NODE_MERGE_STORE_ID_SIZE_MAX` already triggered path A.

### 4. Setting and ceiling

Keep the gate on `ANALYSIS_NODE_MERGE_STORE_ID_SIZE_MAX` (current default 1000) and the hard ceiling already enforced inside `get_parent_pks`. Consider renaming the setting to drop the `MERGE_` prefix once the optimisation applies broadly — leave the rename for a follow-up so this change is small.

### 5. Tests

- Unit test the safety predicate against representative parents:
  - Source nodes with `count <= threshold` and no descendant alias references → substitution applies.
  - `CohortNode` parent feeding a `ZygosityNode` / `AlleleFrequencyNode` child → predicate returns False (alias usage downstream).
  - `CohortNode` parent feeding a plain `FilterNode` (no alias references) → substitution applies.
- Integration test: build a small analysis (sample with ~200 variants → filter chain), assert `get_queryset()` SQL contains `pk IN (literal list)` rather than the parent's full filter chain. Compare row results between substitution-on and substitution-off to confirm equivalence.
- Regression test: verify the existing `MergeNode` path A still produces identical results after refactor.

### 6. Profiling / measurement

Re-run the `/tmp/prof_random_page_cost` bundle (or equivalent) on an analysis where source parents have `count <= 1000`. Expect:

- Total `load_seconds` for descendant nodes drops in proportion to how much of their time was previously spent re-executing the parent's filter chain.
- EXPLAIN plans for affected children show `Bitmap Index Scan on snpdb_variant_pkey` with the literal IN list, rather than `Hash Semi Join` + FilteredRelation joins for the parent.

### 7. Configuration tuning (optional)

Once correctness is established, evaluate raising `ANALYSIS_NODE_MERGE_STORE_ID_SIZE_MAX` to 5000–10000 in the env settings file for diagnostics deployments. Postgres handles literal IN lists of that size fine; the cost is one Python-side materialisation per cached parent version, which the `cache_memoize` already amortises across a 15 min window.

## Adjacent follow-up (not part of this issue)

When the parent is too large for the literal-IN substitution but downstream nodes still pay a full cohort scan repeatedly, the existing `NodeCache` / `VariantCollection` infrastructure (`analysis_node.py:503-528`, `1187-1203`) already supports SQL-level substitution via partition JOIN — opt-in per node class via `use_cache=True`. Worth a separate ticket to evaluate enabling it for `CohortNode`-shaped parents above some count. Reference: notes already captured in `claude/plans/optimise_merge_node_notes.md` §3.

## Files touched (expected)

- `analysis/models/nodes/analysis_node.py` — add safety predicate; modify `get_single_parent_arg_q_dict`.
- `analysis/models/nodes/filters/merge_node.py` — refactor `_get_arg_q_dict_from_parents_and_node` to share the predicate.
- `analysis/tests/` — new tests for predicate and end-to-end substitution.
