# Issue #546 — Speed up small jobs by using explicit variant IDs

GitHub: https://github.com/SACGF/variantgrid/issues/546

## What the issue is asking

When a parent node holds only a small number of variants (think: 200 from a diagnostics sample), the parent's full filter chain is re-executed every time a child or sibling node composes its own queryset. PostgreSQL's planner makes poor choices when these small results are wrapped in `pk IN (subquery)` against the 40M-row variant table.

The shortcut: when the parent's `count` is small, materialise its PKs once and substitute the parent's contribution to the child's `arg_q_dict` with a literal `Q(pk__in=[…])`. The child sees a parent that filters by PK only, and Postgres plans a tight bitmap-or over `snpdb_variant_pkey`.

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

Backed by a `cache_memoize`-d helper keyed on `(parent.pk, parent.version)` (`analysis/models/nodes/analysis_node.py:530-539`), with a 15-minute TTL and a defensive ceiling that raises if `count` is None or above the same setting.

Setting: `ANALYSIS_NODE_MERGE_STORE_ID_SIZE_MAX = 1000` (`variantgrid/settings/components/default_settings.py:333`).

## What is left

Generalise the small-parent substitution from `MergeNode` to all single-input nodes that consume a parent through `get_single_parent_arg_q_dict` (`analysis/models/nodes/analysis_node.py:329-341`).

## Key findings (settle why this is safe)

These were established by tracing the code and confirmed by the baseline test below.

1. **Two independent mechanisms.** `get_queryset` (`analysis/models/nodes/analysis_node.py:569-624`) builds the query from two separate sources:
   - `get_annotation_kwargs()` produces the `FilteredRelation` joins — the `cohortgenotype_alias`, the per-sample `zygosity_alias` `Substr`, internal-count collections, etc. It walks `get_non_empty_parents()` and collects every ancestor's aliases, gated only on `has_input() and uses_parent_queryset` — **independent of `arg_q_dict`** (lines 382-399).
   - `get_arg_q_dict()` produces the WHERE clauses, keyed by which alias they apply to.

   So substituting a parent's `arg_q_dict` contribution to `{None: {pk__in}}` never removes a join. The aliases the **grid** needs for display (zygosity, allele depth, read depth, AF) and the aliases that AF/Zygosity filter nodes reference still come from the annotation path. Line 602-603 raises only if a Q references an alias missing from the annotations — which does not happen here, because substitution only ever removes an alias-keyed Q, never adds one.

2. **The cohortgenotype join is 1:1 per variant within a collection.** Splitting a parent's alias-Q (encoded as `pk__in`) from a descendant's alias-Q on the same join therefore yields the same rows as applying both on one join. The only shape where splitting could differ is a 1:many relation, and every `FilteredRelation` alias in this system is 1:1 per variant. This is why substitution is correct **even when a descendant `ZygosityNode`/`AlleleFrequencyNode` filters the same alias** — and beneficial, since the join then runs over the pk-restricted rows.

3. **Child→ancestor signalling already exists** if it is ever needed: `PopulationNode._get_kwargs_for_parent_annotation_kwargs` pushes `common_variants=False` up the chain (`analysis/models/nodes/filters/population_node.py:67-71`), and `CohortGenotypeCollection.get_annotation_kwargs` reads it to drop the common partition (`snpdb/models/models_cohort.py:511-533`). The composition is not "upward-only / descendants invisible".

4. **`MergeNode` is the safe baseline.** It already substitutes each small parent independently and isolates alias-bearing parents into full subqueries (`merge_node.py:74-84`), and a single `MergeNode` can have several CGC parents each with their own distinct `cohortgenotype_alias`. The single-input path composes alias-Qs with AND (no OR-across-siblings), so it is simpler than the merge case.

## Approach

Substitute the parent in `get_single_parent_arg_q_dict` whenever `parent.count` is in `(0, ANALYSIS_NODE_MERGE_STORE_ID_SIZE_MAX]`, reusing the same `get_parent_pks` helper and setting `MergeNode` already uses. The `count == 0` early-exit (lines 333-335) stays. Sibling and downstream substitutions for the same parent version share the `get_parent_pks` cache hit.

### 1. Wire substitution into `get_single_parent_arg_q_dict`

When the parent is ready and `0 < parent.count <= ANALYSIS_NODE_MERGE_STORE_ID_SIZE_MAX`, return:

```python
cached_pks = AnalysisNode.get_parent_pks(parent)
q = Q(pk__in=cached_pks)
return {None: {q: q}}
```

otherwise return `parent.get_arg_q_dict()` as today.

### 2. Share one helper with `MergeNode`

Factor the "small parent → `{None: {pk__in}}`" decision into a single helper on `AnalysisNode` and call it from both `get_single_parent_arg_q_dict` and `MergeNode._get_arg_q_dict_from_parents_and_node`, so the gate and the cache source are uniform. `MergeNode` behaviour stays identical for inputs that already triggered its path A.

### 3. Setting and ceiling

Keep the gate on `ANALYSIS_NODE_MERGE_STORE_ID_SIZE_MAX` (default 1000) and the hard ceiling already enforced inside `get_parent_pks`.

## Equivalence is the acceptance criterion

The optimisation must change only the SQL, never the results. The baseline test below passes on current `master`; it must keep passing with identical counts and PK sets after the change. If any assertion changes for a specific alias, that alias is the exception — add a narrow guard that falls back to `parent.get_arg_q_dict()` for that alias only, and document it.

### Baseline test (already written, passing on master)

`analysis/tests/test_explicit_pk_substitution_baseline.py` builds real `Variant` + `CohortGenotype` rows across an UNCOMMON CGC and a linked COMMON CGC (`common_collection` + `CohortGenotypeCommonFilterVersion`), makes nodes `ready` with `count` set, and runs the real `get_queryset()` chain. It pins:

- `SampleNode` source spanning both partitions (count = 4).
- `ZygosityNode` HET / HOM_ALT children filtering on the cohortgenotype alias.
- `AlleleFrequencyNode` child filtering the per-sample AF array.
- Grid genotype display columns (zygosity `Substr`, allele depth, read depth, AF) resolved through the alias join for variants in both partitions.

These are the cases the substitution must hold equal — the ZygosityNode/AlleleFrequencyNode tests prove the 1:1-join equivalence, the grid test proves the display joins survive.

### Add: rare-`PopulationNode` partition-pruning equivalence

The baseline covers the common/rare **data** path but not the `common_variants=False` **pruning** path (a rare `PopulationNode` downstream prunes the parent join to the uncommon partition via #1119). `get_parent_pks` materialises the parent standalone (`common_variants=True`, both partitions); a downstream rare filter then removes the common ones via `VariantAnnotation`, so the final result matches the pruned path when the annotation gnomAD version matches. Add a test that builds the minimal `VariantAnnotation` fixture for a couple of variants, places a rare `PopulationNode` below the small source, and asserts the result PK set is identical with the substitution gate on vs off (toggle `ANALYSIS_NODE_MERGE_STORE_ID_SIZE_MAX` via `override_settings`, e.g. 0 to disable, 1000 to enable).

## Files touched (expected)

- `analysis/models/nodes/analysis_node.py` — add the shared small-parent helper; use it in `get_single_parent_arg_q_dict`.
- `analysis/models/nodes/filters/merge_node.py` — call the shared helper from `_get_arg_q_dict_from_parents_and_node`.
- `analysis/tests/test_explicit_pk_substitution_baseline.py` — extend with the rare-`PopulationNode` pruning equivalence test (parameterised on the gate setting).

## Profiling / measurement

Re-run the `/tmp/prof_random_page_cost` bundle (or equivalent) on an analysis where source parents have `count <= 1000`. Expect:

- Descendant `load_seconds` drops in proportion to time previously spent re-executing the parent's filter chain.
- EXPLAIN plans for affected children show `Bitmap Index Scan on snpdb_variant_pkey` with the literal IN list, in place of the parent's `Hash Semi Join` + FilteredRelation joins.
