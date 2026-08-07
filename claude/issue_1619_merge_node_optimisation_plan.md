# Issue #1619 — Analysis MergeNode optimisation

## Goal

When multiple arms that filter over the **same cohort** (same `CohortGenotypeCollection`, i.e. the same
`cohortgenotype_<cgc_pk>` join) feed into a `MergeNode`, produce the merged result in a **single pass** over
the cohortgenotype partition instead of running each arm as a separate full pass, materialising each arm's PKs
into a Python list, and OR-combining them as `Variant WHERE id IN (listA) OR id IN (listB)`.

The optimisation must:
- produce **byte-for-byte identical** node counts and result PK sets to the current behaviour, and
- be **measurably a single pass** (one join to the cohortgenotype partition) for same-cohort arms, while
- leaving cross-cohort arms on the existing materialise-then-OR path (changing those would multiply LEFT JOINs
  under the OR and regress).

## Background — how the query protocol works today

### `arg_q_dict`: AND-of-keyed-filters

Every node exposes its filter as `get_arg_q_dict() -> dict[Optional[str], dict[str, Q]]`. The structure is:

```
{
  annotation_alias_or_None: { q_hash: Q, ... },   # the inner dict's Qs are AND-ed
  ...
}
```

`AnalysisNode.get_queryset()` (`analysis/models/nodes/analysis_node.py:587`) consumes it like this:

```python
a_kwargs = self.get_annotation_kwargs()        # alias -> annotation expression (incl. FilteredRelation joins)
for k, v in a_kwargs.items():
    qs = qs.annotate(**{k: v})                 # add the annotation/join
    if q_and_list := list(arg_q_dict.pop(k, {}).values()):
        qs = qs.filter(reduce(operator.and_, q_and_list))   # immediately AND-filter on that alias
# anything left under None is annotation-independent and is applied LAST:
if q_dict := arg_q_dict.pop(None, {}):
    q_list.extend(q_dict.values())
...
qs = qs.filter(reduce(operator.and_, q_list))
```

So the protocol expresses **AND of per-alias filters**, each applied right after its alias is annotated (the
annotate-then-filter cadence deliberately forces an inner query so the same partition is not joined twice).
There is **no representation for "OR across two different annotation-dependent aliases"**.

### Why same-cohort arms share a join but not a key

For a `SampleNode` over sample `S` in cohort whose collection is `cgc`, `get_annotation_kwargs()` produces:

- `cohortgenotype_<cgc_pk>` → `FilteredRelation('cohortgenotype', condition=...)` — the **shared join**
  (`CohortGenotypeCollection.get_annotation_kwargs`, `snpdb/models/models_cohort.py:521`).
- `sample_<S_pk>` → `Substr("cohortgenotype_<cgc_pk>__samples_zygosity", i, 1)` — a **per-sample** alias over
  that shared join (`Sample.get_annotation_kwargs`, `snpdb/models/models_vcf.py:470`; alias = `zygosity_alias`,
  `models_vcf.py:467`).

The zygosity filter is keyed under `sample_<S_pk>` (`SampleNode._get_node_arg_q_dict`,
`analysis/models/nodes/sources/sample_node.py:69`). Two different `SampleNode` arms over the same cohort key
their zygosity filters under **distinct** `sample_<pk>` aliases, so they never share a `q_dict` key — even
though both aliases are `Substr`s over the **same** `cohortgenotype_<cgc_pk>` `FilteredRelation`.

### What MergeNode does today

`MergeNode._get_arg_q_dict_from_parents_and_node` (`analysis/models/nodes/filters/merge_node.py:134`):
1. For each non-empty parent, take either the small-PK substitution (`get_small_parent_arg_q_dict`, only for
   arms ≤ `ANALYSIS_NODE_STORE_ID_SIZE_MAX`, default 1000) or the full `parent.get_arg_q_dict(disable_cache=True)`.
2. `_get_merged_q_dict` (`merge_node.py:107`) extracts filters common to **all** parents (applied once, AND-ed),
   then calls `_split_common_filters` on the remainder.
3. `_split_common_filters` (`merge_node.py:39`) OR-combines only the `None`-keyed (annotation-independent)
   filters. **Any parent that still has a non-`None` key** — which is every same-cohort `SampleNode`/`CohortNode`
   arm, because its zygosity filter is keyed under `sample_<pk>`/`cohortgenotype_<cgc>` — falls back to:

```python
qs = parent.get_queryset(disable_cache=True)          # a full standalone pass over the partition
variant_ids = list(qs.values_list("pk", flat=True))   # materialise PKs
or_list.append(Q(pk__in=variant_ids))                 # OR a literal pk list
```

That is the duplicate-pass cost. For *N* large same-cohort arms it is *N* full scans of the partition plus *N*
PK materialisations, OR-ed by literal id list.

### Key insight that makes the fix safe and simple

In `MergeNode.get_queryset`, `get_annotation_kwargs()` already collects annotation kwargs from **all** parents
(`AnalysisNode.get_annotation_kwargs`, `analysis_node.py:376`). Because the `FilteredRelation` is de-duplicated
via `existing_annotation_kwargs` while each arm's `sample_<pk>` alias is distinct, the merged queryset already
annotates **one** `cohortgenotype_<cgc_pk>` join plus **every** arm's `sample_<pk>` Substr over it.

Therefore a single `Q` of the form

```python
Q(sample_<A>__in=['E','O']) | Q(sample_<B>__in=['E','O'])
```

references only aliases that are already annotated by the time the `None`-bucket filters are applied (they are
applied **last**, after the whole `a_kwargs` annotate loop). Placing the OR-Q in the `None` bucket yields exactly
the desired SQL — one LEFT JOIN to the partition, two Substr expressions over it, one `WHERE (A) OR (B)`:

```sql
SELECT ... FROM variant
  LEFT JOIN cohortgenotype_partition cg
    ON (cg.collection_id IN (...) AND cg.variant_id = variant.id)
WHERE (SUBSTR(cg.samples_zygosity, idxA, 1) IN ('E','O'))
   OR (SUBSTR(cg.samples_zygosity, idxB, 1) IN ('E','O'))
```

This is approach **2 (general filter rewriting)** from the issue: annotate the shared-join aliases up front,
then apply a single OR `.filter()`. It generalises to ancestors automatically (the `sample_<pk>` /
`cohortgenotype_<cgc>` aliases propagate up unchanged through intermediate filter nodes), and it degrades to
today's materialise path for any arm whose filters reference a *different* join.

## Proposed change

### Where

All changes are in `analysis/models/nodes/filters/merge_node.py`, plus a small helper to identify the join an
alias belongs to. No change to the `get_queryset` consumer is required because the OR-Q rides in the existing
`None` bucket. (If introspection of annotation kwargs proves awkward, an optional alternative is a new dedicated
bucket key understood by `get_queryset`; see "Alternative representation" below. Prefer the `None`-bucket form —
it needs no protocol change.)

### Mechanism

Replace the per-parent fallback in `_split_common_filters` (the block at `merge_node.py:70-90`) with
join-aware grouping:

1. **Identify each arm's cohortgenotype joins.** For each `parent` still on the non-combine path, inspect
   `parent.get_annotation_kwargs()` and collect the set of cohortgenotype `FilteredRelation` alias keys it
   depends on. Helper:

   ```python
   def _filtered_relation_aliases(annotation_kwargs) -> frozenset[str]:
       return frozenset(k for k, v in annotation_kwargs.items()
                        if isinstance(v, FilteredRelation))
   ```

   For a same-cohort `SampleNode`/`CohortNode`/`ZygosityNode`/`AlleleFrequencyNode`/intermediate-filter arm this
   is `{ "cohortgenotype_<cgc_pk>" }`. (The `FilteredRelation` set alone does **not** classify an arm — step 2
   does, by resolving each *filter key* back to one of these aliases. The set is just the universe of valid
   resolution targets for this arm.)

2. **Resolve each non-`None` filter key to its owning cohortgenotype join.** This is the load-bearing helper and
   must handle every key shape these nodes actually emit — not just bare `FilteredRelation` aliases. Filter keys
   come in three forms, all of which sit over a single `cohortgenotype_<cgc>` join:
   - **the bare join alias** `cohortgenotype_<cgc>` — `SampleNode` quality/AF array filters, vcf-locus filters,
     and `CohortNode` COUNT-mode ref/het/hom counts (`CohortMixin.{ref,het,hom}_count_annotation_arg` return the
     cgc alias for non-sub-cohorts, `cohort_mixin.py:103-124`). Resolves to itself.
   - **a `sample_<pk>` Substr alias** — `SampleNode` zygosity (`Sample.zygosity_alias`, a `Substr` over
     `cohortgenotype_<cgc>__samples_zygosity`, `models_vcf.py:467-480`). Resolves to the cgc alias named in its
     `Substr` source expression.
   - **a derived cohortgenotype column/annotation key of the form `cohortgenotype_<cgc>__<suffix>`** —
     `CohortNode` COUNT-mode het-or-hom (`non_ref_call_count_column` = `cohortgenotype_<cgc>__non_ref`, an
     `F()+F()` annotation, `cohort_node.py:166-168`) and SIMPLE_ZYGOSITY-mode keys (`het_count_column` etc. =
     `cohortgenotype_<cgc>__het_count`, real columns on the join, `cohort_mixin.py:80-87`,
     `cohort_node.py:211-218`). Resolve by matching the `count_column_prefix` (`cohortgenotype_<cgc>__`) — i.e.
     split on `__` and keep the leading `cohortgenotype_<cgc>` segment **iff** it is one of this arm's
     `_filtered_relation_aliases`.

   Implement resolution as: for a key `k`, if `k` is in the arm's `_filtered_relation_aliases`, return `k`; elif
   `k` is a `sample_<pk>` alias whose `Substr` source references `cohortgenotype_<cgc>`, return that cgc alias;
   elif `k` has the form `<alias>__<rest>` where `<alias>` is in `_filtered_relation_aliases`, return `<alias>`;
   else return `None` (un-resolvable → forces materialise, safe). Capture this in a single helper
   `_arm_join_signature(parent) -> Optional[str]` returning the cohort genotype alias **iff** *every* non-`None`
   filter key resolves to the **same** single cohortgenotype join, else `None`. Resolving against the
   `FilteredRelation` set (not a bare string-prefix guess) is what keeps sub-cohort arms — whose any-sample-called
   VC alias is a *second* `FilteredRelation` and is not a `cohortgenotype_<cgc>__…` key — out of the combine path.

3. **Group the non-combine arms by join signature.** For each group whose signature is a single shared
   cohortgenotype join and which has **> 1 arm**:
   - For each arm, flatten its `arg_q_dict` into one AND-Q: `arm_q = reduce(and_, all Qs across all keys)`
     (every referenced alias — `cohortgenotype_<cgc>`, `sample_<pk>`, and any `None`-keyed real-column filters
     like gene lists — is annotated before the `None` bucket runs, so a single combined Q is valid).
   - OR the arms: `group_q = reduce(or_, arm_qs)`.
   - Append `group_q` to `or_list` (it joins the existing OR composition in `_split_common_filters`).
   - Mark these arms handled (do **not** materialise them).

4. **Everything else keeps today's behaviour.** Arms whose signature is `None` (cross-cohort, multi-join, or
   anything we cannot prove shares a single cohortgenotype join) and singleton groups fall through to the
   existing `parent.get_queryset(disable_cache=True)` → `Q(pk__in=...)` materialise path. The annotation-
   independent (`None`-keyed common) path is unchanged.

5. The combined OR-Q is returned exactly as today (`_split_common_filters` returns a single `Q`), which
   `_get_merged_q_dict` stores under `arg_q_dict[None] = {self._get_node_q_hash(): q}` (`merge_node.py:130`).
   Because it lands in the `None` bucket, `get_queryset` applies it after all aliases are annotated. No consumer
   change needed.

### Gating (must-haves from the issue)

- **Gate on shared join, not blanket OR.** Only arms that provably sit over the *same single*
  `cohortgenotype_<cgc>` `FilteredRelation` are combined. Cross-cohort / multi-join arms stay on materialise.
- **Only large arms matter.** Arms ≤ `ANALYSIS_NODE_STORE_ID_SIZE_MAX` are already collapsed to literal PK
  lists upstream (`get_small_parent_arg_q_dict`, `analysis_node.py:543`) — they never reach this path and need
  no rewrite.
- **Applies to ancestors, not just direct parents.** Detection keys on the shared annotation join (which
  propagates up through intermediate filter nodes), not on "are my parents SampleNodes", so an arm that filters
  (Zygosity, AF, quality, gene list, Tissue, etc.) *before* the merge is still detected and combined.
- **disable_cache / picklability invariants preserved.** Arms are still read via `get_arg_q_dict(disable_cache=
  True)` (the `#240 / ad35a7fb1` stale-cache protection). The combined OR-Q references only annotation aliases
  and literal values — no embedded `QuerySet` — so `arg_q_dict` stays picklable for the q-cache `cache.set`
  (the `#546` concern that motivated materialising PKs in the *old* path does not apply, because we are not
  embedding a lazy queryset).

### Alternative representation (only if `None`-bucket proves insufficient)

If any consumer path needs the OR-group annotated/inner-queried distinctly from generic `None` filters, add an
explicit bucket, e.g. a sentinel key understood by `get_queryset` that means "annotate these aliases, then apply
this one OR-Q". This is a protocol change and should be avoided unless a concrete failure (e.g. an unexpected
double-join shown in SQL inspection) forces it. The default plan is the no-protocol-change `None`-bucket form.

## Correctness concerns to handle explicitly

1. **`modifies_parents` / single-parent fast path** (`merge_node.py:16-37`) is unaffected — it triggers only
   when all parents have identical `arg_q_dict`s; combining is a no-op there.
2. **Common-filter extraction first.** `_get_merged_q_dict` extracts all-parent-common filters before
   `_split_common_filters` runs. Combined arms must still have those common filters applied — they are applied
   once under their own key (AND), and the per-arm flatten happens on the *remaining* (post-extraction) dict, so
   no filter is dropped or double-applied. Add a test with a filter common to all arms plus per-arm differences.
3. **`None`-keyed per-arm filters inside a combined arm** (e.g. a `SampleNode` gene-list restriction, keyed
   `None`). Flattening must fold these into that arm's AND-Q so the OR keeps them arm-local — otherwise a
   gene-list restriction on arm A would wrongly gate arm B. Explicit test required.
3a. **Filter-key shapes the resolver must cover.** The non-`None` filter key is **not** always the bare
    `cohortgenotype_<cgc>` `FilteredRelation` alias. `SampleNode` zygosity keys under `sample_<pk>`
    (`sample_node.py:77-79`); `CohortNode` keys under the bare cgc alias (COUNT ref/het/hom), under
    `cohortgenotype_<cgc>__non_ref` (COUNT het-or-hom `F()+F()` annotation), or under `cohortgenotype_<cgc>__<col>`
    (SIMPLE_ZYGOSITY). A resolver that only matches `isinstance(v, FilteredRelation)` aliases silently drops the
    last two shapes to materialise — correct results, but the CohortNode optimisation never fires and test #9
    fails its single-pass assertion. The step-2 resolver (bare alias → `sample_<pk>` Substr → `<cgc_alias>__…`
    prefix) is what closes this; test #9 (a)/(b)/(c) is the guard.
4. **Sub-cohort arms** (`CohortMixin`, `is_sub_cohort`) introduce an extra join — the pre-computed
   any-sample-called `VariantCollection` alias (`#1551`, `cohort_mixin.py:150`). That is a *second* join, so a
   sub-cohort arm's signature is multi-join → it stays on materialise (safe). Test that sub-cohort arms are not
   combined and still correct.
5. **Mixed arms** (one same-cohort large arm + one cross-cohort large arm + one small arm): the same-cohort
   group needs ≥ 2 members to combine; a lone same-cohort large arm with no same-cohort sibling stays on
   materialise (combining one arm buys nothing). Decide and test: combine only when group size ≥ 2.
6. **LEFT JOIN semantics under OR.** With `(subA OR subB)` over a LEFT-joined partition, rows with no
   cohortgenotype produce NULL substrings that fail both predicates — correct (those variants are not in the
   cohort). Verify via PK-set equivalence.
7. **`queryset_requires_distinct`.** Merge inherits parent distinct-ness (`analysis_node.py:404`). The combined
   single-pass query must apply `.distinct()` identically where required. Include a distinct-requiring arm
   (e.g. gene-list / multi-row join) in tests.
8. **Node counts.** `node_counts` (`analysis_node.py:928`) runs `get_node_counts_and_labels_dict` on the merged
   queryset and also asserts single-parent count sanity (not applicable to merge with >1 unique parent). Assert
   the TOTAL count matches the legacy path.

## Implementation steps

1. Add `_filtered_relation_aliases(annotation_kwargs)` and `_arm_join_signature(parent)` helpers in
   `merge_node.py` (import `FilteredRelation` from `django.db.models`).
2. Refactor `_split_common_filters` so the "non-combine parents with non-`None` keys" branch first groups arms
   by `_arm_join_signature`, builds one OR-Q per shared-single-join group of size ≥ 2, and only materialises the
   ungrouped remainder. Keep the recursion for nested common filters intact.
3. Keep `_get_merged_q_dict` storing the result under `None` (unchanged).
4. Confirm no change needed in `analysis_node.get_queryset`; if SQL inspection shows a double-join, fall back to
   the explicit-bucket representation (above) and update `get_queryset` accordingly.
5. Add an off-switch for benchmarking/regression: a settings flag (e.g. `ANALYSIS_MERGE_NODE_SINGLE_PASS`,
   default `True`) so tests can force the legacy materialise path and assert equivalence. This is also the
   safety valve for rollout.

## Testing plan

New test module: `analysis/tests/test_merge_node_single_pass.py`. Model the fixtures on
`test_explicit_pk_substitution_baseline.py` (it already builds a trio cohort with UNCOMMON + COMMON partitions,
inserts genotypes with explicit `samples_zygosity` strings, and has `_ready()` / `_child()` helpers). Run every
equivalence test **twice** — once with `ANALYSIS_MERGE_NODE_SINGLE_PASS=True` (new path) and once `False`
(legacy materialise path) — and assert identical results, exactly as the existing baseline toggles
`ANALYSIS_NODE_STORE_ID_SIZE_MAX`. Set `ANALYSIS_NODE_STORE_ID_SIZE_MAX=0` in the equivalence tests so arms are
genuinely large (forced onto the non-small path), which is where the optimisation applies.

### A. Correctness / equivalence (single-pass == legacy == union-of-arms)

For each scenario assert three things are equal: (a) merged PK set on the new path, (b) merged PK set on the
legacy path, (c) the union of each arm's standalone PK set.

1. **Two SampleNode arms, same cohort, disjoint zygosity** — proband HET arm ∪ proband HOM_ALT arm.
2. **Two SampleNode arms, same cohort, different samples** — proband-on-sampleA ∪ mother-on-sampleB (distinct
   `sample_<pk>` aliases over one shared join). This is the core case the issue targets.
3. **Three+ same-cohort arms** — verify N-way OR collapses to one pass.
4. **Overlapping arms** — arms whose PK sets intersect (union must dedupe; assert distinct).
5. **Arm with quality filters** — `min_ad`/`min_dp`/`min_gq`/`max_pl` (keyed on `cohortgenotype_<cgc>` array
   index) OR-ed with a plain zygosity arm.
6. **Arm with a gene-list restriction** (`restrict_to_qc_gene_list`, `None`-keyed) OR-ed with a non-restricted
   arm — proves `None`-keyed per-arm filters stay arm-local (concern #3).
7. **Arm with VCF locus filters** (`get_vcf_locus_filters_arg_q_dict`, keyed on the cgc alias).
8. **Intermediate filter nodes between source and merge** — `SampleNode → ZygosityNode → Merge` on one arm and
   `SampleNode → AlleleFrequencyNode → Merge` on the other, same cohort. Proves ancestor detection (concern:
   "applies to ancestors").
9. **CohortNode arms** — two `CohortNode`s over the same cohort, same-cohort → combined. Cover **every** filter
   key shape so the step-2 resolver is exercised end-to-end: (a) COUNT mode ref/het/hom (keyed on the bare cgc
   alias), (b) COUNT mode het-or-hom (keyed on `cohortgenotype_<cgc>__non_ref`, an `F()+F()` annotation), and
   (c) SIMPLE_ZYGOSITY (keyed on `cohortgenotype_<cgc>__het_count` etc.). Assert all three combine (single join,
   no `pk IN`) — a resolver that only matched bare `FilteredRelation` aliases would silently drop (b) and (c) to
   the materialise path, so this test is the guard against that regression.
10. **A filter common to all arms plus per-arm differences** — proves common-filter extraction + per-arm OR
    compose correctly (concern #2).
11. **Distinct-requiring arm** — an arm that sets `queryset_requires_distinct` (e.g. multi-gene-list) merged
    with another; assert no duplicate PKs and count equals legacy (concern #7).
12. **Variants spanning UNCOMMON + COMMON partitions** — reuse the baseline's split-partition fixture so the
    single OR pass still unions both partitions correctly.
13. **Node count equality** — for several of the above, also assert `node.get_queryset().count()` and the
    `BuiltInFilters.TOTAL` count via the counts path match between new and legacy.

### B. Cases that must STAY on the legacy materialise path (no regression)

14. **Cross-cohort arms** — two SampleNodes over two *different* cohorts (different CGCs). Assert correctness AND
    assert (via SQL inspection) the query did **not** fold them into one OR over a single join (each cohort
    keeps its own join; combining is not attempted).
15. **Sub-cohort arm** — a sub-cohort arm (extra any-sample-called VC join) is not combined with a sibling;
    result still correct (concern #4).
16. **Single same-cohort large arm + unrelated arm** — the lone same-cohort arm is not "combined with itself";
    result correct (concern #5).
17. **Small arms** — arms ≤ `ANALYSIS_NODE_STORE_ID_SIZE_MAX` still take the upstream literal-PK substitution
    (covered today by `test_explicit_pk_substitution_baseline.test_merge_node_equivalent_with_gate_on_vs_off`;
    add a mixed small+large same-cohort merge).

### C. "Is actually a single pass" (the "and is faster" requirement)

Prove the structural win deterministically rather than by wall-clock timing (timing is flaky in CI):

18. **Join-count assertion.** Render the merged queryset SQL (`str(qs.query)` or
    `library.utils.database_utils.queryset_to_sql`) and assert the cohortgenotype partition / `cohortgenotype_
    <cgc>` join appears **exactly once** for a same-cohort multi-arm merge on the new path, versus **once per
    arm** (or a `pk IN` per arm) on the legacy path. This is the direct, stable proxy for "single pass".
19. **No `pk IN (large list)` on the new path.** Assert the new-path SQL for same-cohort arms contains the
    Substr-OR predicate and does **not** contain a materialised `pk IN (...)` arm (which the legacy path emits).
20. **Optional opt-in benchmark** (marked/skipped by default, e.g. gated behind an env flag like the existing
    `VG_QUERY_PROFILE` harness): build a cohort with many variants, time the merged query new vs legacy, assert
    new ≤ legacy. Keep it out of the default suite to avoid flakiness; the join-count assertion (18) is the
    authoritative "faster" check.

### Running

```bash
python3 manage.py test --keepdb analysis.tests.test_merge_node_single_pass
# plus the existing baseline must still pass unchanged:
python3 manage.py test --keepdb analysis.tests.test_explicit_pk_substitution_baseline
```

## Risks & rollout

- **Primary risk:** an unexpected second join (Django planning the `FilteredRelation` twice) turning the
  "single pass" into a self-join. Mitigated by the SQL join-count assertions (tests 18–19); if it occurs, switch
  to the explicit-bucket representation.
- **Secondary risk:** a filter wrongly hoisted across arms (breaking arm-locality). Mitigated by tests 3, 6, 10.
- **Safety valve:** `ANALYSIS_MERGE_NODE_SINGLE_PASS=False` restores the exact current behaviour without a code
  revert, and is the toggle the equivalence tests use.
- The change is confined to `merge_node.py` (plus the settings flag); the `get_queryset` protocol is unchanged
  in the preferred design, so blast radius is small.

## Relevant code references

- `analysis/models/nodes/filters/merge_node.py` — `_split_common_filters` (`:39`), `_get_merged_q_dict` (`:107`),
  `_get_arg_q_dict_from_parents_and_node` (`:134`).
- `analysis/models/nodes/analysis_node.py` — `get_queryset` annotate-then-filter loop (`:587`), `get_arg_q_dict`
  (`:436`), `get_small_parent_arg_q_dict` (`:543`), `get_annotation_kwargs` (`:376`).
- `analysis/models/nodes/cohort_mixin.py` — `get_cohort_and_arg_q_dict` (`:150`), sub-cohort join (`:150-175`).
- `analysis/models/nodes/sources/sample_node.py` — `_get_node_arg_q_dict` (`:69`).
- `analysis/models/nodes/sources/cohort_node.py` — `_get_node_arg_q_dict` (`:172`), zygosity arg-q builders.
- `snpdb/models/models_vcf.py` — `Sample.get_annotation_kwargs` (`:470`), `get_cohort_genotype_alias_and_field`
  (`:482`), `zygosity_alias` (`:467`).
- `snpdb/models/models_cohort.py` — `CohortGenotypeCollection.cohortgenotype_alias` (`:504`),
  `get_annotation_kwargs` (`:521`), `get_zygosity_q` (`:571`).
- `analysis/tests/test_explicit_pk_substitution_baseline.py` — fixture + equivalence-toggle pattern to reuse.
- Settings: `ANALYSIS_NODE_STORE_ID_SIZE_MAX` (`variantgrid/settings/components/default_settings.py:360`),
  `ANALYSIS_NODE_CACHE_Q` (`:356`).
