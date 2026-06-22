# Analysis Node / Variant Query Performance Gains

Analysis of the CohortGenotype data model and Analysis node query pipeline to identify
opportunities for making large Variant queries faster.

---

## Current Architecture Summary

Variant queries are built by composing Django Q objects through a DAG of analysis nodes.
The base queryset starts with `Variant.objects.all()` (or from annotation version), then
nodes annotate and filter it. The key join is from `Variant` → `CohortGenotype` via a
`FilteredRelation` that restricts to specific `CohortGenotypeCollection` partitions.

`CohortGenotype` packs all sample genotype data into one row per variant:
- `samples_zygosity` — packed 1-char string per sample (R/H/O/.)
- `samples_allele_depth`, `samples_read_depth`, etc. — PostgreSQL arrays, 1 entry per sample
- `het_count`, `hom_count`, `ref_count`, `unk_count` — pre-aggregated zygosity counts

---

## 1. Missing GIN Indexes on Array Fields

**Problem:** Array field filters like `samples_allele_depth[i] >= min_value` and
`samples_genotype_quality[i] >= min_gq` are applied in many node types (SampleNode,
CohortNode, TrioNode via AbstractCohortBasedNode). There are no GIN or expression indexes
on these array fields in the migrations.

Without an index, every filter on `samples_allele_depth[N]` requires a sequential scan
of the relevant partition of `snpdb_cohortgenotype`.

**Opportunity:** For the most commonly-filtered fields (`samples_genotype_quality`,
`samples_allele_depth`), a GIN index on the whole array would make `ANY(array) >= value`
fast. However since the filter is on a *specific index* (`array[i]`), a functional index
would be more appropriate. Since the sample index isn't fixed, the realistic options are:

- **GIN index** on `samples_genotype_quality` and `samples_allele_depth` — would accelerate
  range queries that the PostgreSQL planner can rewrite as `value = ANY(array)` but NOT
  `array[i] >= value` at a specific position.
- **Partial index** per partition: if cohorts are small (1–3 samples), functional expression
  indexes like `CREATE INDEX ON snpdb_cohortgenotype_<pk> ((samples_genotype_quality[0]))`
  would be directly usable. This doesn't scale well for large cohorts.
- **BRIN index** on the partition table to at least prune block ranges.

**Realistic gain:** Since most partitions are already small (one VCF = one partition), the
sequential scan cost may be acceptable. The biggest win would be for the `samples_zygosity`
regex filter (see §3).

---

## 2. samples_zygosity Regex Is Anchored But Not Index-Assisted

**Problem:** Zygosity filtering uses:
```python
Q(cohortgenotype_alias__samples_zygosity__regex=r'^H..')  # example 3-sample cohort
```
And for excludes, a negative lookahead:
```python
Q(cohortgenotype_alias__samples_zygosity__regex=r'^((?!H..))')
```

PostgreSQL's `~` regex operator cannot use a B-tree index. Even though the pattern is
anchored (`^`), the planner must scan all rows in the partition.

**Opportunity — trigram (pg_trgm) GIN index:**
For longer cohorts (>3 samples), a GIN trigram index on `samples_zygosity` via
`CREATE INDEX ... USING GIN (samples_zygosity gin_trgm_ops)` would let PostgreSQL prune
rows for specific character patterns using trigram similarity. However, regex with `^` and
wildcards often can't use trigrams well.

**Better opportunity — replace regex with direct column filters:**
For single-sample nodes (the most common case), zygosity filtering is done via
`Sample.get_cohort_genotype_alias_and_field("zygosity")` which returns a `Substr` annotation
of the zygosity string. The Q filter is then `zygosity_alias__in=[zygosities]`. This is
efficient for single samples.

For multi-sample cohort/trio nodes, `get_zygosity_q()` builds a full regex across the whole
string. A potentially faster approach: instead of one regex across all positions, use per-sample
`Substr` annotations for the relevant samples only and filter on `__in`. This avoids regex
entirely. The trade-off is more annotate() calls, but each one is an O(1) substring and
the `IN` filter is index-friendly. This is already done for single-sample SampleNode — the
same pattern could be extended to trio/cohort zygosity queries.

---

## 3. VariantGeneOverlap Subquery Is Materialised In-Process

**Problem:** `GeneListNode._get_node_q()` and CompHet's `get_parent_genes()` both call:
```python
VariantGeneOverlap.objects.filter(version=v, gene__in=gene_ids_qs)
    .values_list("variant_id", flat=True)
```
This is then wrapped in `Q(pk__in=subquery_qs)`. Django sends this as a correlated subquery.
The `VariantGeneOverlap` table has `unique_together = ('version', 'variant', 'annotation_run', 'gene')`,
which gives an index on `(version, variant, annotation_run, gene)` — but the query filters on
`(version, gene)` which is not a prefix of this index.

**Opportunity:**
Add a composite index on `(version, gene)` to `VariantGeneOverlap`:
```python
class Meta:
    unique_together = ('version', 'variant', 'annotation_run', 'gene')
    indexes = [
        models.Index(fields=['version', 'gene'], name='variantgeneoverlap_version_gene_idx'),
    ]
```
This would make gene-list lookups (filtering by version + gene set) significantly faster.
The resulting `variant_id` list subquery then uses the variant FK index. This is one of the
most impactful single-index additions possible given how frequently GeneListNode is used.

---

## 4. CompHet Executes 3 Separate Queries Sequentially

**Problem:** `CompHet._get_comp_het_q_and_two_hit_genes()` in `trio_node.py:126-149`:
```python
# Query 1: variants matching mum_but_not_dad
common_genes = set(get_parent_genes(mum_but_not_dad))
# Query 2: variants matching dad_but_not_mum
                & set(get_parent_genes(dad_but_not_mum))
# Query 3: count variants per gene (two-hit genes)
two_hits = parent_genes_qs.annotate(gene_count=Count("pk")).filter(gene_count__gte=2)
```

Each `get_parent_genes(q)` call materialises a distinct gene queryset (with `.distinct()`
applied). Since `varianttranscriptannotation__gene` can have many entries per variant,
`.distinct()` is expensive on large parents.

**Opportunity — push to SQL:**
Queries 1 and 2 both touch `VariantGeneOverlap` and could be combined into a single SQL
query using `INTERSECT` or by filtering with a `GROUP BY gene HAVING COUNT(DISTINCT zygosity_pattern) = 2`
style aggregation. A single SQL approach would reduce round-trips and allow the DB planner
to optimise. However this requires bypassing the Q-based abstraction.

**Short-term opportunity — `values_list` without `distinct`:**
The current code uses:
```python
qs.values_list("varianttranscriptannotation__gene", flat=True).distinct()
```
`DISTINCT` is applied in SQL. Instead, using `VariantGeneOverlap` directly (as it already
deduplicates per variant per gene) would avoid the expensive `DISTINCT` on transcript annotation:
```python
VariantGeneOverlap.objects.filter(version=vav, variant__in=qs).values_list("gene_id", flat=True)
```
This is already what `get_overlapping_genes_q` does downstream — using it upstream in comphet
would be consistent and faster.

---

## 5. MergeNode Falls Back to In-Process ID Materialisation

**Problem:** In `merge_node.py:76-81`, when any parent has non-None annotation args (i.e., uses
a FilteredRelation like CohortGenotype), the merge falls back to:
```python
qs = parent.get_queryset(disable_cache=True)
variant_ids = qs.values_list("pk", flat=True)
or_list.append(Q(pk__in=variant_ids))
```

This executes the parent queryset eagerly and stores the list of PKs. For large parents
(millions of variants), this can cause significant memory pressure and a slow round-trip.

**Note:** There is `ANALYSIS_NODE_MERGE_STORE_ID_SIZE_MAX` to avoid this for large parents,
but that path only applies when `parent.count <= threshold` and the fallback for parents
*above* that threshold runs `parent.get_arg_q_dict(disable_cache=True)` on line 138. The
issue is specifically for parents that *have* non-None annotation args but fall in the
`<= threshold` bucket — those are materialised eagerly.

**Opportunity:** Use a proper `UNION` via `qs1.union(qs2)` at the queryset level rather than
`pk__in` with an in-memory list. Django's `.union()` produces a SQL `UNION ALL` which is
often much faster than `pk__in (SELECT ...)` for large sets because PostgreSQL can stream
rather than materialise. This needs careful handling of annotations being incompatible
across branches.

---

## 6. Sub-cohort Count Annotations Use String Functions Per Row

**Problem:** `CohortNode._get_annotation_kwargs_for_node()` (cohort_node.py:96-113) handles
sub-cohorts by building counts from substrings of `samples_zygosity`:
```python
remove_hom = Replace(sub_cohort_zygosity, Value(Zygosity.HOM_ALT), Value(''))
hom_count = Length(sub_cohort_zygosity) - Length(remove_hom)
```
This runs `SUBSTR`, `REPLACE`, and `LENGTH` in SQL for every row in the partition. For large
cohorts (many samples), `sub_cohort_zygosity = Concat(*sample_substrings)` concatenates N
`SUBSTR` calls. For a 10-sample sub-cohort this is 10+ SQL function calls per row.

**Opportunity:** Pre-compute sub-cohort CohortGenotype records. The existing two-partition
model (UNCOMMON/COMMON) already handles partitioning by variant frequency. A similar approach
of pre-computing a separate `CohortGenotypeCollection` per commonly-used sub-cohort would
replace runtime string manipulation with a direct column lookup. The downside is storage
cost and the need to rebuild on cohort changes — but this already happens for the full cohort.

**Alternative:** Use a `GENERATED ALWAYS AS` (computed column) in PostgreSQL for the
sub-cohort character position masks, rather than computing in every query. Not directly
expressible in Django ORM but achievable via migration RunSQL.

---

## 7. Contig Filtering Not Applied as Early as Possible

**Problem:** `node_queryset_filter_contigs` applies contig filtering in `get_queryset()`
*after* all annotation/filter steps are applied (analysis_node.py:586-588). For nodes
like `GeneListNode`, which already knows the set of contigs (via `_get_node_contigs()`),
the contig filter could be applied before the CohortGenotype join to reduce the join size.

**Opportunity:** Move contig filtering earlier — apply it as an initial `WHERE` clause on
the base queryset before `.annotate()`. Since `Variant.locus.contig` is a FK chain with
`Locus.position` indexed, restricting to specific contigs early can eliminate large fractions
of the variant table before any JOIN.

GeneListNode already implements `_get_node_contigs()` correctly. The issue is that
`get_queryset()` always applies it last. Exposing a method to get "pre-annotation filters"
that the base `get_queryset()` applies before annotations would allow this.

---

## 8. Debug print() Statement in Hot Path

**Problem:** `cohort_mixin.py:145`:
```python
if q_and:
    print(q_and)
```
This is a debug print in `get_cohort_and_arg_q_dict()` which is called for every cohort/trio
query. It will spam stdout and, in some environments, flush stdout synchronously, adding
latency in the hot path. Should be removed or converted to `logging.debug()`.

---

## 9. PopulationNode keep_internally_classified_pathogenic Causes Extra Query

**Problem:** `PopulationNode._get_node_arg_q_dict()` (population_node.py:185-204): when
`keep_internally_classified_pathogenic=True`, it calls:
```python
classified_variant_ids = self._get_parent_classified_variant_ids(parent)
or_q.append(Q(pk__in=classified_variant_ids))
```
`_get_parent_classified_variant_ids` materialises a full list of variant IDs (classified
as P/LP in the parent's queryset). The result is cached via `@cache_memoize` for 5 minutes.

For the first run or after cache expiry, this triggers a full execution of the parent
queryset filtered by classification — which may itself be expensive. The `pk__in` with a
Python list then embeds potentially thousands of IDs into the SQL statement.

**Opportunity:** Use a subquery expression instead of materialising IDs:
```python
classified_qs = parent.get_queryset().filter(q_classified).values('pk')
or_q.append(Q(pk__in=classified_qs))  # subquery, not list
```
This keeps everything in SQL. The `cache_memoize` approach is a workaround for what could
be a persistent cached subquery (or a NodeCache write). Using the persistent NodeCache for
"classified variants from parent" would be more durable than a 5-minute in-memory cache.

---

## 10. NodeCache Not Used by Default for Most Nodes

**Problem:** `AnalysisNode.use_cache` (analysis_node.py:492) returns `False` by default.
Only a small number of expensive nodes (VennNode, etc.) override this. Most filter nodes
recompute their queryset on every page load and count request.

However there is already an in-memory Q-object cache via `_cache_node_q` and Django's
cache framework. The issue is that for very deep analysis graphs (10+ node chains), each
page request traverses the full parent chain to reconstruct the Q dict.

**Opportunity:** Make the Q-dict cache (Redis) more aggressively used. Currently it is only
cached when `_cache_node_q = True`. For nodes that are frequently-read output nodes (i.e.,
they have `output_node=True` and many downstream users), automatically enabling `_cache_node_q`
would reduce repeated Q construction. The cache key already includes `node_version.pk` so
invalidation on edit is correct.

---

## 11. AllVariantsNode Skips the CohortGenotype Join

**Note — this is actually correct behaviour.** `AllVariantsNode` returns all variants and
has no cohort, so no `FilteredRelation` is added. Downstream nodes that do have cohort context
add their own joins. This is architecturally sound and shouldn't change.

However, when `AllVariantsNode` is used as a base and followed by a `PopulationNode` (common
in analysis templates), the query joins `Variant → VariantAnnotation` without restricting by
cohort first. On a large database this can be a very wide scan. Adding a `LIMIT` or sampling
hint to `AllVariantsNode` count queries (not grid queries) could help responsiveness.

---

## Priority Summary

| # | Issue | Effort | Impact |
|---|-------|--------|--------|
| 3 | Add `(version, gene)` index on `VariantGeneOverlap` | Very Low (1 migration) | High — every GeneListNode query |
| 8 | Remove `print(q_and)` in `cohort_mixin.py:145` | Trivial | Low-Medium (reduces I/O noise in prod) |
| 2 | Replace cohort/trio zygosity regex with per-sample Substr+IN filters | Medium | Medium-High — regex is unsupported by indexes |
| 4 | Use `VariantGeneOverlap` in CompHet instead of `varianttranscriptannotation__gene.distinct()` | Low-Medium | Medium — speeds up CompHet first-run |
| 7 | Apply contig filter before CohortGenotype JOIN in `get_queryset()` | Medium | Medium — depends on gene list sizes |
| 5 | MergeNode: use `.union()` queryset instead of `pk__in` list | High (careful testing) | Medium-High for large merges |
| 9 | PopulationNode: use subquery not list for keep_internally_classified_pathogenic | Low | Medium — removes list embedding in SQL |
| 1 | Add GIN/functional indexes on array fields | Medium | Low-Medium — partitions usually small |
| 6 | Pre-compute sub-cohort CGC records | High | Medium — only relevant for sub-cohort analyses |
| 10 | Broaden `_cache_node_q` usage for output nodes | Low-Medium | Low-Medium — reduces Q-rebuild on page loads |
