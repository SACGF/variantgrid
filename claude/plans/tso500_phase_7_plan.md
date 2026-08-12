# TSO 500 Phase 7 — analysis grouping node (variantgrid_private#223)

Phase 7 of [`tso500_overall_plan.md`](tso500_overall_plan.md). One entry point that gathers every VCF
sample belonging to an extraction (later a specimen), so the DNA arm's small-variant, CNV and exon-CNV
calls — and after Phase 5 the RNA arm's splice and fusion calls — land in one analysis without anyone
concatenating files outside VG.

Consumes Phase 1's `Specimen → Extraction → Sample` directly, and Phase 4 for how `Sample.extraction`
gets set. The grouping level was settled on the issue: **Extraction first, then Specimen** — Patient is
too coarse (same patient sequenced at multiple timepoints, wanted independently) and `SequencingSample`
too fine (one library/index, so it cannot unite the DNA and RNA arms).

---

## The two decisions

**Resolve the samples at query time; do not pre-build a cohort.** The node holds a reference to an
Extraction (or Specimen), finds its samples when the query runs, and ORs a per-sample subquery — the
mechanism `MergeNode` already falls back to.

**One node, not four.** `SampleNode` grows a level and the extra FKs rather than gaining
`ExtractionNode`/`SpecimenNode`/`PatientNode` siblings.

Both reverse the leaning in the earlier comments on #223, which expected `Cohort` to be the vehicle
because `Cohort.vcf` is already nullable and spans VCFs. It does — but reading what
`cohort_genotype_task` actually writes changes the answer.

## Why not a pre-built cross-VCF cohort

**It throws away the data this project exists to carry.** `_get_insert_cohort_genotype_sql`
(`snpdb/tasks/cohort_genotype_tasks.py:167-173`) hardcodes

```python
empty_json_list = "'[]'::jsonb"
columns = {
    ...
    "format": empty_json_list,
    "info":   empty_json_list,
}
```

and never populates `filters` at all — it packs only the seven `COLUMN_IS_ARRAY_EMPTY_VALUE` columns
(`snpdb/models/models_cohort.py:722`). So a custom multi-VCF cohort drops `INFO/CN`, `SEGID`,
`FORMAT/SM`, the `VG_UNDECLARED_FILTERS` values Phase 2 restores at insert, and the record-level FILTER.
Phase 2 existed to preserve those and Phase 6 exists to display them.

**A `CohortNode` over such a cohort cannot filter on FILTER either.** `has_filters` goes through
`_get_vcf()` → `cohort.get_vcf()` (`analysis/models/nodes/cohort_mixin.py:24-28`, `:242`), which is None
once `Cohort.vcf` is null. No filter UI, no filters column, and the copied rows have `filters IS NULL`
— read as PASS — anyway, so `LowUniqueAlignments` silently becomes a passing call.

**Build cost is proportional to the variant table, not the cohort.** The insert is
`FROM snpdb_variant LEFT OUTER JOIN <each partition> … WHERE (ref+het+hom+unk) > 0`, with every term
wrapped in `coalesce(…)`. Nothing is null-rejecting, so Postgres cannot convert any of those joins to
inner joins and has to scan all of `snpdb_variant` (then `SELECT DISTINCT`) to emit the ~130 rows of a
TSO 500 extraction. Per extraction. The trade isn't "pre-process once, query fast" — preparation is
O(variant table) and querying is O(tiny partition) either way.

**It goes stale against Phase 4's own design.** Adding a sample bumps `cohort.version` and rebuilds the
CGC (`Cohort.increment_version`, `snpdb/models/models_cohort.py:110`). Phase 4 is built around files and
link calls arriving in any order with a reconcile task, so the RNA arm landing a day late — or reconcile
setting `Sample.extraction` — invalidates every cohort it touches. Query-time resolution has no
staleness state to keep.

**And it still wouldn't give sapath#301 what it needs.** `AbstractCohortBasedNode._get_q_and_list`
(`analysis/models/nodes/sources/cohort_node.py:52-72`) applies one `min_ad`/`min_dp` across all packed
indices with an ANY/ALL group operation; there is no per-sample threshold. Per-caller cutoffs would
still need N× `SampleNode` → `MergeNode`.

**Where a cohort would win, for the record.** Genuine cross-sample zygosity questions — "called in both
arms", "called in ≥2 of these samples" — are one indexed regex over a packed row
(`CohortGenotypeCollection.get_zygosity_q`), and only work within a single CGC. The group node also
gives up `get_cached_label_count_for_cohort`, so its node counts are live queries. Both are acceptable
at 3-4 samples over ~100-row partitions, and either can be added *behind* this node later without
changing its identity: a `VariantCollection` (`snpdb/models/models_variant.py:930`) is the vehicle if
live counts ever get slow — it is a partitioned PK set with its own `get_arg_q_dict` and archive
handling, so it caches the result without touching genotype data the way a cohort copy would — and a
cohort if a genuine cross-sample zygosity requirement appears.

## How the query is built

Per sample, build the arg_q_dict `SampleNode` builds today, render it to SQL, and OR the results on
`pk` — exactly `MergeNode._split_common_filters` (`analysis/models/nodes/filters/merge_node.py:91-94`):

```python
pk_qs = parent.get_queryset(disable_cache=True).values_list("pk", flat=True)
sql, params = pk_qs.query.sql_with_params()
or_list.append(Q(pk__in=RawSQL(sql, params)))
```

Factor that out of `MergeNode` into a shared helper so there is one implementation. Read the comment
block above it first: it documents why `RawSQL` rather than `list(qs.values_list("pk"))` — picklability
for the q-cache `cache.set`, and not materialising an unbounded PK list into Python.

**Not** an OR across the `cohortgenotype_<pk>` aliases. That *would* compile — `get_queryset`
(`analysis/models/nodes/analysis_node.py:626-646`) applies every annotation before the `None`-keyed
filters — but it is an OR over LEFT JOINs, which reproduces the same unprunable full scan as the cohort
build above. The per-alias filter loop exists precisely to keep each join in its own subquery.

Two consequences worth stating:

- **No `distinct()` needed.** `pk__in` subqueries don't fan out rows the way joins do, so
  `queryset_requires_distinct` stays False.
- **One sample degenerates to today's query.** A single-sample group should short-circuit to the
  existing alias path (`Sample.get_cohort_genotype_alias_and_field`, `snpdb/models/models_vcf.py:473`,
  and the `Substr` zygosity annotation at `:461`). An extraction with one arm and one caller then
  produces byte-for-byte the query `SampleNode` produces now.

VCF FILTER handling comes out well. `NodeVCFFilter` stores filter *ids* on the node and translates them
per VCF at query time (`get_filter_codes(node, vcf)`, `analysis/models/nodes/analysis_node.py:1352-1361`),
so one node-level selection ("PASS", "LowUniqueAlignments") resolves independently into each VCF's code
set. `get_vcf_locus_filters_arg_q_dict` (`cohort_mixin.py:~200`) just needs to loop the samples' VCFs
instead of calling `self._get_vcf()` once.

**Cache key trap.** `CohortMixin._get_cache_key` (`cohort_mixin.py:30-41`) folds in a single CGC pk so
the node's cached arg_q_dict invalidates when a VCF is reloaded. At group level it must fold in every
sample's CGC pk *and* cohort version, or a reloaded VCF leaves a stale cached query behind.

## One node with a level, not four node classes

The only thing that differs between levels is "object → set of samples". Zygosity/AD/DP/GQ, the gene-list
restriction, the locus filters, the Q construction, `_get_cohorts_and_sample_visibility_for_node`, the
grid columns — all identical once the sample set exists. Four classes would share every one of those
through a mixin and differ in one resolver method.

The decisive argument is templates. `AnalysisTemplateRun.populate_arguments`
(`analysis/models/models_analysis.py:610`) clones a fixed DAG and sets fields *by name*, with
`AnalysisVariable` keyed `(node, field)` and its `class_name` derived from the FK's related model. A
separate `ExtractionNode` means every template that starts from a `SampleNode` needs a duplicate built
around it — the whole downstream DAG, per level, per deployment. As a field on the same node, the
template author exposes `extraction` as the variable and nothing else in the template changes, which
makes adding `extraction`/`specimen` to `analysis_templates_tag`'s `single_model_args`
(`analysis/templatetags/related_analyses_tags.py:134`) the two-line change the overall plan predicted,
since `hidden_inputs` is already keyed on field name.

A fair amount of code also keys on `SampleNode` by class and keeps working:
`analysis/tasks/analysis_grid_export_tasks.py:183` (`SampleNode.objects.get_subclass(analysis=analysis,
output_node=True)`), `analysis/related_analyses.py:23`, `analysis/forms/forms.py:84`.

And Patient becomes cheap rather than a whole node. The overall plan ruled it out partly because
`Cohort.genome_build` forces a single build — an objection that dies with the cohort, since query-time
resolution just restricts to the analysis build. Patient is then one more choice that a deployment can
leave switched off.

### What the class split would have bought, and how to get it anyway

Discoverability. The add-node menu is built from `get_node_types_hash()`
(`analysis/models/nodes/node_types.py:19`), keyed by class, so one class is one menu entry — a user
hunting for "Specimen" would have to add a Sample node and change a dropdown.

Fix it in the menu rather than the models: let a node class declare extra entries carrying initial field
values, so `node_create` (`analysis/views/views_json.py:115-127`, currently
`NODE_TYPES_HASH[node_type]` → `objects.create(...)`) can split `SampleNode:EXTRACTION` into class plus
defaults. Four menu entries, four labels and icons, one table.

### Shape

- **`source_level` field** (SAMPLE / EXTRACTION / SPECIMEN / PATIENT) alongside nullable `extraction`,
  `specimen` and `patient` FKs beside the existing `sample`. Explicit rather than inferred from which FK
  is set, because `node_create` stamps it from the menu entry and the editor must choose a widget before
  any value exists. Defaults to SAMPLE, so every existing row is unchanged.
- **Sample-level only:** `sample_gene_list`, `restrict_to_qc_gene_list`, `GeneCoverageMixin`. Hide them
  at group levels rather than inventing a meaning over N samples.
- **`_get_proband_sample_for_node`** returns the DNA arm's sample at extraction level. That is what makes
  downstream `AncestorSampleMixin` nodes (Zygosity, MOI, GeneList, AlleleFrequency) auto-populate
  sensibly instead of giving up on an ambiguous ancestor set.
- **`_stats_cache.py:99`** maps node class → handler, so the handler branches on level. A group-level
  node returns `None` rather than silently reusing a single-cohort count.
- **`handle_sample_pre_delete`** (`analysis/signals/source_data_invalidation.py:55`) needs a second Q —
  `Q(samplenode__extraction=instance.extraction)` — because deleting a sample changes a group node's
  result without touching its FK. That work exists under either design.
- **Per-sample thresholds are a child table**, keyed `(node, sample)`, carrying the same
  `min_ad`/`min_dp`/`min_gq`/`max_pl` set the node already holds
  (`analysis/models/nodes/sources/sample_node.py:24-27`). `NodeAlleleFrequencyFilter` →
  `NodeAlleleFrequencyRange` (`analysis/models/nodes/analysis_node.py:1364`) is the precedent for a
  per-node child table feeding `_get_node_arg_q_dict`. The node's own fields stay as the value applied
  to any sample with no row — so a sample that `reconcile_pending_extractions` attaches after the node
  was configured gets the node default rather than nothing, and the editor shows a row only where the
  user has overridden one.
- **`SampleNodeForm`** is where the complexity actually lands. It already deletes fields conditionally
  for `has_genotype` and `lock_input_sources` (`analysis/forms/forms_nodes.py:783-793`); it gains
  level-conditional widgets and the per-sample threshold rows.

## Configuring the node

### Autocompletes — mostly already built

Phase 1 shipped the forwarding chain. `SpecimenAutocompleteView` and `ExtractionAutocompleteView`
(`patients/views_autocomplete.py:26`, `:40`) already narrow on whichever of patient/specimen/extraction
the form has set. The gap is the far end: `SampleAutocompleteView`
(`snpdb/views/views_autocomplete.py:114`) forwards genome build and `exclude_archived` but not
extraction — one clause completes the chain in both directions. On the form side, the `forward=` wiring
the `sample_gene_list` widget already demonstrates (`analysis/forms/forms_nodes.py:781`), plus
`genome_build_fields` extended past `["sample"]` so the specimen and extraction lists are
build-restricted too.

### A preview endpoint, so the node's reach is visible before it runs

The counts already exist: `CohortGenotypeStats` (`snpdb/models/models_cohort_stats.py:24`) holds
per-sample `variant_count` keyed on `(cohort_genotype_collection, sample, passing_filter)`. So

`GET …/extraction/<pk>/samples` → `[{sample, vcf, genome_build, variant_count, has_genotype}]` + totals

is a join over stats rows with no variant-table access. Keyed on the *object* pk rather than the node,
since it has to answer before the node is saved — which also makes it reusable by Phase 3's specimen
page and Phase 4's reconcile UI.

Two things decide whether it helps or misleads:

- **Apply exactly the filters the node will** — `Sample.filter_for_user`, the analysis genome build,
  archived VCFs. A preview that counts samples the node then drops is worse than no preview.
- **Report exclusions rather than omitting them.** "4 samples, 3 in GRCh38 — 1 not included" is the
  message. Silent narrowing is the failure mode that bites here, precisely because the user didn't
  hand-pick the samples.

An analysis has exactly one genome build, so restricting to it is a fact rather than a policy choice —
which is what makes Patient level cheap even where timepoints span builds. The warning is the whole of
the work, and it belongs on the node as well as in the preview, since the exclusion happens at query
time: a build the node silently drops is the same failure as a VCF that was archived, and
`NodeStatus`/the node's warning text is where the user will actually see it.

Where stats haven't been calculated yet, return null and let the UI say so.

## Seeing what you got

### Which VCF a row came from

The PK-driven filter carries no provenance, but the grid doesn't display from the filter.
`VariantGrid._get_grid_only_annotation_kwargs` (`analysis/grids.py:165-176`) passes **all** the node's
cohorts to `get_variantgrid_zygosity_annotation_kwargs`, which loops them
(`snpdb/grid_columns/grid_sample_columns.py:44-64`), adds a FilteredRelation per CGC and coalesces each
packed column to a per-cohort placeholder:

```python
empty_data = empty_value * cohort.sample_count       # e.g. "..."
annotation_kwargs[packed_column] = Coalesce(f"{alias}__{column}", Value(empty_data), ...)
```

So every row already arrives carrying, per VCF, either real packed genotype data or the missing-value
placeholder. "Which VCF is this from" is `zygosity != '.'` per cohort — a client-side formatter over
data the grid already fetches, no extra join, no change to the filtering path. It is also why
`format`/`info` survive per VCF, which the cohort-copy route would have destroyed.

Render it as one derived "Source" column of compact badges rather than making users read it off N sets
of packed columns. Column count grows with sample count, which is nothing at 3-4 samples per extraction
and worth watching at specimen level with many timepoints.

### The node badge

`SampleNode.get_rendering_args()` already returns patient JSON and `sampleNodeUpdateState`
(`variantgrid/static_files/default_static/js/samplenode.js:82-103`) draws the male/female/deceased
pedigree symbol from it. A small VCF glyph with a ×N badge is the same mechanism with another key in the
dict; `appearance_dirty` / `appearance_version` (`analysis/models/nodes/analysis_node.py:1061`) already
force the re-render when configuration changes. Hover is a `title` listing the VCF names; a richer
popover is more work.

- `get_rendering_args()` runs for every node on every analysis render, so resolving samples → VCFs there
  wants the treatment `get_sample_ids` gets — `@cache_memoize(DAY_SECS, args_rewrite=lambda s: (s.pk,
  s.version))` (`analysis/models/nodes/analysis_node.py:219`).
- The badge and the preview endpoint must come from **one** server-side helper on the node. If the ×3 on
  the canvas and the "3 VCFs" in the editor can disagree, the badge stops being worth glancing at.

The badge answers a question the endpoint can't: six months later, why does this specimen node return
fewer rows than it used to — because a VCF was archived.

## Order of work

1. Extraction level: `source_level` + FKs, the resolver, the shared subquery-OR helper extracted from
   `MergeNode`, the single-sample short-circuit, per-VCF locus filters, cache key.
2. Form and autocomplete forwarding (the one `SampleAutocompleteView` clause), then the preview endpoint.
3. The grid Source column.
4. The node badge, sharing the endpoint's helper.
5. Specimen level — the same node with a wider sample query. Patient after that, if wanted.

Per-sample threshold rows land with (1) or (2) depending on how the form shakes out; they are what
sapath#301 actually asked for, so they are not optional at extraction level.

**Done when** one extraction node on the TSO 500 test run returns the DNA arm's small-variant, CNV and
exon-CNV calls in a single grid, each row showing which VCF it came from and carrying that VCF's
`INFO`/`FORMAT` values, with different `min_ad` per caller, and the node showing ×3 on the canvas.

## Settled, with what was decided against

| Question | Decision |
|---|---|
| Per-sample threshold UI | **Repeated rows in the editor**, over a per-VCF-source default. VG's usual shape is a global default you may then override per object, so the rows are needed either way — and the defaults want to key on panel as well as caller, which is its own piece of work |
| Patient level and multiple genome builds | **Restrict to the analysis build and warn** which VCFs that excluded. An analysis is single-build by definition, so this is a message rather than a decision |
| `Sample.extraction` is a single FK (`snpdb/models/models_vcf.py:340` — "TODO: A sample may have >1 extractions (eg tumor/normal subtraction)") | **Leave it.** If it becomes real the query-time resolver absorbs it; a materialised cohort would not have |
| Cached label counts at group level | **Live counts, no cache.** `VariantCollection` is the next step if a specimen ever spans enough timepoints to hurt |

Deployment-wide threshold defaults per VCF source — and per panel, which is the part that makes it
bigger than a settings dict — are #1717 rather than part of this phase. The per-node values this phase
builds stay the thing that runs either way; #1717 only changes what a new node starts from.
