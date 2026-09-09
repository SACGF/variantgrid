# analysis — research notes

Verified against 7c4408c62 on 2026-09-06

The analysis app is the interactive filter: an `analysis/models/models_analysis.py:Analysis` is a DAG of
`analysis/models/nodes/analysis_node.py:AnalysisNode` subclasses over one genome build and one pinned AnnotationVersion,
where source nodes (sample, cohort, trio, duo, quad, pedigree, all variants) start the graph and filter nodes narrow it,
each contributing a Django Q rather than a stored result set. Around that core sit the version/lease machinery that
recounts nodes in celery after every edit, the node grid and its exports, templates that are cloned and parameterised per
sample or cohort (and auto-launched on VCF import), variant tags, and the trio karyomapping side-feature. This document is
the story behind `analysis/CLAUDE.md`: how an edit becomes a count, why a node is a Q and not a table, and what has gone
wrong before. Model fields, URLs, commands, tasks, signals and settings are in the generated maps
([models](../maps/models.md#analysis), [urls](../maps/urls.md#analysis), [commands](../maps/commands.md),
[tasks](../maps/tasks.md#analysis), [signals](../maps/signals.md), [settings](../maps/settings.md));
`claude/domain.md` has the vocabulary.

## Flows

### A node edit, end to end

An editor POST lands in `analysis/views/nodes/node_view.py:NodeView.form_valid`, which sets `queryset_dirty`, saves and
calls `analysis/models/nodes/node_utils.py:update_analysis`; a drag or a connection change goes through
`analysis/views/views_json.py:NodeUpdate` and does the same (a pure move saves without dirtying).
`analysis/models/nodes/analysis_node.py:AnalysisNode.save` takes a `select_for_update` on the Analysis row and hands off
to `analysis/models/nodes/analysis_node.py:AnalysisNode._save`, which recomputes validity and shadow colour, and - only
when `parents_changed` or `queryset_dirty` - calls `analysis/models/nodes/analysis_node.py:AnalysisNode.bump_version`
(version+1, DIRTY, count/errors/cloned_from cleared) and then recursively saves every child with `queryset_dirty=True`,
so a whole subtree bumps in one transaction. The save ends by `get_or_create`ing the
`analysis/models/nodes/analysis_node.py:NodeVersion` for the new number; everything scoped to a version (the Redis Q
cache key, `analysis/models/nodes/analysis_node.py:NodeCache`, `analysis/models/nodes/analysis_node.py:NodeTask`, column
summaries) hangs off that row and cascades away when
`analysis/tasks/node_update_tasks.py:delete_analysis_old_node_versions` deletes every NodeVersion that is not its node's
latest.

`update_analysis` fires that cleanup and `analysis/tasks/analysis_update_tasks.py:create_and_launch_analysis_tasks` on
`scheduling_single_worker`. The scheduler is state-driven (#346): it does not build a chain of tasks mirroring the
graph, it calls `analysis/tasks/analysis_update_tasks.py:lease_ready_nodes`, which locks the analysis' unsettled nodes
(`skip_locked`), loads the graph once via `analysis/tasks/analysis_update_tasks.py:_load_graph`, and for each DIRTY
node whose parents are all settled (`analysis/tasks/analysis_update_tasks.py:_node_ready_to_lease` - an error parent
counts as settled, so the child runs and fails fast with ERROR_WITH_PARENT) writes a NodeTask lease (`leased_by`,
`lease_expires`, `attempt_count`, `run_after`) and flips the node to QUEUED. At most
`ANALYSIS_NODE_DISPATCH_MAX_NODES_PER_ANALYSIS` go out per dispatch; each leased node becomes the signature from
`analysis/tasks/analysis_update_tasks.py:_node_launch_signature`, an `update_node_task` on `analysis_workers`, chained
behind a `node_cache_task` only when the node builds its own cache.

`analysis/tasks/node_update_tasks.py:update_node_task` loads the subclass at the leased version, then
`analysis/models/nodes/analysis_node.py:AnalysisNode.claim_for_load` does a conditional UPDATE from a claimable status to
LOADING and restarts the lease window (a task that sat in a backlog should not spend its lease queued). After a reclaim two
tasks can exist for one node/version; the UPDATE is the arbiter and the loser exits without touching anything. The winner
runs `analysis/models/nodes/analysis_node.py:AnalysisNode.load` (below), maps exceptions to statuses
(`NodeOutOfDateException` means the user edited meanwhile - exit quietly, a newer version is coming; `OperationalError`
backs off through `analysis/tasks/node_update_tasks.py:_backoff_node`, which sets DIRTY with a future `run_after`;
`MemoryError` under the worker's RLIMIT_AS perma-fails and re-raises as `NodeOutOfMemoryException` so Rollbar sees it),
clears the lease in a `finally`, and calls `analysis/tasks/node_update_tasks.py:_trigger_rescheduling` - the dispatcher
immediately and again in 3 s, covering the race where a concurrent dispatcher read statuses just before this commit. A
node finishing is what leases its children. `analysis/tasks/analysis_update_tasks.py:dispatch_analysis_backlog` on beat
is the catch-all: it tops the in-flight count up to `ANALYSIS_NODE_DISPATCH_BACKLOG_IN_FLIGHT_TARGET` from analyses with
un-leased loading nodes, and it is where dead workers are noticed. Both dispatchers honour the `JobsControl` pause.

### Counts, stored pks and the grid

`AnalysisNode.load` runs the subclass `_load` (returning fields to persist), then
`analysis/models/nodes/analysis_node.py:AnalysisNode.node_counts`, then writes everything through
`analysis/models/nodes/analysis_node.py:AnalysisNode.update`, a conditional UPDATE on (pk, version). `node_counts` first
records `live_data_sources` on the NodeVersion, then asks
`analysis/models/nodes/analysis_node.py:AnalysisNode._get_cached_label_count` for each configured label (TOTAL plus the
analysis' `analysis/models/models_analysis.py:AnalysisNodeCountConfiguration` records, which include one label per tag
when `node_count_auto_add_tags` is on): a clone reuses `cloned_from`'s counts, a pass-through node reuses its single
parent's, all-zero parents short-circuit to zero, and a SampleNode with no count-affecting filters reads the precomputed
`CohortGenotype*Stats` row via `analysis/models/nodes/stats_cache.py:get_cached_label_count_for_cohort`. Whatever is
left is one aggregate query in `analysis/models/nodes/node_counts.py:get_node_counts_and_labels_dict`. For a node at or
under `ANALYSIS_NODE_STORE_ID_SIZE_MAX` (1000) `_get_variant_ids_to_store` also pulls the exact pk list, and the list
length becomes the total - the count and the pks came from the same load. Counts, `tag_counts` from `_get_load_data` and
the pks land on the NodeVersion as `load_data` / `variant_ids` in one `.update()` that stamps `modified`, because the
client's `analysis/views/views_json.py:nodes_status` poll uses `modified` to know a recount has landed. `analysis/models/nodes/analysis_node.py:AnalysisNode._raise_or_warn_count_mismatch`
then checks for a label count above the total or a single-parent node bigger than its parent: an error for a
deterministic node, a warning for one whose `live_data_sources` say its tables move under it (#235).

The grid reads the same queryset. `analysis/views/views_node.py:node_load` redirects to errors, the grid or an async-wait
page by status; `analysis/views/views_grid.py:NodeGridConfig` and `analysis/views/views_grid.py:NodeGridHandler` are
`cache_page`d for a week under the node version in the URL, and the handler holds a per-user, per-URL cache lock and a
`major_operation` slot so a double-click cannot run a minutes-long query twice.
`analysis/views/views_grid.py:NodeGridHandler._get_redirect` sends a pass-through node to its parent's URL
(`analysis/models/nodes/analysis_node.py:AnalysisNode.get_grid_node_id_and_version`) so the two share one cache entry.
`analysis/grids.py:VariantGrid` builds columns from the analysis' `CustomColumnsCollection` (every viewer sees the same
grid, permissions checked as the analysis owner), adds per-node extra columns (`_get_node_extra_columns` - cohort counts,
sample genotype cells, VCF FILTER), and uses `analysis/grids.py:VariantGrid.known_count` so the datatable never re-counts
what the load already counted; sorting is disabled above `ANALYSIS_GRID_SORT_MAX_ROWS`. Exports
(`analysis/views/views_grid.py:node_grid_export`) create a `CachedGeneratedFile` keyed on the hash from
`analysis/tasks/analysis_grid_export_tasks.py:get_node_grid_downloadable_file_params_hash` (node, version, user, filters,
transcript collection) and run `analysis/tasks/analysis_grid_export_tasks.py:export_node_to_downloadable_file` in celery
(#1257), which retries itself while the output node is still loading rather than blocking a worker;
`analysis/grids.py:ExportVariantGrid.iter_export_rows` walks pks a contig at a time so the annotation joins never see a
full-table sort.

### Materialised nodes: Venn and Intersection

Two nodes cannot be a Q and so write a `VariantCollection` instead. `analysis/models/nodes/filters/venn_node.py:VennNode`
keeps its own `analysis/models/nodes/filters/venn_node.py:VennNodeCache` keyed on the two parent NodeVersions and a
region (A-only, intersection, B-only); `analysis/models/nodes/filters/venn_node.py:VennNode.get_cache_task_args_set`
`get_or_create`s one per region the set operation needs and returns a
`analysis/models/nodes/filters/venn_node.py:venn_cache_count` task for any region not yet SUCCESS, which
`_node_launch_signature` chains ahead of the node's own update. `venn_cache_count` pulls both sides' pks into Python and
does set arithmetic there, skipping any side an empty parent already decides - SQL EXCEPT/INTERSECT across the partition
tables used to die with "too many range tables". It truncates partial records first and is idempotent, so a re-lease or a
sibling Venn over the same parents reuses an in-flight build. `analysis/models/nodes/filters/venn_node.py:VennNode._get_node_q`
raises: the node's Q always comes from `_get_node_cache_arg_q_dict` as `variantcollectionrecord__variant_collection__in`.

`analysis/models/nodes/filters/intersection_node.py:IntersectionNode` is the last user of the generic
`analysis/models/nodes/analysis_node.py:NodeCache`: `use_cache` is true for a selected BED collection or a long pasted
variant list, `analysis/models/nodes/analysis_node.py:AnalysisNode.get_cache_task_args_set` creates the NodeCache (or
finds the parent's, for a pass-through node) and the chained `analysis/tasks/node_update_tasks.py:node_cache_task` calls
`analysis/models/nodes/filters/intersection_node.py:IntersectionNode.write_cache`, which streams the parent queryset as a
VCF into an `intersectBed` pipe that loads the collection itself. A node whose cache another node is building is held
back by `analysis/tasks/analysis_update_tasks.py:_node_cache_ready` while the collection is PROCESSING; a CREATED
collection with no builder does not block, the node just runs live.

### Templates and auto-analyses on import

A template is an Analysis with `template_type=TEMPLATE` plus `analysis/models/models_analysis.py:AnalysisVariable` rows
naming the node fields a run must supply; the editor only offers the variable widget on source-node fields of a template
(`analysis/views/nodes/node_view.py:NodeView.get_form`). `analysis/models/models_analysis.py:AnalysisTemplate.new_version`
refuses without a source-field variable, clones the analysis via `analysis/models/models_analysis.py:Analysis.clone`
(toposorted `save_clone` of every node with `cloned_from` pointing at the original NodeVersion, edges and variables
re-pointed, the node-count configuration copied so "all counts off" survives) into an invisible SNAPSHOT, and creates an
`analysis/models/models_analysis.py:AnalysisTemplateVersion` with `active=False` - a draft only writers can run until
`analysis/models/models_analysis.py:AnalysisTemplateVersion.activate` makes it the one everyone runs (#1496).
`analysis/models/models_analysis.py:AnalysisTemplateVersion.filter_for_user` returns exactly the active versions the
user can see plus the drafts they can write.

`analysis/analysis_templates.py:run_analysis_template` is the run: `analysis/models/models_analysis.py:AnalysisTemplateRun.create`
clones the snapshot into a visible analysis on the latest validated AnnotationVersion for the build;
`analysis/models/models_analysis.py:AnalysisTemplateRun.populate_arguments` resolves each variable (Guardian
`check_can_view` as the run's user, type-checked against `class_name`) into an
`analysis/models/models_analysis.py:AnalysisTemplateRunArgument`, recording an error rather than raising; then
`analysis/analysis_templates.py:populate_analysis_from_template_run` names the analysis from `analysis_name_template`
(`%(input)s` is the first of pedigree/trio/quad/duo/cohort/sample), sets the fields node by node in toposort order with
`queryset_dirty=True` so `_save` cascades sample changes downstream, hides any node that
`hide_node_and_descendants_upon_template_configuration_error` says should vanish along with its descendants
(`analysis/views/views_json.py:node_reveal_hidden` brings them back), and finishes with
`analysis/models/nodes/node_utils.py:reload_analysis_nodes` - a bulk bump of every node's version and status in two
UPDATEs plus a `bulk_create` of NodeVersions, then `update_analysis`.

Auto-analyses hang off `upload`'s `vcf_import_success_signal`: `analysis/signals/signal_handlers.py:handle_vcf_import_success`
chains `analysis/tasks/auto_analysis_tasks.py:auto_run_analyses_for_vcf` (one run per sample per matching
`analysis/models/models_analysis.py:AutoLaunchAnalysisTemplate`, matched on enrichment kit and a sample-name regex by
`analysis/analysis_templates.py:get_auto_launch_analysis_template_matches`, skipped when the sample already has a related
analysis) and `analysis/tasks/auto_analysis_tasks.py:reload_auto_analyses_for_vcf`, which re-bumps the sample-tab and
cohort-export analyses of a re-imported VCF. A template that `requires_sample_gene_list` is skipped until the QC gene list
arrives, and `analysis/signals/signal_handlers.py:handle_active_sample_gene_list_created` retries then. The sample tab
and the cohort VCF export use `analysis/analysis_templates.py:get_sample_analysis` /
`analysis/analysis_templates.py:get_cohort_analysis`, one hidden run per object per template version recorded in
`SampleAnalysisTemplateRun` / `CohortAnalysisTemplateRun`, named by `ANALYSIS_TEMPLATES_AUTO_SAMPLE` and
`ANALYSIS_TEMPLATES_AUTO_COHORT_EXPORT`.

### Tags

A `analysis/models/models_variant_tag.py:VariantTag` is a user's tag on a Variant, optionally inside an analysis (and
node/node version); permissions delegate to the analysis when there is one. Creating or deleting one runs
`analysis/signals/signal_handlers.py:variant_tag_create` / `variant_tag_delete`: synchronously it dirties the visible
`analysis/models/nodes/filters/tag_node.py:TagNode`s that read that tag (`analysis/tasks/variant_tag_tasks.py:analysis_tag_nodes_set_dirty`,
so the client's next `analysis_node_versions` poll already sees the bump) and adds or removes the tag's node-count label
(`analysis/tasks/variant_tag_tasks.py:update_analysis_tag_node_count_config`, #21); on commit it queues
`analysis/tasks/variant_tag_tasks.py:variant_tag_created_task`, which dirties the hidden "Tagged Variants" node, recounts
the tag labels in place and links the tag to an Allele and a liftover pipeline. Tagging does not bump the other nodes'
versions - the tag count on a node badge is recomputed by
`analysis/models/nodes/node_utils.py:update_analysis_tag_node_counts` against the NodeVersion each READY node is already
on, looking the tagged variants up once per analysis (`analysis/models/nodes/node_counts.py:get_tagged_variant_ids_by_label`)
and merging with `jsonb_set` + `||` so a concurrent load cannot lose its own labels.
`analysis/models/nodes/filters/tag_node.py:TagNode.tagged_variants_q` is the one place local versus global scope is
decided: this analysis' tags by Variant pk (avoiding the Allele race right after tagging), other analyses' by Allele
through `analysis/models/models_variant_tag.py:VariantTag.variants_for_build_q`, with `tagged_within_days` anchored to
the NodeVersion's creation so an old analysis reproduces what it showed. A global node is a snapshot and says so in
`get_warnings`. The editor's tag picker is `analysis/models/nodes/filters/tag_node.py:TagNode.get_tag_counts`, one grouped
count over the taggings in scope - it is a hint beside each pill rather than the node's own count, so it is worked out
when the editor renders and the load never touches it (#1820, #1838).

### Source nodes

Every genotype-aware node mixes in `analysis/models/nodes/cohort_mixin.py:CohortMixin`: `_get_cohort` names the cohort,
`cohort_genotype_collection` finds its packed genotype partition (archived or missing data becomes a configuration error,
not a 500), `_get_annotation_kwargs_for_node` registers the `cohortgenotype_<id>` alias and
`analysis/models/nodes/cohort_mixin.py:CohortMixin.get_cohort_and_arg_q_dict` puts the zygosity, quality, VCF FILTER and
allele-frequency Qs under that alias so `annotate_and_filter_queryset` joins the partition once. Its `_get_cache_key`
folds in the collection pk (and the sub-cohort any-sample-called VariantCollection of #1551, which turns a regex
exclusion into a hash join) so a reloaded VCF invalidates the cached Q.
`analysis/models/nodes/sources/sample_node.py:SampleNode` is one sample or, at extraction/specimen/patient level
(`patients/models_enums.py:SampleSourceLevel`, offered per `ANALYSIS_SAMPLE_NODE_LEVELS`), every sample of a grouping object: a single sample degenerates to the alias path, a group
ORs one `pk IN (subquery)` per sample via `analysis/models/nodes/cohort_mixin.py:get_sample_pk_in_q`,
each annotated with only its own VCF's join. Per-sample overrides live in
`analysis/models/nodes/sources/sample_node.py:SampleNodeSampleFilter`. Downstream nodes that need "the sample" get it
through `analysis/models/nodes/cohort_mixin.py:AncestorSampleMixin.handle_ancestor_input_samples_changed`, which auto-sets
the field when exactly one proband sample is upstream, and the proband patient when the samples make that ambiguous
(#1855) - this is why `_save` cascades `ancestor_input_samples_changed`.

`analysis/models/nodes/sources/cohort_node.py:CohortNode` filters by het/hom counts, a simple zygosity, or per-sample
zygosities (`CohortNodeZygosityFiltersCollection`), and adds the ref/het/hom count columns. Trio, duo and quad share
`analysis/models/nodes/family_inheritance.py:FamilyInheritanceNodeMixin` (inheritance-versus-family checks, waivable as
warnings) and one `AbstractFamilyInheritance` strategy object per mode, built by the node's `_inheritance_factory`;
`analysis/models/nodes/sources/trio_node.py:TrioNode`, `analysis/models/nodes/sources/duo_node.py:DuoNode` (#1829, a
proband and one parent, with Denovo becoming "absent in parent") and `analysis/models/nodes/sources/quad_node.py:QuadNode`
(a sibling too) differ only in the zygosity tuples. Compound het is
`analysis/models/nodes/family_inheritance.py:AbstractCompHetInheritance`: three queries find genes with a maternal-only
and a paternal-only hit (over `VariantGeneOverlap`, not transcript annotation, so a long SV VEP skipped still counts, #940),
memoised per node version, then the Q is "comp-het zygosity AND overlaps one of those genes". The mosaic-parent modes
(#1830) look for low-VAF alt reads in one parent through `analysis/models/nodes/family_inheritance.py:mosaic_evidence_q`.
`analysis/models/nodes/sources/pedigree_node.py:PedigreeNode` is the general case over a PED file: recessive and
dominant Qs from affected/unaffected zygosity sets via `CohortGenotypeCollection.get_zygosity_q`.
`analysis/models/nodes/sources/all_variants_node.py:AllVariantsNode` starts from every variant of the build, optionally
by contig and zygosity-count bounds from `VariantZygosityCountCollection`.

### Karyomapping

`analysis/models/models_karyomapping.py:KaryotypeBins` classifies each trio variant by the (proband, father, mother)
genotype tuple into in-phase/out-of-phase paternal and maternal bins. Creating a Trio queues
`analysis/tasks/karyomapping_tasks.py:create_genome_karyomapping_for_trio` (`snpdb/signals/signal_handlers.py:trio_post_save_handler`),
which streams the trio's genotypes once into `GenomeKaryomappingCounts` plus one `ContigKaryomappingCounts` per contig -
the relatedness summary shown for the trio. A `analysis/models/models_karyomapping.py:KaryomappingAnalysis`
(`analysis/views/views_karyomapping.py:create_and_view_karyomapping_analysis_for_trio`) adds `KaryomappingGene` rows
whose bins over the gene's flanking region drive the Plotly scatter and CSV download.

## Why it is shaped this way

A node is a filter, not a result set. `analysis/models/nodes/analysis_node.py:AnalysisNode.get_arg_q_dict` composes the
parent's `{alias: {hash: Q}}` dict with the node's own (`merge_arg_q_dicts`; the hash lets
`analysis/models/nodes/filters/merge_node.py:MergeNode` factor common filters out of an OR), and
`analysis/models/nodes/analysis_node.py:annotate_and_filter_queryset` applies each annotation alias then the Qs that
depend on it, then the alias-less Qs, so a chain of ten nodes is one SQL statement with each genotype partition joined
once. Materialising every node would mean writing millions of pks per edit on a 7.4M-variant cohort; composing means an
edit costs one count query per dirty node and the grid's page query is the same plan the count ran. The exceptions
(Venn, Intersection) are exactly the operations Postgres could not plan well as subqueries.

Caching is in tiers because each tier fails differently. The Redis Q cache (`ANALYSIS_NODE_CACHE_Q`, keyed on the
NodeVersion pk plus whatever `_get_cache_key` overrides fold in) exists because building the dict for a compound-het or
phenotype node costs seconds and the grid asks for it on every page; it is never invalidated, only outlived. NodeVersion
`load_data` is the count cache the badges read. NodeCache/VennNodeCache are on-disk sets for the two materialising
nodes. `cloned_from` is the fourth tier: a template run's nodes point at the snapshot's NodeVersions and reuse their
counts and grid cache until either side is edited, which is why `bump_version` clears it.

Explicit-pk substitution (#546): `analysis/models/nodes/analysis_node.py:AnalysisNode.get_arg_q_dict` answers with
`Q(pk__in=[...])` from `NodeVersion.variant_ids` when the node stored its pks at load, so a MergeNode or a Venn over
three 200-variant parents becomes a bitmap-or over the pk index instead of three nested subqueries, and the node's own
grid page is a pk lookup rather than its filter chain under 50 joins. The check sits ahead of the Redis Q cache, which
holds the real filter the load ran. Large nodes stay as `analysis/models/nodes/analysis_node.py:queryset_to_pk_in_q`, a RawSQL semi-join that also keeps the
dict picklable for Redis (a live `TransformerQuerySet` inside a Q is not, #240). The list is stored at load, after the
count, so "small" and "these pks" agree by construction; `analysis/tests/test_explicit_pk_substitution.py` pins the
two paths to the same answer.

Scheduling is single-worker and state-driven (#346) because the previous design - a celery group/chord built from a
toposort at edit time - deadlocked: `wait_for_node` tasks occupied every worker waiting for parents that could not get a
worker. Now nothing waits: the only cross-node dependency is "parents settled", read from the database by the one
dispatcher on `scheduling_single_worker`, so leases are written by one process and need no row-lock gymnastics, and a
finished node kicks the dispatcher rather than its children. The lease (`NodeTask.lease_expires`, restarted at claim)
is the only way a lost worker is detected; `MAX_NODE_ATTEMPTS` bounds retries and reclaims together. Workers never
`save()` a node (#431): `AnalysisNode.update` is a conditional UPDATE on (pk, version), so a user edit racing a load
wins by construction and the load discovers it through `NodeOutOfDateException`.

## History

The app began with a jqGrid front end and a toposorted chord of celery tasks per edit; `wait_for_node` and
`wait_for_cache_task` in `analysis/tasks/node_update_tasks.py` survive only for in-flight messages at deploy time. The
MergeNode rewrite (#240) stopped writing parent caches and introduced the hash-keyed Q dicts. Explicit pks (#546) and
the NodeVersion lease record (`analysis/migrations/0105_nodetask_lease_on_node_version.py`, #346) landed in mid 2026;
`analysis/migrations/0115_nodecount_variant_ids_nodeversion_live_data_sources_and_more.py` added the stored pk list and
count provenance (#235), and `analysis/migrations/0132_nodeversion_load_data.py` folded the per-label `NodeCount` rows
into `NodeVersion.load_data` (#1825, #1820). Node exports moved to celery-backed `CachedGeneratedFile`s (#1257). `VariantGrid.known_count` arrived with the
KnownCountPaginator (#1700); grids then moved from jqGrid to native DataTables (#1785, #1815).
The ClassificationsNode was split into Classifications and ClinVar nodes (#1789). Templates gained draft versions
(#1496). Per-tag node counts (#21), the TagNode editor's pill picker (#1820), sample-node grouping levels and per-sample
overrides (`analysis/migrations/0130_samplenodesamplefilter.py`), waivable field errors
(`analysis/migrations/0131_analysisnode_ignore_field_errors.py`), the DuoNode (#1829,
`analysis/migrations/0133_add_duo_node.py`) and the mosaic-parent modes (#1830) are the most recent additions.

## Traps

The Q cache is keyed only on NodeVersion pk (plus what `_get_cache_key` adds), so a test that changes a node and expects
a different queryset without saving needs `ANALYSIS_NODE_CACHE_Q=False`; a node whose query depends on something outside
its fields and its parents must fold that into `_get_cache_key` the way `CohortMixin` folds in the genotype collection,
or deleting and re-importing a VCF leaves a `FieldError` for an alias that no longer exists.
`analysis/signals/source_data_invalidation.py:_bump_nodes` is the pattern for bumping nodes when their source rows are
deleted: bump, cascade, then dispatch on commit so the node settles in ERROR_CONFIGURATION instead of spinning as DIRTY.

`_get_cached_label_count` reuses a parent's label count only when the parent was loaded with that label configured;
adding a node-count type after the fact makes every node run the SQL once. A count bigger than the parent's is a real
bug for a deterministic node: a filter that fans out over a multi-valued join (transcript annotation, gene lists) must
set `queryset_requires_distinct` or use a subquery, and `_get_variant_ids_to_store` will catch the case where the pk
list overruns the count. `AnalysisNode.load` persists only what `_load` returns and what `update()` is given - setting
attributes on `self` inside a task is lost.

`lease_ready_nodes` perma-fails a node found LOADING with a lapsed lease on the theory that it killed its worker
(OOM, SIGKILL); only QUEUED nodes are re-leased. A node whose load legitimately exceeds `LEASE_SECONDS` (10 min) will be
failed by the next backlog sweep even though the worker is still running it - the lease is not heart-beaten. `analysis/views/views_node.py:node_cancel_load` revokes the celery task, `pg_cancel_backend`s the recorded `db_pid`
and saves CANCELLED without bumping the version - a load that survives the revoke can still `update()` READY over it.

`Analysis.VERSION_BUMP_FIELDS` (custom columns, default sort) must bump `Analysis.version`
(`analysis/forms/forms.py:AnalysisForm`): the node editor, grid config and grid data are `cache_page`d under the analysis
version and node version in the URL, so a settings change that forgets the bump serves last week's columns. The grid
handler's cache lock is per user and URL: a second user hitting the same slow node runs the query again, and a gunicorn
worker killed inside the lock leaves it held for up to 10 minutes.

`VennNode.get_cache_task_args_set` sets `errors` and ERROR directly on the node with a `save()` if creating the cache rows
fails - the one place a node is saved from dispatch - and that status is never cleared automatically; the user must
re-save the node. `venn_cache_count` holds both sides' pk sets in worker memory, so a Venn over two whole-cohort parents
is the analysis app's most likely `MemoryError`. Cloning is `save_clone` per node in toposort order and each subclass
that owns satellite rows (`FilterNodeItem`, `TagNodeTag`, contigs, sample filters) must copy them itself; a new node
model with a related table that forgets `save_clone` silently loses its configuration in every template run.
