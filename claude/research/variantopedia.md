# variantopedia — research notes

Verified against a96540a68 on 2026-09-25

variantopedia owns no models (`claude/maps/models.md#variantopedia` is empty). It is the read side of the variant
database: the variant and allele pages, search, the All Variants and tag work-list pages, the tag stats page, the
nearby-variants summaries, and the superuser server status pages. Everything it shows lives in snpdb, annotation,
classification, analysis (VariantTag) and seqauto; this app composes it. URLs, the one task, the one signal receiver and
the two commands are in the generated maps ([urls](../maps/urls.md#variantopedia), [tasks](../maps/tasks.md#variantopedia),
[commands](../maps/commands.md), [signals](../maps/signals.md)); `claude/domain.md` has the vocabulary.

## Flows

### Variant page

`variantopedia/views.py:view_variant` picks a build (URL, else the user's default build if the variant is in it, else
`Variant.any_genome_build` - contigs such as MT sit in several builds), sets it on
`snpdb/genome_build_manager.py:GenomeBuildManager`, and hands off to
`variantopedia/views.py:variant_details_annotation_version` with the build's latest AnnotationVersion. That view is also
reachable directly with an explicit annotation version id, with the bare `variant_details.html` template - analysis
grids load it that way from `grid.js`. It builds
`annotation/transcripts_annotation_selections.py:VariantTranscriptSelections` for the transcript picker, the canonical
transcripts per enrichment kit (`variantopedia/views.py:get_genes_canonical_transcripts`), the VariantAllele/ClinGen data
(left out when `needs_clingen_call()` so the page fetches it asynchronously), tag presence across the allele's builds,
and the short label the page title uses. Gene-level variants (fusion, copy number, splice event) are resolved by
`_get_gene_fusion` / `_get_gene_copy_number_event` / `_get_splice_event_variant` in `variantopedia/views.py` and shown
in place of a coordinate. Tabs load separately: nearby variants, sample information
(`variantopedia/views.py:variant_sample_information` is a shell; the grid and counts are drawn client side from the
`variant_sample_genotypes` API), and the tag counts grid (`variantopedia/grids.py:VariantTagCountsColumns`) whose rows
expand to `variantopedia/grids.py:VariantTagDetailColumns`. `variantopedia/views.py:variant_grid_row_detail` is the
expanded row under any variant grid row (transcripts, classifications, ClinVar, build variants) - `grid.js` calls it.

### Allele page

`variantopedia/views_allele.py:view_allele` (also at `/a<id>`) shows the allele across builds: the
`classification/variant_card.py:AlleleCard`, and the valid classification `Overlap` rows for the allele split into
single-context and cross-context (`classification/models/overlaps_model.py:Overlap`, which replaced the discordance
reports / allele origin grouping descriptions). `view_allele_from_variant` redirects variant to allele only when
`PREFER_ALLELE_LINKS` is on. `create_variant_for_allele` (POST) re-runs liftover from the allele's imported (non-liftover)
build into the requested build, retrying every conversion tool; with no imported VariantAllele it silently just
redirects.

### Search

`variantopedia/views.py:search` is a thin shell over `snpdb/search.py:search_data` - the search machinery lives in snpdb.
Every call is logged with `create_event` (even an empty search, which is made so the page can list the enabled search
types), and a single preferred result redirects straight to it unless `mode=preview`. `?classify` swaps in
`variantopedia/forms.py:SearchAndClassifyForm`, which search passes on so results offer "create classification".

### All Variants page

`variantopedia/grids.py:AllVariantsGrid` over the build's latest annotation version. The chromosome / gene / variant type
/ min sample count controls are applied on a Search click (`variantopedia/templates/variantopedia/variants.html`); a
direct grid hit (CSV export, bookmark) falls back to the user's last saved filters
(`snpdb/variant_filters.py:get_all_variants_filters`). With no selective filter
(`snpdb/variant_filters.py:is_selective`) the grid matches nothing rather than scanning the table. Nothing is sortable:
every page is in `(contig, position, pk)` order, which streams off the `snpdb_locus` unique index (#1663).
`AllVariantsGrid.paging` takes the page's pks first and only then joins the annotation columns for those rows (#1887);
a gene filter is bounded by the gene's locus span (`snpdb/variant_filters.py:get_gene_bounds_q`) so it seeks.
`AllVariantsGrid.known_count` reports the planner's row estimate (`~1.2M`) instead of a COUNT(*) once it passes
`APPROXIMATE_COUNT_MIN`.

### Tags: work list, stats

`variantopedia/views.py:variant_tags` pairs two grids that share one filter state: `VariantTagsColumns` (one row per
tagging) and `TaggedVariantGrid` (one row per variant, with a tag-events count), both in `variantopedia/grids.py`. Both
start from `analysis/models/models_variant_tag.py:VariantTag.get_for_build`, so a tagging shows in a build once its allele
has a variant there, or immediately in the build it was made in (before liftover assigns the allele). The build's
coordinate is picked with a `FilteredRelation` on the allele's VariantAllele, coalesced to the tag's own variant
(`VariantTagsColumns._build_variant_field`), and annotation is joined through the build's partition
(`_annotation_version_queryset`). Resolved taggings (a classify-queue tag whose classification exists) are hidden unless
the tag config says show them or the page flips `show_resolved` for one view
(`variantopedia/grids.py:show_resolved_variant_tags`, #1901). Without `show_group_data` a user sees only their own
taggings; an explicit user filter overrides that but stays permission checked. The CSV export
(`variantopedia/views.py:variant_tags_export`) streams straight off the grid's queryset rather than paging it.

`variantopedia/views_tag_stats.py` (#1751) is a skeleton page whose cards each fetch a JSON endpoint cached in Redis per
user and parameters for a day (`_cached_json`). It counts on allele, so there is no build switch; tags with no allele yet
are a separate count. Only tags the user can see are counted, and an allele origin radio narrows every card.

### Nearby variants

`variantopedia/interesting_nearby.py:get_nearby_qs` gives querysets for the same codon, exon, InterPro domain, a
`VARIANT_DETAILS_NEARBY_RANGE` window and (if `VARIANT_DETAILS_NEARBY_SHOW_GENE`) each overlapping gene;
`get_nearby_summaries` / `interesting_counts` aggregate each into a text summary (variant count, classification and
ClinVar counts by clinical significance, DB zygosity sums) plus tag counts. `variantopedia/grids.py:NearbyVariantsGrid`
shows one region's queryset.

### Server status (superuser)

`variantopedia/views_server_status.py:server_status` pings celery workers by the names in `CELERY_WORKER_NAMES` (plus
analysis / heavy workers when enabled), lists long-running SQL (`library/django_utils/database_utils.py:long_running_sql`),
checks reference FASTA access, whether the highest variant is annotated, disk and integration status, and handles POST
actions (Test Slack, Health Check, Test Rollbar, run-integration, kill-pid via `signal_backends`). The tabs
(`server_status_activity`, `server_status_settings`, `health_check_details`, `database_statistics`) are loaded
separately. Health data comes from `library/health_check.py`: apps answer `health_check_signal` (recent activity) and
`health_check_overall_stats_signal` (totals); this app contributes
`variantopedia/signals/integration_health_check.py:integration_health_check`. The nightly Slack digest is
`variantopedia/tasks/server_status_tasks.py:notify_server_status` (beat at 19:00, `variantgrid/celery.py`), gated on
`HEALTH_CHECK_ENABLED`.

The two commands (`deployment_check`, `download_annotation_data`) are deploy tooling that live here by history.

## Why it is shaped this way

The app is a presentation layer so that the heavy logic stays next to its data: search in snpdb, annotation selection in
annotation, overlaps in classification. What is left here is page composition and grid query shaping, and the grid
decisions are all about the size of `snpdb_variant`: no user sorting, pk-then-columns paging, estimated counts, and a
match-nothing default on All Variants, each introduced after a page blew the statement timeout (#1279, #1651, #1663,
#1887). The tag grids key on allele rather than variant so a tagging follows the variant across builds, and the tag stats
page goes further and drops builds entirely, because counting per build double-counts lifted-over tags.

## History

The app dates from the 2020 repo "blank slate"; its last model-ish code, a `SearchTypes` constants class, went when
search types moved to the preview coordinator (search itself is `snpdb/search.py`). The
All Variants page gained its at-least-one-filter rule and genomic ordering in #1663 (2026-07) and two-phase paging in
#1887. Tag stats and tag merging arrived in #1751 (2026-08), then allele-origin-aware tags. Resolved-tag visibility
moved from a per-user grid setting to the tag config in #1901, and tag colours became the tag config (#1892). The allele
page switched from discordance reports to `Overlap` in the somatic squash (2026-09). Long-running query and backend
signalling helpers moved out of the server status view into `library/django_utils/database_utils.py` (#404).

## Traps

Viewing a page writes: `view_allele` calls `snpdb/clingen_allele.py:link_allele_to_existing_variants` (may create
VariantAlleles) and both pages call `annotation/models/models.py:ClinVarRecordCollection.set_allele_for_variants`, which
backfills `allele` on ClinVar record collections so the ClinVar AJAX finds them. A read-only replica would break these
pages.

`variant_details_annotation_version` swallows any exception from transcript selection with `log_traceback` (an
unannotated variant is expected to fail there), so a real bug in that path shows up as a page missing its transcript
panel, not a 500 - check the logs.

Confirmed bugs (not fixed):
- Nearby tag counts are inflated. `variantopedia/interesting_nearby.py:interesting_counts` puts the tag `StringAgg`
  in the same aggregate as the classification count join, so each tag is repeated once per classification on the allele
  (verified read-only: allele 54409, 2 tags, 2 classifications, counted 4). It also counts every user's tags, not only
  those visible to the viewer (`VariantTag.filter_for_user` is not applied), unlike the tag grids and tag stats.
- `variantopedia/interesting_nearby.py:get_transcripts_and_codons`, `get_transcript_and_exons` and
  `get_transcript_and_domains` read the source variant's `varianttranscriptannotation_set` with no version filter, so
  rows from historical VAVs feed the codon/exon/domain lookups and the method summaries.
- `variantopedia/views_server_status.py:server_status`: the "Not annotated, no AnnotationRun!" branch catches
  `AnnotationRun.DoesNotExist`, which a `filter()` never raises; an unannotated highest variant with no run shows a
  warning reading `AnnotationRuns: ` instead. The POST actions also render rather than redirect (a TODO in the code), so
  a browser refresh repeats Test Slack / kill-pid.
- `variantopedia/tasks/server_status_tasks.py:notify_server_status` is a bare `celery.shared_task` with no queue; it
  lands on the default `db_workers`, against the rule that every task names its queue.
- The `gene` parameter of `variantopedia/grids.py:VariantTagsColumns.filter_queryset` is sent by no page, and it filters
  on the tag's own variant rather than the build variant, so it would drop taggings made in another build.
