# #1433 — Analysis tags: out of date tags ("tag freshness")

[#1433](https://github.com/SACGF/variantgrid/issues/1433): medical scientists periodically review
tags (e.g. a rare MPL variant tagged `NOTSomaticReportable` years ago) and want to see at a glance
whether a tag still reflects current group opinion. apgmp's comment frames the winning design:
show **total plus fresh** counts, rendered as `SomaticReportable x200 (42 fresh)` — 200 tag events
all-time, 42 within the freshness window — plus a visual cue when a tag's most recent event is
older than a user-configurable cutoff. Total leads so the familiar number keeps its position; the
bracketed part is the freshness annotation.

Related: `claude/plans/1751_tag_stats_plan.md` — its SA Path export analysis showed 45% of tag
events restate an existing (variant, tag) pair, and that different-user re-confirmation is signal
(group consensus), not noise. That directly supports the core decision below.

## Decisions

- **No new curation/review model.** Re-tagging *is* the curation act: tagging the variant again in
  a current analysis bumps freshness naturally and preserves history. `VariantTag` is a
  `TimeStampedModel` (`analysis/models/models_variant_tag.py:39`) so every tag event already
  carries `created` — no backfill needed. The issue's "bump date / add review audit" options are
  rejected in favour of this.
- **The staleness cutoff is a `SettingsOverride` field**, so it cascades
  Global → Organization → Lab → User with later levels overriding when non-null
  (`snpdb/models/models_user_settings.py:186`). A lab sets its policy (e.g. 730 days); individual
  users can override on their settings page.
- **Null cutoff = feature off** — displays stay exactly as today, so deployments that don't tag
  (e.g. Shariant) see no change.
- **Freshness is computed from `VariantTag.created`** at display/query time. No denormalisation:
  tag volumes per variant are small and the grid already aggregates per row via subquery.

## Deliverables (independently shippable, in order)

### A. `variant_tag_stale_days` setting

1. Add to `SettingsOverride` (`snpdb/models/models_user_settings.py:186`):
   `variant_tag_stale_days = models.IntegerField(null=True, blank=True, help_text=...)` —
   "Tag events older than this many days are considered stale: grids show fresh vs total counts
   and mark tags whose most recent event is older. Blank inherits the next level up / disables."
   Follow the `node_grid_auto_load_max_variants` pattern (models_user_settings.py:240) — it was the
   most recent field added the same way.
2. Add the field to the `UserSettings` dataclass (`snpdb/models/models_user_settings.py:403`) —
   `variant_tag_stale_days: Optional[int]`. The cascade in `UserSettings.get_for()` picks it up by
   field name automatically.
3. `SettingsOverrideForm` (`snpdb/forms.py:518`): add a label ("Variant Tag Stale Days") and gate
   visibility on `settings_config.analysis_enabled` in `_hide_unused_fields`, same as `tag_colors`.
   `fields = "__all__"` means the org/lab/user forms pick it up automatically.
4. snpdb migration.
5. Convenience accessor: `UserSettings.variant_tag_stale_date` property returning
   `timezone.now() - timedelta(days=...)` or `None`, so callers share one definition of "stale".

### B. Fresh/total counts + stale styling in the analysis grid

1. `get_variantgrid_extra_annotate` (`snpdb/grid_columns/custom_columns.py:78`): change the
   `tags_global` `StringAgg` payload from `tag_id` to `tag_id` + `:` + ISO date, i.e.
   `StringAgg(Concat("tag_id", Value(":"), Cast(TruncDate("created"), TextField())), delimiter=Value("|"))`
   producing `Artefact:2024-03-01|Artefact:2019-07-12|SomaticReportable:2023-05-06`.
2. `tagsGlobalFormatter` (`variantgrid/static_files/default_static/js/grid.js:509`): split each
   entry on `:`, count fresh (date ≥ cutoff) vs total per tag. Render total-first:
   - cutoff set: `Tag x200 (42 fresh)` with a tooltip like
     "42 of 200 tag events within the last 2 years, most recent 2023-05-06"
   - no cutoff configured: current rendering (`Tag x 200`)
   - tag entirely stale (newest event < cutoff): reduced opacity on the tag span plus a clock
     glyph, tooltip showing the most recent date.
   The cutoff reaches JS the same way `variantTagsReadOnly` does (read off the analysis window at
   grid.js:502) — set a `variantTagStaleDays` JS variable in the analysis template from the view
   context (`analysis/views/views.py:176` already builds `user_settings` for that context).
3. Export path `format_items_iterator` (`analysis/grid_export.py:224`): update the parser for the
   new `tag:date` payload — count on the tag part, and emit fresh counts when the analysis user's
   cutoff is set: `Artefact x 200 (42 fresh)`. This function feeds both CSV and VCF exports, and its
   test (`analysis/tests/test_grid_export.py:176`) pins the old format — update it.
4. `snpdb/grids.py:501` (`tags_global` column def) and `snpdb/vcf_export_columns.py:156` consume
   the same field — verify both render via the updated formatter/iterator rather than raw values.

### C. Stale display on the variant page

`VariantTagCountsColumns` (`variantopedia/grids.py:362`) already returns `count` and
`last_created` (via `get_variant_tag_counts_qs`, `analysis/models/models_variant_tag.py:115`).
1. Annotate a `fresh_count` (`Count("id", filter=Q(created__gte=stale_date))`) when the requesting
   user's cutoff is set, and add it as a column.
2. `tagRenderer` client renderer: same opacity + clock treatment as the grid when
   `last_created` is older than the cutoff (pass the cutoff into the page context).

### D. Date filter on the Tags node

The user-facing input is a rolling window in days (so analysis templates stay meaningful), but the
cutoff is anchored to the **node's save date, not query-execution time**: the node already records
when each version was saved (`NodeVersion.created` — `get_warnings` at
`analysis/models/nodes/filters/tag_node.py:28` looks it up for the snapshot warning). Anchoring to
save date makes the filter deterministic per node version — re-running an old analysis reproduces
what it showed when saved (audit of what the analysis was at that time), matching the existing
"global tags are a snapshot … press save to refresh" semantics. No extra snapshot field or
migration column is needed beyond the days value itself.

1. `TagNode` (`analysis/models/nodes/filters/tag_node.py:16`): add
   `tagged_within_days = models.IntegerField(null=True, blank=True)` + analysis migration.
2. `_get_node_q()` (tag_node.py:44): when set, look up this version's `NodeVersion.created`
   (fall back to `timezone.now()` if missing), compute
   `cutoff = node_version_created - timedelta(days=tagged_within_days)`, and filter both
   querysets — `variants_with_tags.filter(created__gte=cutoff)` for this-analysis tags, and
   `tags_qs.filter(created__gte=cutoff)` before `variants_for_build_q` in
   `TagNodeMode.ALL_TAGS`.
3. `TagNodeForm` (`analysis/forms/forms_nodes.py:873`): expose the field with label
   "Only tags added within (days)"; hide it for the invisible analysis-tags node alongside
   `mode`/`node_input`.
4. `get_node_name()` / `_get_method_summary()`: append e.g. `≤ 730d` when set; the method summary
   should also state the resolved cutoff date so grids/exports record the actual date filtered on.

### E. Staleness filter on the tag stats page (#1751)

Depends on the stats page from `claude/plans/1751_tag_stats_plan.md` landing (page of cards, each
lazy-loading a JSON endpoint, responses cached in Redis keyed by a hash of endpoint parameters +
genome build).

1. Page-level "freshness window (days)" control, defaulting to the requesting user's
   `variant_tag_stale_days` (blank when the user has none set = no filtering, today's behaviour).
2. The value is passed as a request parameter to each card's JSON endpoint and **included in the
   parameter hash that keys the Redis cache**, so different windows cache independently and a
   user's default window gets warm-cache behaviour like any other parameter. No cache-version bump
   needed — new parameter, new hash.
3. Endpoints apply it as `VariantTag.created >= now - window` on their base queryset. Cards where
   it earns its keep: headline strip (fresh totals alongside all-time), genes, mega-artefacts, and
   the tagged-variants grid link-outs. Time-series cards (tags over time, tag's genes over time)
   already display the date axis and ignore the parameter.
4. Changing the control re-fetches the cards (same lazy-load path with the new parameter), showing
   the per-card spinner/Recalculate pattern the stats plan already defines.

## Out of scope

- `VariantTag.get_summary_for_variant` (`analysis/models/models_variant_tag.py:122`) feeds
  classification evidence autopopulation (`classification/autopopulate_evidence_keys/evidence_from_variant.py:228`)
  — its text stays unchanged so classification evidence wording isn't perturbed.
- The variant preview tag counts (`analysis/templates/analysis/tags/tag_counts.html`, via
  `snpdb/signals/variant_zygosity_preview_extra.py`) — could get the same treatment later, but the
  hover-preview is transient and the main grids cover the review workflow.
- Any per-`Tag` staleness config (e.g. different cutoffs for `Artefact` vs `SomaticReportable`).
  Start with one cutoff per user/lab; revisit if requested.

## Tests

- Node: `TagNode.tagged_within_days` filters old tag events out (this-analysis and global modes),
  and the cutoff anchors to the `NodeVersion` save date, not query time — extend the existing tag
  node tests.
- Export: `test_grid_export.py` — new `tag:date` payload parses to `Tag x 5 (2 fresh)` with a
  cutoff and `Tag x 5` without.
- `UserSettings.variant_tag_stale_date` — None when unset, correct date when set (the cascade
  itself is existing generic machinery; don't re-test it per level).
