# #1892 - TagColor becomes TagConfig, owned by the analysis

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-22
Status: landed 2026-09-22 (master)

## Why

A `TagColorsCollection` row per tag now holds colour, sort order and the 1-click flag (#343, #1888), so the
collection is a lab's "how we use tags" configuration, and the name misleads: it made the stale-days window
(#1433) look like it belonged on user settings and the analysis rather than here. Renaming to **TagConfig**
matches `CustomColumnsCollection` / `CustomColumn`, moves stale days onto the collection, and lets an analysis
point at one collection so everyone opening it sees the same colours, order, quick tags and staleness. Custom
columns already work this way (`analysis/models/models_analysis.py:Analysis` has `custom_columns_collection`,
copied from user settings at creation).

## Models

`snpdb/models/models_user_settings.py`

```python
class TagConfigCollection(GuardianPermissionsAutoInitialSaveMixin, TimeStampedModel):   # was TagColorsCollection
    user = models.ForeignKey(User, null=True, blank=True, on_delete=CASCADE)
    name = models.TextField()
    version_id = models.IntegerField(null=False, default=0)
    # Moved from SettingsOverride / Analysis - one window per collection, so it is set in the same place
    # as the colours and quick tags and shared by every user of an analysis that uses the collection
    variant_tag_stale_days = models.IntegerField(null=True, blank=True, choices=VARIANT_TAG_STALE_DAYS_CHOICES,
                                                 help_text="Tag events older than this are considered stale: grids show "
                                                           "fresh vs total counts and mark tags whose most recent event is "
                                                           "older. Blank disables staleness.")


class TagConfig(TimeStampedModel):   # was TagColor
    collection = models.ForeignKey(TagConfigCollection, null=True, on_delete=CASCADE)
    tag = models.ForeignKey(Tag, on_delete=CASCADE)
    rgb = models.CharField(max_length=7)
    sort_order = models.IntegerField(null=True, blank=True)
    quick_tag = models.BooleanField(default=False)

    class Meta:
        unique_together = ('collection', 'tag')


class SettingsOverride(models.Model):
    ...
    tag_config = models.ForeignKey(TagConfigCollection, on_delete=SET_NULL, null=True, blank=True,   # was tag_colors
                                   help_text="Tag colours, sort order, 1-click tags and staleness (modify/create these in 'Tag settings'). "
                                             "Initial tag config when creating an analysis")
    # variant_tag_stale_days removed - lives on TagConfigCollection


@dataclass
class UserSettings:
    ...
    tag_config: TagConfigCollection   # was tag_colors
    # variant_tag_stale_days removed; variant_tag_stale_date property reads tag_config
```

`VARIANT_TAG_STALE_DAYS_CHOICES` moves from `SettingsOverride` to module level next to `TagConfigCollection`
(the analysis form and the collection page both need it).

`analysis/models/models_analysis.py`

```python
class Analysis(...):
    VERSION_BUMP_FIELDS = ["custom_columns_collection", "default_sort_by_column", "tag_config_collection"]
    ...
    tag_config_collection = models.ForeignKey(TagConfigCollection, null=True, blank=True, on_delete=SET_NULL)
    # variant_tag_stale_days removed - comes from tag_config_collection
```

Resolution order everywhere a page needs tag config: the analysis' `tag_config_collection` when the page is
inside an analysis and it is set, otherwise the viewer's `UserSettings.tag_config`. `variant_tag_stale_date`
on both `Analysis` and `UserSettings` derive from the resolved collection's `variant_tag_stale_days`
(None when there is no collection or the field is blank).

## Migrations

All existing migrations stay as they are (`0211`, `0212`, `0272` and `analysis` `0147` are pushed).

1. `snpdb/migrations/0273_rename_tagcolor_tagconfig.py` - `RenameModel` TagColorsCollection → TagConfigCollection,
   `RenameModel` TagColor → TagConfig, `RenameField` SettingsOverride.tag_colors → tag_config,
   `AddField` TagConfigCollection.variant_tag_stale_days, then `RunPython` that copies stale days onto collections:
   for every `SettingsOverride` with both `tag_config` and `variant_tag_stale_days` set, in override order
   (GlobalSettings, then Organization, Lab, User overrides), set the collection's value so the most specific
   override wins, and leave a collection alone once a user-level override has written it. Then
   `RemoveField` SettingsOverride.variant_tag_stale_days. Precedent for the rename:
   `snpdb/migrations/0091_rename_usertagcolors_tagcolor.py`.
2. `analysis/migrations/0148_analysis_tag_config_collection.py` - `AddField` `tag_config_collection`, then
   `RunPython` that sets it for every existing analysis from its owner's resolved settings (walk
   GlobalSettings → the owner's default lab's organization override → lab override → user override, last
   non-null `tag_config` wins - the same order as `UserSettings.get_settings_overrides`, written against the
   historical models). Then `RemoveField` `variant_tag_stale_days`.

Bump `CACHE_VERSION` in `variantgrid/settings/components/default_settings.py`: `UserSettings` is cached in
Redis and carries the collection.

## Rename surface

Everything named for colours becomes tag config. Names to change, with the new spelling:

- URLs in `snpdb/urls.py`: `tag_color_collections_datatable` → `tag_config_collections_datatable`,
  `view_tag_colors_collection` → `view_tag_config_collection`, `clone_tag_colors_collection` →
  `clone_tag_config_collection`, `set_tag_color` stays (it sets a colour). Path segment
  `settings/tags/collection/...` stays. Kwarg `tag_colors_collection_id` → `tag_config_collection_id`.
- `snpdb/views/views_tag.py` (`view_tag_colors_collection`, `set_tag_color`, `tag_settings`),
  `snpdb/views/views_json.py:clone_tag_colors_collection`, `snpdb/grids.py:TagColorsCollectionColumns`.
- `snpdb/templates/snpdb/settings/view_tag_colors_collection.html` → `view_tag_config_collection.html`; its
  copy says "columns" in two places (the clone button and the permission warning in the view) - say tag config.
  Add a stale-days select (the `VARIANT_TAG_STALE_DAYS_CHOICES` plus blank) saved through the same POST handler
  as the name, and `increment_version` on change.
- `snpdb/templates/snpdb/settings/tag_settings.html`: heading "Tag colours & sorting" → "Tag config", copy
  mentions colours, sort order, 1-click tags and staleness.
- `snpdb/utils.py`: `_resolve_tag_colors_collection` → `_resolve_tag_config_collection`, the
  `tag_colors_collection=` kwarg → `tag_config_collection=` on `get_tag_sort_order_by_tag`,
  `get_tag_quick_tags`, `get_all_tags_and_user_colors`, `get_tag_styles_and_colors`.
- `analysis/templatetags/user_tag_color_tags.py` → `tag_config_tags.py` (every `{% load user_tag_color_tags %}`
  follows - 17 templates, listed by `grep -rl user_tag_color_tags`). `tag_colors_collection_link` →
  `tag_config_collection_link` and its template `analysis/templates/analysis/tags/tag_colors_collection_link.html`
  is renamed to `tag_config_collection_link.html`.
- `snpdb/forms.py` label `"tag_colors": "Tag Colours"` → `"tag_config": "Tag Config"`; the
  `variant_tag_stale_days` label and visibility entries go.
- `analysis/fake_variant_tags.py`, `snpdb/tag_operations.py` (`TagForeignKey("tag config", TagConfig, ...)`),
  `snpdb/models/models_user_settings.py` module docstring.
- `snpdb/management/commands/fix_tag_colors_collection_permissions.py` is a 2022 one-off that says it can be
  removed once past migration 0089: delete it.
- Changelog `variantgrid/templates/default_templates/changelog.html`: add a line for #1892 (rename, stale days on
  the collection, analysis-level tag config). Existing changelog lines keep their historical wording.
- `claude/domain.md`: add TagConfigCollection / TagConfig entries if tags are listed there.
- Tests: `snpdb/tests/test_tag_operations.py`, `snpdb/tests/test_tag_sort_order.py`, `snpdb/tests/test_urls.py`,
  `snpdb/tests/test_variant_tag_stale_setting.py`, `variantopedia/tests/test_tagged_variant_grid.py`,
  `analysis/tests/test_models.py` (the "copied from user settings" test now asserts `tag_config_collection`).

## Analysis uses its own collection

- `Analysis.set_defaults_and_save` copies `user_settings.tag_config` into `tag_config_collection` (next to
  `custom_columns_collection`).
- `analysis/forms/forms.py:AnalysisForm` lists `tag_config_collection` in place of `variant_tag_stale_days`;
  the analysis settings tab (`analysis/templates/analysis/analysis_settings_details_tab.html`) appends a
  "Manage Tag Config" crosslink after the select, the same way it does "Manage Custom Columns".
- `analysis/views/views.py:get_analysis_settings` keeps the `variant_tag_stale_days` key (the JS reads it from
  `ANALYSIS_SETTINGS`) and sources it from the analysis' resolved collection.
- The template tags that render tag CSS, order and quick tags take the collection from the analysis when the
  page has one: `render_tag_styles_and_formatter`, `render_variant_tag_order`, `render_variant_quick_tags` and
  `render_node_count_styles` accept an optional `analysis` argument and resolve
  `analysis.tag_config_collection` first. Analysis pages pass it: `analysis/templates/analysis/analysis.html`,
  `analysis_editor_and_grid.html`, `analysis_settings_node_counts_tab.html`, `classify_report_tab.html`.
  `analysis/models/nodes/node_counts.py:get_tag_node_count_colors` takes the same optional collection.
- The Classify & Report tab (`analysis/views/views_classify_report.py`) has no analysis of its own: when every
  tagging in the queue came from one analysis (`analysis/classify_report.py:tag_config_analysis`) it draws with
  that analysis' collection, otherwise the viewer's, and `tag_config_collection_link` says which under the
  "Tags to action" card.
- `analysis/grid_export.py` already reads `node.analysis.variant_tag_stale_date`; that property now resolves via
  the collection. Outside an analysis (`variantopedia/grids.py`, `variantopedia/views.py`) the viewer's
  `UserSettings.variant_tag_stale_date` is used as now.

## Crosslinks on user settings

`snpdb/templates/snpdb/settings/view_user_settings.html` appends, after the `columns` and `tag_config` selects,
a link to the selected collection's page (`get_absolute_url`), updated on change, and a link to the
listing page (`custom_columns`, `tag_settings`) when nothing is selected. Reuse the `createCrossLink` /
`setCrossLink` pattern from `variantgrid/static_files/default_static/js/global.js`. Lab and organization
settings pages use the same `settings_override` template tag, so put the JS in
`snpdb/templates/snpdb/tags/settings_override.html` so all three levels get it.

## Tests that earn their keep

- Analysis created via `set_defaults_and_save` gets the user's `tag_config_collection` (replaces the stale-days
  copy test in `analysis/tests/test_models.py`).
- `Analysis.variant_tag_stale_date` is None without a collection, and derives from the collection's days when set.
- Template tag resolution prefers the analysis' collection over the viewer's user settings.
- Migration data copy: a `SettingsOverride` chain where the lab and user both set stale days ends with the user's
  value on the collection (use `TestMigrations`-style or a plain unit test of the helper the migration calls).

## Definition of done

`scripts/vg tests --explain --run` passes, `scripts/vg docs check` passes, `scripts/vg map` refreshed,
`vg page /analysis/<id>` and `/settings/tags` render, `Status:` above updated.
