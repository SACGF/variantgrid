# One-click tagging for the tags a lab uses most

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-21
Status: in progress

[#1888](https://github.com/SACGF/variantgrid/issues/1888) (Tagging - very common tags). Today every tagging in the
analysis grid is two steps: click the (+) in the tags column, then pick the tag in the select2. Tag use is Pareto
shaped - a somatic analysis is nearly always SomaticReportable / NOTSomaticReportable / Artefact - so the issue asks
for a per-tag "1 click" option: keep the (+) for the full autocomplete, and beside it draw one (+) per quick tag, in
the tag's colour and in sort order, that tags on a single click.

---

## Data

Quick tags are a property of a tag colours collection, next to colour and sort order. That is where tags are already
customised per user, lab or deployment: `snpdb/models/models_user_settings.py:TagConfigCollection` resolves through
user settings (user → lab → organization → global), is shareable via guardian permissions, and is cloned to be
customised. A somatic lab's collection marks its three tags; a germline lab's marks its own; someone working both
sides picks the collection. A `TagConfig` row already exists purely to hold `sort_order` with no colour set, so a row
that exists purely to hold `quick_tag` follows the same rule.

```python
class TagConfig(TimeStampedModel):
    collection = models.ForeignKey(TagConfigCollection, null=True, on_delete=CASCADE)
    tag = models.ForeignKey(Tag, on_delete=CASCADE)
    rgb = models.CharField(max_length=7)  # '#rrggbb', '' when the row only holds sort_order / quick_tag
    sort_order = models.IntegerField(null=True, blank=True)
    # Drawn as its own (+) in the analysis grid's tags column, so tagging with it is one click (#1888)
    quick_tag = models.BooleanField(default=False)

    class Meta:
        unique_together = ('collection', 'tag')
```

Migration: snpdb/migrations/0272_tagcolor_quick_tag.py (new), an `AddField` with the default. No data migration - nothing
is quick until someone ticks it.

A merged-away tag's `TagConfig` rows are deleted by `snpdb/tag_operations.py` (`TagForeignKey("tag colour settings", ...)`
is unique with the collection so cannot be repointed), so a quick flag dies with the tag it was on; the survivor is
ticked again by hand if wanted. A retired tag is excluded when the list is built (below), so its button disappears
without touching the row - reinstating brings it back.

## What reaches the grid

`TagConfigCollection.get_quick_tags() -> list[str]`: the collection's tag ids with `quick_tag=True` whose tag is
live (`Tag.live_qs()`), ordered by `(sort_order or 0, tag id)` - the same rule `snpdb/utils.py:get_all_tags_and_user_colors`
and `sortVariantTags` in `variantgrid/static_files/default_static/js/grid.js` use, so the buttons sit in the same order
the pills do.

`snpdb/utils.py:get_tag_quick_tags(user, tag_config_collection=None) -> list[str]` beside
`snpdb/utils.py:get_tag_sort_order_by_tag`, resolving the collection the same way and returning `[]` with no collection.

A new template tag `render_variant_quick_tags` in `analysis/templatetags/tag_config_tags.py`, the twin of
`render_variant_tag_order`, emitting the JSON list. `analysis/templates/analysis/analysis.html` sets
`variantQuickTags = {% render_variant_quick_tags %};` next to `variantTagOrder`. The sample variants tab
(`snpdb/templates/snpdb/data/sample_variants_tab.html`) is read only, so it gets nothing - the formatter draws no
buttons when `readOnly` is set anyway.

## Grid

`VariantGridFormat.tags` in `variantgrid/static_files/default_static/js/variantgrid_formats.js`: when the grid is not
read only, after the existing (+) and before the pills, one button per entry of `aWin.variantQuickTags`:

```html
<a class='quick-tag' href='javascript:void(0)' variant_id='123' tag_id='Artefact' title='Tag as Artefact'>
    <span class='grid-tag tagged-Artefact'><span class='user-tag-colored quick-tag-button'></span></span>
</a>
```

The inner markup is the same as a pill's, so the existing per-tag CSS (`.tagged-<tag> > .user-tag-colored`, from
`render_tag_styles_and_formatter`) colours it with no new per-tag rules. `.quick-tag-button` in
`variantgrid/static_files/default_static/css/global.scss` (in the `table.variantgrid-datatable` block beside
`.add-variant-tag`) makes it a 16px circle with the Font Awesome `\f055` plus glyph, sized and aligned like the plain (+),
with a grey fallback background for a quick tag that has no colour set - the tag colour and its contrasting text colour
then override it. It is an `<a>`, so the row-click guard in `variantgrid/static_files/default_static/js/datatable_definition.js`
leaves it alone the same way it does the plain (+).

Tag ids are alphanumeric (the tag settings page enforces it on creation), but the button carries the tag as an
attribute rather than inlined into a `javascript:` href, and `escapeHtml` is used for the title as `getVariantTagHtml`
does.

Click handling in `variantgrid/static_files/default_static/js/grid.js`: `gridCompleteExtra` binds a click handler on
`.quick-tag` alongside `tagClickHandler`, reading `variant_id` / `tag_id` from the anchor and the node from
`closest("table.grid").attr("node_id")` exactly as `showTagAutocomplete` does. The success path - draw a new pill in
the cell only when `response.created`, via `variantTaggingPillOptions` and `getVariantTagHtml`, and bind
`tagClickHandler` on it - is lifted out of `showTagAutocomplete` into one function both callers use, so a quick tag
and an autocomplete tag produce the identical pill. `addVariantTag` / `setVariantTag` in
`variantgrid/static_files/default_static/js/analysis.js` are unchanged: the POST to `set_variant_tag` already
does `get_or_create` on the node's proband, so a second click on a tag the row already has is a no-op
(`created` false, no second pill) and the tag counts still refresh.

The buttons are always drawn, whether or not the row already has the tag: a tagging's identity includes whose it is
(sample / patient), so "already tagged" is a per-proband question the formatter would have to re-derive, and a cell
whose buttons come and go would shift the pills about. The cost is one more click that does nothing.

## Tag colours page

`snpdb/templates/snpdb/settings/view_tag_config_collection.html` gets a third column, "1-click", a checkbox per tag
row (disabled without write permission, like the colour picker). Ticking posts `quick_tag=<tag>&value=true|false` to
the page's own URL; `snpdb/views/views_tag.py:view_tag_config_collection` gains a branch beside `name` and `tag_order`
that `update_or_create`s the `TagConfig` row with `defaults={"quick_tag": ...}`, ignoring an unknown tag the same way
`tag_order` does. The page text below the heading says what the column does: "Tick 1-click to draw the tag as its own
(+) in the analysis grid, so tagging with it takes one click".

The tag settings page (`snpdb/templates/snpdb/settings/tag_settings.html`) already links each user to their collection;
nothing changes there.

## Tests

`snpdb/tests/test_tag_sort_order.py` already covers the collection's ordering and the view's POST branches; add to it
(or a sibling `TagQuickTagTest` in the same file):
- `get_quick_tags` returns only ticked tags, in sort order then name, and drops a retired tag.
- posting `quick_tag` via the view creates the row (with `rgb=''`), and posting `false` clears it; the row with no
  colour still emits no colour (the existing `test_sort_order_only_rows_emit_no_color` pattern).
- `get_tag_quick_tags` with no collection returns `[]`.

`scripts/vg tests --explain` after the change names the rest. `vg page /analysis/<id>` and the tag colours page before
and after, for query counts.

## Out of scope

- The variant details page's "Add Tag" (`variantopedia/templates/variantopedia/variant_details.html`) keeps its
  autocomplete - the issue is about the grid column, and that page tags once per visit rather than row after row.
- Filtering the quick tags by the analysis's allele origin: an analysis carries no allele origin bucket, and the
  collection is already the per-lab choice. If a mixed lab wants it, `Tag.allele_origin_bucket` is there to filter on
  later.

## Decisions

- Quick tags live on the tag colours collection, not on `Tag` or `UserSettings` - the collection is already how tags
  are customised per user / lab / deployment, and it is shared and cloned as one thing.
- Drawn as coloured (+) circles rather than mini named pills, as the issue asks; the tag name is the tooltip.
- Always drawn, never hidden for a row that already has the tag (see Grid).
