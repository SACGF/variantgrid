# #1751 — variant tag stats page, plus tag merging

[#1751](https://github.com/SACGF/variantgrid/issues/1751): tagging is one of the biggest workflows
in the system (~366 tags/day on VG3 prod) but has no aggregate view. Planning was informed by a full
export of SA Path's GRCh38 tags (`claude/variant_tags_export_2026-08-18.csv`, 387,810 rows,
May 2019 – May 2026, exported from the VG4 test clone). Related:
[variantgrid_private#1805](https://github.com/SACGF/variantgrid_private/issues/1805) (artefact
management / duplication).

Two independently shippable deliverables:

- **A. Tag stats** — a stats page linked from `variantopedia/variant_tags/`, plus small additions to
  the existing variant tags page, the tagged-variants grid, and the user page.
- **B. Tag merging** — admin merge UI on `snpdb/settings/tags`, case-insensitive creation
  validation, and an SA Path data migration merging `artefact` → `Artefact`.

## What the data showed (drives the design)

- 41 tag names, dominated by dismissal/triage tags: `PopFreqTooCommon` 68k, `NOTSomaticReportable`
  55k, `Artefact` 49k, `SomaticReportable` 42k, `SingleVariantARGene` 41k. Two clear workflows —
  germline (trio/singleton) and somatic haem — with somatic growing from ~10% (2021) to ~36% (2025).
- 79 users, heavily specialised by workflow; top 10 users ≈ 60% of tags.
- **Duplication is the headline.** Only 214,265 distinct (variant, tag) pairs in 387,810 tag events
  — 45% of tagging restates a (variant, tag) that already exists. The top variant
  (`20:32434638 AG>A`, the ASXL1 c.1934dup homopolymer) has been tagged `Artefact` **4,599 times**;
  JAK2 V617F tagged `SomaticReportable` 1,807 times. Every count on the stats page must therefore
  say which identity it counts: *tag events*, *distinct (variant, tag)*, or *distinct variants*.
- Within a single analysis, the same (variant, tag) recurs 47,994 times — but only 3,224 of those
  are the same user (true dupes); 44,770 are a *different* user independently re-confirming, which
  is signal (inter-user agreement), not noise.
- `TagNode` already supports excluding tagged variants (`TagNodeInput.PARENT_NOT_TAGGED`,
  `analysis/models/nodes/filters/tag_node.py:63` — renders as "Exclude Global/Analysis <tags>").
  Labs are either unaware of it or their workflow is to deliberately re-confirm each artefact per
  case. No node work is needed for #1805; the stats page surfaces the numbers that let labs decide
  (see the mega-artefacts card below).
- Exactly one case-collision pair exists at SA Path: `artefact` (404) vs `Artefact` (48,705). No
  Levenshtein-distance-1 pairs — but other deployments have them (the #1805 reporter's tag list
  contains a literal `Aretefact` typo), so merge suggestions earn their keep.
- The somatic team classifies somatic variants outside VariantGrid, so `SomaticReportable` tags
  having no `Classification` is expected — any tag→classification comparison must scope to
  germline tags only. (Deferred entirely; not in this plan's widgets.)

## A. Tag stats

### Architecture

Sizing reality: 400k `VariantTag` rows is small. Tag/user/month group-bys are sub-second;
gene-level and co-occurrence joins land in the seconds range. The design guards against stacking a
dozen queries into one synchronous view, not against any single slow query:

- **The page is a skeleton of cards; each card lazy-loads its own JSON endpoint** on page load
  (spinner per card). First paint is instant.
- **Every endpoint caches its JSON in Redis** (`django.core.cache`), keyed by a hash of its
  parameters + genome build, ~24h timeout. Responses include a `calculated` ISO timestamp and each
  card shows it subtly ("as of 3 hours ago"), so staleness is visible rather than mysterious.
- Charts render client-side with **plotly** — copy the pattern from
  `variantopedia/templates/variantopedia/database_statistics_detail.html:113`
  (`js/lib/plotly-latest.min.js` + `js/plotly_helper.js`).
- Interactive cards (gene picker, co-occurrence, tag-by-gene-over-time) recompute on a
  "Recalculate" button with a spinner; results cached per parameter hash.
- Genome build handling follows the existing page: `UserSettings.get_genome_build_or_default` +
  `VariantTag.get_for_build` (`variantopedia/views.py:479`).
- New views in `variantopedia/views_tag_stats.py` (`views.py` is already ~1,900 lines), URLs under
  `variantopedia/urls.py` next to the existing `variant_tags` entries: page at
  `variant_tags/stats/<genome_build_name>` plus one JSON endpoint per card.
- Gene attribution for a tag joins
  `VariantTag → variant → variantannotation` (latest `VariantAnnotationVersion` for the build)
  `→ transcript_version → gene_version → gene_symbol` — the same path
  `VariantTag.gene_symbol` (`analysis/models/models_variant_tag.py:70`) walks per-object, done as
  one ORM `values()/annotate()` aggregate.

### The cards

1. **Headline strip** — total tag events, distinct (variant, tag) pairs, distinct tagged variants,
   distinct analyses, active taggers this year. States the three identities once, up top, so every
   other card can just label which one it uses.
2. **Tags over time** — stacked bar per month, one series per tag (small tags grouped into
   "other"), with a pie of overall proportions beside it. Tag events.
3. **Your tagging** — the current user's counts: total, top tags, last-30-days sparkline, plus a
   note: *"Search for other users to see their stats on their user page"* (see user page below).
4. **By lab** — tag events rolled up per lab via each tagger's lab memberships; a user in multiple
   labs counts into each, footnoted on the card. Primary table is per-user (unambiguous), lab
   rollup secondary.
5. **Genes** — default: top 10 gene symbols as horizontal stacked bars, segments = that gene's top
   5 tags + "other". Below the chart, gene-symbol and tag autocompletes (both exist —
   `snpdb/autocomplete/Tag` and the genes app's symbol autocomplete) + Recalculate, so users can
   swap in their own gene/tag sets.
6. **Tag co-occurrence** — tag × tag heatmap counting distinct alleles carrying both tags (self-join
   `VariantTag` on `allele`, `tag_a < tag_b`). Conflict pairs (`Artefact`+`SomaticReportable`,
   `SomaticReportable`+`NOTSomaticReportable` — 62 and 55 variants in the export) glow here.
   Below it, a tag multi-select → link to the tagged-variants grid filtered to variants carrying
   **all** selected tags, with preset buttons for the known conflict pairs.
7. **Mega-artefacts** — top re-tagged `Artefact` variants with per-variant event counts and the
   aggregate ("these 40 variants account for N Artefact tags"), with copy pointing at the existing
   Tags node Exclude option for analysis templates. This is the #1751/#1805 bridge: it gives labs
   the number that tells them whether re-confirming each time is worth 4,599 clicks.
8. **A tag's genes over time** — pick a tag (default `SomaticReportable`) → stacked area per
   quarter of its top N genes + "other". The "how often do we report TP53/JAK2/KRAS" chart.

### Additions to existing pages

- **`variantopedia/variant_tags/`** (`variantopedia/views.py:479`): a small summary strip above the
  existing per-tag counts (totals + tags-this-month + link to the stats page).
- **Tagged-variants grid** (`TaggedVariantGrid`, `variantopedia/grids.py:292`): add a sortable
  tag-event count column (aggregate per variant), so "most re-tagged variants" is a one-click sort
  on the page people already use, rather than a separate widget.
- **User page** (`snpdb.views.view_user`, `snpdb/urls.py:102`): a "Tags & classifications" section —
  tag event count, top tags, classifications count. Same lazy-load endpoint the stats page's
  "Your tagging" card uses, parameterised by user.

## B. Tag merging

`Tag.id` is a CharField primary key (`snpdb/models/models.py:52`) — the name *is* the pk — so a
merge repoints foreign keys and deletes the row. The complete FK set:

| Model | Location | On collision after repoint |
|---|---|---|
| `VariantTag.tag` | `analysis/models/models_variant_tag.py:45` | delete rows that become exact dupes (same variant/tag/analysis/user) |
| `TagColor.tag` | `snpdb/models/models_user_settings.py:95` | keep the surviving tag's color, delete the other |
| `TagNodeTag.tag` | `analysis/models/nodes/filters/tag_node.py:141` | delete if node already references the survivor; bump affected nodes' versions so they re-run |
| `VCFTag.tag` | `snpdb/models/models_vcf.py:324` | delete dupes |
| `SampleTag.tag` | `snpdb/models/models_vcf.py:594` | delete dupes |

### Case policy

Keep display case (forcing lowercase would turn `NOTSomaticReportable` into mush) and **no DB
constraint** (case-insensitive uniqueness at the database level has been more trouble than it's
worth here). Instead, validate at the single creation choke point: `CreateTagForm` used by
`tag_settings` (`snpdb/views/views.py:1331`) rejects a name whose lowercase form matches an
existing tag, with a message naming the existing tag. Also replace the bare
`Tag.objects.create` + `except:` there with the form raising a proper validation error.

### Admin merge UI

Lives on `snpdb/settings/tags` (`snpdb/urls.py:84`), superuser-only.

- Each tag row/detail shows usage counts (per FK table) and **suggested merges**: tags at
  Levenshtein distance ≤ 1 case-insensitively (computed in Python — tag counts are tiny).
- **One flow only, started from the tag that will die**: "Merge `artefact` into → [tag picker]".
  There is deliberately no reverse "absorb other tags into this one" flow — a single mental model,
  always stated from the victim's side, is what stops an admin merging the wrong way.
- Confirmation states the direction in record terms and requires typing the dying tag's name
  (GitHub-style), since it cannot be undone:
  > All **404** variant tags, **2** analysis nodes and **3** color settings using `artefact` will
  > be changed to `Artefact`, and `artefact` will be deleted. This cannot be undone.
- The merge runs in a transaction: repoint each FK table per the collision rules above, bump
  affected `TagNode` versions, delete the tag, and report what changed.

### SA Path data migration

Per the manual-migration pattern (CLAUDE.md): a management command wrapping the same merge routine
the UI uses, plus a migration containing a `ManualOperation` with a `test=` callable so it only
registers on deployments that actually have case-colliding tags. For SA Path it merges
`artefact` → `Artefact` and, as part of the same pass, deletes true duplicate `VariantTag` rows —
same (variant, tag, analysis, user), keeping the earliest (~3,224 rows). Different-user re-tags are
kept: they are agreement data.

## Order of work

1. **B** first (merge + creation validation + migration) — small, self-contained, and the stats
   page then reports on clean data.
2. **A** stats endpoints + page, then the three existing-page additions (variant_tags strip, grid
   count column, user page section).
