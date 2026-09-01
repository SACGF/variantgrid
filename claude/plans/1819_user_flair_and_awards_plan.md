# User flair and awards (crowns, trophies, badges) #1819

Written by Claude Fable 5 (claude-fable-5), 2026-09-01


## Design

Three things, kept deliberately separate:

1. **Flair** — every user picks an emoji (🐱 etc.) on their settings page. It shows on their own
   user page and settings page only. It is self-expression, not a prize.
2. **Awards** — three kinds:
   - **Title**: held until someone takes it. 👑 all-time, 🥇 this month, 🏆 today. Computed
     nightly/hourly; loses `active` when lost.
   - **Badge**: permanent once earned, tiered bronze/silver/gold like GitHub achievements
     (Cold Case, Night Owl, milestones…). Never revoked. Profile shows progress to the next tier.
   - **Kudos**: the existing hand-given `UserAward`. Permanent.
3. **Grid rendering** — names on grids and in `{% user %}` are plain, with one exception: while a
   user holds a title, grids render `👑 Evelyn 🐱` — title icon before the name, *their own flair*
   after it. Badges and kudos never appear on grids. Rarity is the point: one crown and one
   trophy per subject across the whole site.

Off switches:
- `settings.USER_AWARDS_ENABLED` (deployment-wide, default `True`; clinical deployments can turn it
  off — no computation, no rendering, settings UI hidden).
- `SettingsOverride.show_user_awards` — the global/org/lab/user cascade. Controls whether the
  *viewer* sees title decorations on grids. Flair picker and the award cabinet on the profile stay.

Awards reward chores the system needs done (Cold Case, condition matching, wiki) as much as
volume. Classification-count titles are included but every deployment can disable individual
definitions via `settings.USER_AWARDS_DISABLED_KEYS`.

## What already exists

- `UserAward` (`snpdb/models/models.py:872`): `user`, `award_text`, `award_level`
  (`UserAwardLevel` G/S/B, `snpdb/models/models_enums.py:11`), `active`. Admin-created only
  (`snpdb/admin.py:336`).
- `UserAwards` helper (`snpdb/models/models.py:890`) sorts them and renders `html` — a single
  trophy icon overlaid on the avatar, tooltip listing all awards.
- Rendered in three places: `snpdb/templates/snpdb/tags/user.html:7` (avatar overlay in the
  `{% user %}` tag), `uicore/templates/uicore/page/base.html:201` (navbar), and the list on
  `snpdb/templates/snpdb/settings/view_user_settings.html:75`.
- CSS `.user-award`, `.user-award-gold/silver/bronze`
  (`variantgrid/static_files/default_static/css/global.scss:260-290`).
- `AvatarDetails` (`snpdb/models/models_user_settings.py:399`) — `preferred_label`, `awards`,
  `url`, `background_color`. Single choke point for "how a user is displayed".
- `UserSettings` cascade: model field on `SettingsOverride`, dataclass field on `UserSettings`,
  widget in `SettingsOverrideForm` (`snpdb/forms.py:523`), migration. `variant_tag_stale_days`
  (commit `git log -S variant_tag_stale_days`) is the template to copy.
- 16 grids render users as raw `RichColumn(key="user__username")`: `analysis/grids.py:495,553,716,
  867(actor),927,1091`, `snpdb/grids.py:82,146,352,391,432,455,752`, `pathtests/grids.py:25`,
  `pedigree/grids.py:17,37`. These bypass `AvatarDetails` entirely.
- Per-user tag stats already computed for the user page (`tag_stats_for_other_user`,
  `variantopedia/views.py`) — reuse its query shape for the tagger titles.
- `auditlog` is registered on `Analysis` and node sub-models (`analysis/models/models_analysis.py:804`,
  `analysis/models/nodes/analysis_node.py:1655`) — `LogEntry.actor`/`timestamp` gives "analyses
  modified by user" without new tracking.
- Beat schedule lives in `variantgrid/celery.py` (`app.conf.beat_schedule`); scheduled tasks run on
  `scheduling_single_worker`. Note the comment there: `crontab()` has had timezone trouble, raw
  seconds are reliable.
- User search already exists (`snpdb/signals/user_search.py`) → `view_user`
  (`snpdb/views/views.py:1087`, template `snpdb/templates/snpdb/settings/view_user.html`).

## Part A — Model changes (`snpdb`)

### `UserAward`

Add:

```python
kind = models.CharField(max_length=1, choices=UserAwardKind.choices, default=UserAwardKind.KUDOS)
definition_key = models.TextField(null=True, blank=True)   # null for kudos
subject = models.TextField(null=True, blank=True)          # e.g. tag name for per-tag titles
period = models.CharField(max_length=1, choices=AwardPeriod.choices, null=True, blank=True)  # titles only
count = models.IntegerField(null=True, blank=True)         # score (titles) / raw progress (badges)
```

- `UserAwardKind`: `TITLE = "T"`, `BADGE = "B"`, `KUDOS = "K"`. `AwardPeriod`: `ALL_TIME = "A"`,
  `MONTH = "M"`, `DAY = "D"`.
- Unique constraint on `(user, definition_key, subject, period)` where `definition_key` is not null
  (`UniqueConstraint(condition=Q(definition_key__isnull=False))`).
- Computed rows are reused, not recreated: a title regained flips `active` back on; `created`
  is first-earned, `modified` is last change. `award_text` is filled from the definition so the
  existing `__str__`/tooltip keep working.
- `icon_class` property: titles → `fa-crown` / `fa-medal` / `fa-trophy` by period; badges →
  the definition's icon; kudos → `fa-trophy` as now. Level colour classes unchanged.
- Existing rows: migration sets `kind=KUDOS` (the default) — no data migration needed.

### `UserSettingsOverride.flair`

```python
flair = models.CharField(max_length=16, null=True, blank=True, choices=USER_FLAIR_CHOICES)
```

`USER_FLAIR_CHOICES` is a curated list (~60 emoji: animals, food, space to match the avatar
theme, a few objects) in `snpdb/models/models_enums.py`. A curated list avoids grapheme-cluster
validation and guarantees consistent rendering. Field lives on `UserSettingsOverride` only (not the
cascade — it is personal). Form: `RadioSelect` with a `flair-picker` CSS class laying the radios
out as a grid of buttons (`global.scss`).

Admins can set another user's flair too: `UserSettingsOverride` is already in Django admin
(`snpdb/admin.py:31`) so the field appears there for free, and `view_user` gets a small
superuser-only flair form (same picker widget, POSTs to `view_user`) so an admin can hand someone
a 🐐 without leaving their page. "Flair" is the name used everywhere — model field, form label,
CSS, docs.

### `SettingsOverride.show_user_awards`

`BooleanField(null=True, blank=True)`, cascade field like the others: add to `SettingsOverride`,
the `UserSettings` dataclass, `BlankNullBooleanSelect` widget in `SettingsOverrideForm`, and the
`GlobalSettings` default (`True`). Migration `snpdb/migrations/0231_*.py`.

### `settings`

`default_settings.py`: `USER_AWARDS_ENABLED = True`, `USER_AWARDS_DISABLED_KEYS = set()`.

## Part B — Award definitions registry (`snpdb/user_awards.py`)

```python
@dataclass(frozen=True)
class AwardDefinition:
    key: str                         # "top_tagger", "cold_case"
    kind: UserAwardKind
    title: str                       # "Top tagger"
    description: str                 # how to earn it - shown on locked badges
    icon: str                        # font-awesome class, e.g. "fa-tag"; profile + subject icon on grids
    periods: tuple[AwardPeriod, ...] = ()   # titles only
    tiers: tuple[int, ...] = ()             # badges only: (bronze, silver, gold) thresholds

    def compute(self, since: Optional[datetime]) -> dict[Optional[str], dict[int, int]]:
        """ {subject: {user_id: count}} - subject None for un-subjected awards.
            since=None means all time. Badges ignore since. """
```

- `register_award(definition)` / `get_award_definitions()` — module-level registry. Each app
  defines its awards in `<app>/user_awards.py` and imports the module from `AppConfig.ready()`
  (the same pattern the apps use to load signal receivers).
- `admin_bot()` is always excluded from counts.
- Definitions listed in `settings.USER_AWARDS_DISABLED_KEYS` are skipped by the task and their
  existing awards deactivated.

### Title computation

For each definition, period and subject: take the user(s) with the max count (count ≥ 1; ties
all hold it). Upsert `UserAward(kind=TITLE, active=True, count=…)`; deactivate every other active
row for that `(definition_key, subject, period)`. `award_text` e.g.
`Top tagger of 'Pathogenic' (all time)`.

### Badge computation

For each definition and user: `count` → tier = highest threshold ≤ count. Upsert with
`award_level`=tier and `count`. Rows with count below bronze are still upserted with
`active=False` so the profile can show progress ("Cold Case: 3 / 5 for bronze"). Badges are never
deactivated once earned; `award_level` only ever moves up.

## Part C — Scheduled recompute

- `snpdb/tasks/user_award_tasks.py`: `@shared_task(queue=SCHEDULING_SINGLE_WORKER)
  update_user_awards(periods: list[str])`. Exits immediately when `USER_AWARDS_ENABLED` is off.
- Beat entries in `variantgrid/celery.py` (raw seconds, per the timezone note there):
  - `user-awards-daily`: `[DAY]` every hour — the trophy for "today" should move during the day.
  - `user-awards-nightly`: `[ALL_TIME, MONTH]` + all badges, every 24h.
- `manage.py update_user_awards` runs the same thing synchronously, with a `ManualOperation`
  migration so existing deployments get a backfill on upgrade.
- Post-save signals are deliberately not used: the queries are cheap in aggregate and a title
  changing hands mid-day via signals would be a distraction, not a feature.

## Part D — Rendering

### `AvatarDetails` (`snpdb/models/models_user_settings.py`)

- `flair: Optional[str]` — from `UserSettingsOverride.flair`.
- `titles: list[UserAward]` — active TITLE awards, ordered ALL_TIME → MONTH → DAY.
- `title_icon_html` — the highest-period title's icon (crown beats medal beats trophy) with a
  tooltip listing every held title. Empty when none.
- `grid_label_html(viewer_settings: UserSettings) -> SafeString` — the single rule:

  ```
  if settings.USER_AWARDS_ENABLED and viewer_settings.show_user_awards and self.titles:
      f"{title_icon} {escape(preferred_label)} {flair or ''}"
  else:
      escape(preferred_label)
  ```

- `awards` (existing `UserAwards`) gains `badges`, `kudos`, `titles` partitions and a
  `progress(definition)` helper for the profile page.

### Grids

- New factory in `snpdb/views/datatable_view.py`: `user_column(key="user", label="Created by",
  **kwargs) -> RichColumn` — sorts/filters on `<key>__username` (as now) but renders via a server
  `renderer` that looks up `AvatarDetails` for `<key>__id` (add to `extra_columns`) and returns
  `grid_label_html(self.viewer_settings)`. `AvatarDetails` per user id is cached on the
  `DatatableConfig` instance for the request so a 100-row grid does one settings lookup per user.
- Replace the 16 `RichColumn(key="…__username")` sites with `user_column(...)`. `actor__username`
  at `analysis/grids.py:867` is the same shape with `key="actor"`.
- jqgrid variant grids (`library/jqgrid/`) show tag usernames inside the tag cell; leave them —
  the tag hover already shows the user via `{% user %}`.

### `{% user %}` tag and navbar

- `user.html:7` and `base.html:201`: swap `awards.html` for `title_icon_html` — the avatar overlay
  becomes titles-only, so a kudos or badge stops decorating the user everywhere. Add the flair
  after the name when `title_icon_html` is set, mirroring the grid rule. Respect the viewer's
  `show_user_awards` via `user_details.user_settings` (already cached in the context).

### Profile pages

`view_user.html` and `view_user_settings.html` share a new include `snpdb/tags/award_cabinet.html`
(inclusion tag `{% award_cabinet user %}` in `user_tags.py`):

- Name rendered with flair (always, regardless of titles — this is the one place flair is
  unconditional).
- **Titles held**: icon, text, count, since `modified`.
- **Badges**: earned ones with tier colour, count, earned date (`created`); unearned ones greyed
  with description and a progress bar to the next tier (`count / next_threshold`). Hidden badges
  (fun ones like Night Owl) are listed only once earned — `AwardDefinition.hidden: bool`.
- **Kudos**: the existing list.
- Settings page additionally has the flair picker (`UserSettingsOverrideForm.flair`) and the
  `show_user_awards` cascade row.
- Everything under `{% if settings.USER_AWARDS_ENABLED %}` — expose via `settings_tags`.

## Part E — Initial award definitions

### Titles (`periods = (ALL_TIME, MONTH, DAY)` unless noted)

| key | app | subject | count | icon |
|---|---|---|---|---|
| `top_tagger` | analysis | per `Tag` and `None` (overall) | `VariantTag` created by user | `fa-tag` |
| `top_classifier` | classification | None | published `ClassificationModification` by user, distinct classifications | reuse the classification menu icon |
| `top_analyst` | analysis | None | distinct `Analysis` with an `auditlog.LogEntry` by actor (Analysis + node content types) | `fa-project-diagram` |
| `top_wiki` | genes/variantopedia | None | `Wiki` subclasses `last_edited_by=user`, `modified` in period (no edit history exists — approximation, documented) | `fa-book` |

### Badges (`tiers = (bronze, silver, gold)`)

| key | app | count | tiers | icon | hidden |
|---|---|---|---|---|---|
| `cold_case` | classification | classifications the user modified where the previous modification was > 2 years older | 1/5/25 | `fa-user-secret` | no |
| `matchmaker` | classification | `ConditionTextMatch.last_edited_by=user` with non-empty `condition_xrefs` | 10/100/1000 | `fa-link` | no |
| `wiki_scribe` | genes | wikis `last_edited_by=user` | 5/25/100 | `fa-feather` | no |
| `tagger` | analysis | tags created | 100/1000/10000 | `fa-tags` | no |
| `classifier` | classification | classifications created | 10/100/1000 | classification icon | no |
| `analyst` | analysis | analyses created (`template_type` null) | 10/100/1000 | `fa-project-diagram` | no |
| `founder` | snpdb | years since `date_joined` | 1/3/5 | `fa-landmark` | no |
| `night_owl` | analysis | tags created 22:00–05:00 in the user's `UserSettings.tz` | 10/100/1000 | `fa-moon` | yes |
| `early_bird` | analysis | tags created 05:00–07:00 in the user's tz | 10/100/1000 | `fa-sun` | yes |

Deferred until the data exists to compute them (note in follow-ups, not built now):
- `peacemaker` (discordances resolved) — `DiscordanceReport` records `resolution` but the
  implementer must confirm who resolved it is recorded; otherwise attribute via the concordant
  `ClassificationModification`.
- `janitor` (stale tags cleaned) — `VariantTag` deletion leaves no record.
- `fixer` (`ImportedAlleleInfoValidation` issues resolved) — needs a resolved-by field.

## Part F — Tests

- `snpdb/tests/test_user_awards.py`: title upsert/deactivate on change of holder, ties, badge
  tier only moves up, count-below-bronze row inactive, disabled keys deactivated,
  `USER_AWARDS_ENABLED=False` short-circuits.
- `grid_label_html`: plain when no titles; decorated with flair when held; plain when viewer has
  `show_user_awards=False`; flair-less title when user has no flair.
- One definition test per app for the query (e.g. `cold_case` picks up the > 2 year gap and
  ignores a 1 year gap; `night_owl` respects the user's timezone).
- `URLTestCase` entries for the settings page with the picker and a user page with awards.

## Phases

1. Model + migration + settings flags + `AvatarDetails` additions + registry + task + management
   command (no definitions yet). Existing kudos rendering keeps working.
2. Grid `user_column` and the 16 replacements; `user.html`/`base.html` titles-only overlay.
3. Flair picker, `show_user_awards` cascade row, award cabinet on both profile pages.
4. Titles: `top_tagger`, `top_classifier`, `top_analyst`, `top_wiki`.
5. Badges: `cold_case`, `matchmaker`, `wiki_scribe`, milestones, `founder`, `night_owl`,
   `early_bird`.

## Follow-ups (out of scope)

- Lab-head kudos UI (give an award from the user page) — the model already supports it.
- Lab-level ratio awards (lowest stale-tag %, highest evidence completeness) on the lab page.
- Weekly email line ("3 tags behind X for the Pathogenic crown") via the existing
  `email_weekly_updates`.
- Deferred badges above once their events are recorded.
