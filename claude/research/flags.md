# flags — research notes

Verified against a96540a68 on 2026-09-25

The flags app hangs typed, commentable, resolvable "flags" off any model that mixes in
`flags/models/models.py:FlagsMixin` - today `snpdb/models/models_variant.py:Allele`,
`classification/models/classification.py:Classification` and
`classification/models/clinical_context_models.py:ClinicalContext`. It is small (one models module, one REST views
module, no tasks, commands or tests of its own), but classification leans on it for workflow state: submitted, unshared,
withdrawn, discordant, pending changes, significance change, condition resolution, internal review and suggestions are
all flag types (`classification/models/flag_types.py:ClassificationFlagTypes`). `flags/__flags_readme.md` records the
team's position: new development has moved away from flags, but existing uses stay. Fields and URLs are in the maps
([models](../maps/models.md#flags), [urls](../maps/urls.md#flags)); receivers are in [signals](../maps/signals.md).

## Flows

### Lookup tables, then collections

`FlagTypeContext` (classification, allele, clinical_context), `FlagType`, `FlagResolution` and the
`FlagTypeResolution` join that says which resolutions a type allows are reference data created by data migrations -
`flags/migrations/0002_initial_data.py` for the bulk, later ones (and some in classification's migrations) for types
added since. A resolution carries a `flags/models/enums.py:FlagStatus` (open, closed, rejected); a flag is "open"
when its resolution's status is OPEN (`FlagCollection.Q_OPEN_FLAGS`). `flags/models/models.py:FlagType.default_resolution`
picks the type's OPEN resolution, falling back to CLOSED. Each consuming app wraps the type ids in a helper object
(`classification/models/flag_types.py:ClassificationFlagTypes`, `snpdb/models/flag_types.py:AlleleFlagTypes`, whose
allele types are all commented out now).

A model gets a `flags/models/models.py:FlagCollection` lazily: `FlagsMixin.flag_collection_safe` creates one in the
model's `flag_type_context()` and saves the FK. The collection is the unit the UI and API address; it has no FK back to
its owner, so `flags/models/models.py:FlagCollection.source_object` finds it either from the extra-info signal (below)
or by probing every reverse `*_set` accessor on the collection until one returns a row.

### Raising, commenting, resolving

Everything goes through `flags/models/models.py:Flag.flag_action`, which is atomic: it optionally changes the resolution
(checked against `FlagType.permission`), writes a `FlagComment` (whose `resolution` is set only when this comment
changed it) and sends `flag_comment_action`. `flags/models/models.py:FlagCollection.add_flag` creates the Flag already
at its resolution and calls `flag_action(first_comment=True)` so every flag starts with a comment recording its opening
resolution - `flags/models/flag_health_check.py:flag_chanced_since` depends on that to reconstruct a flag's state at a
past date. No `user` means `admin_bot` with permission checks off.

Server code rarely calls `add_flag` directly. The idempotent entry points are
`flags/models/models.py:FlagCollection.get_or_create_open_flag_of_type` (reuse the open flag; reopen a closed one when
`reopen`, the type is `only_one`, or `reopen_if_bot_closed` and the last comment was admin_bot's; match or close others
by `data` keys via `old_data` / `close_other_data`), `FlagCollection.ensure_resolution` (only_one types only: move the
latest flag to a resolution, creating it only if that resolution is open - used by classification import, significance
change and ClinVar exclusion) and `close_open_flags_of_type` (bot-closes every open flag of a type, optionally matching
`data`). `FlagCollection.filter_for_flags` (and its deprecated alias `filter_for_open_flags`, still the one callers use)
turns "has an open flag of these types" into a subquery on any FlagsMixin queryset - the classification dashboard and
`snpdb/models/models_user_settings.py:UserSettings.classification_issue_count` are built from it.

### Permissions

`flags/models/models.py:FlagPermissionLevel` is a totally ordered str enum: NO_PERM < USERS < OWNER < ADMIN < SYSTEM.
`FlagCollection.permission_level` gives admin_bot SYSTEM and superusers ADMIN, and otherwise asks the source object's
`flag_user_permission` - OWNER if it `can_write`, USERS if it `can_view`
(`classification/models/classification.py:Classification.flag_user_permission` checks view on the last published
version). A type's `raise_permission` gates raising, `permission` gates changing resolution; `user_private` flags are
hidden from USERS other than the raiser (`flags/views/views.py:FlagHelper.is_viewable_flag`).

### The UI and the REST API

All flag UI is `variantgrid/static_files/default_static/js/flags.js`, a client-side store fed by two DRF views.
`flags/views/views.py:FlagsView` GET takes comma-separated collection ids and returns flag types, resolutions,
collections (with the viewer's `user_permission` level and the extra info), open and just-closed flags; `?history=<id>`
returns one collection's full comment history and `?since=<ts>` the comments since a poll. POST raises a flag.
`flags/views/views.py:FlagView` returns or acts on one flag. Serialisation is `flags/views/views.py:FlagHelper`, which
calls `flags/models/models.py:fetch_flag_infos`: it sends `flag_collection_extra_info_signal` with a
`flags/models/models.py:FlagInfos`, and the owning apps' receivers
(`classification/models/classification.py:get_extra_info`, `snpdb/models/models_variant.py:get_extra_info`) attach a
label, links and the `source_object` for each collection. ClinicalContext has no receiver, so its collections fall back
to the `*_set` probe and a generic label.

### Reacting to flags

`flag_comment_action` is how other apps follow flag changes: `classification/models/condition_text_matching.py:check_for_withdrawn`
resyncs condition text matching when the withdrawn flag moves (the receiver in `classification/models/clinvar_export_models.py`
is a stub). `flags/models/flag_health_check.py:flag_chanced_since` counts flags opened and resolved since a date per type
for `classification/signals/classification_health_checks.py`.

## Why it is shaped this way

A generic collection per object, rather than per-model flag tables, means a new workflow state is a data migration
rather than a schema change, and gave every object the same comment thread and dialog in
`flags.js`. The cost is that flags carry little structure - hence `data` (JSON, added so flag events carry the value
that triggered them, `bb69fffc7`) and `FlagType.attributes` (queried by `FlagsMixin.has_open_flag_with_attribute` and
`FlagCollection.get_flag_of_type`). That ceiling is why the readme says new features should not be built on flags.

`FlagType` uses `library/django_utils/django_object_managers.py:ObjectManagerCachingImmutable` because types are looked up
by id all over classification code and practically never change. The collection deliberately knows nothing about its
owner so the app has no dependency on classification or snpdb; the signal is the seam in the other direction.

## History

The app dates from 2020 as part of the classification (then "variant classification") work;
`flags/migrations/0003_rename_variant_classification.py` renamed the context and type ids when the model was renamed.
Later data migrations added the "not public" ClinVar type (0005), pending changes (0008) and condition resolution
(0009); `ensure_resolution`-style upkeep of only_one types followed. Watching / unseen activity was never finished: `FlagWatch` and `set_watcher` were removed in `486ec74f4` (2022,
which also added the admin raise permission) but the view still calls it (Traps), sub-flags survive only as
`FlagInfos.record_sub_flag` with the caller commented out, and allele flags (missing 37/38, `0004`) are no longer
raised. Recent changes are housekeeping: the DRF endpoints lost format suffixes (#1375), the flag event data landed
(2025), and the 2026 sweeps (#1795 dead code, #1773 import cycles, #1421 `library/django_utils` move) touched imports only.

## Traps

**The REST API trusts the browser for permissions (confirmed, not fixed).** `flags/views/views.py:FlagHelper.add_flag`
calls `FlagCollection.add_flag` without `permission_check`, whose default is False, so `FlagsView` POST lets any
logged-in user raise any flag type - SYSTEM-only ones included, even from another context - on any collection id;
`raisableFlags()` in `flags.js` is the only gate. `FlagsView` GET and `?history=` likewise return non-private flags
and comments of any collection to a user whose level is NO_PERM, and `Flag.flag_action` lets a NO_PERM user comment
on a flag (including a private one: its private check only looks at USERS).

**`add_flag(permission_check=True)` would crash (confirmed, latent).** It compares
`current_level < flag_type.raise_permission`, a plain str; `FlagPermissionLevel.__lt__` reads `other.level`, so it
raises AttributeError. Every caller passes `permission_check=False` or no user today, which is why it has not been seen;
compare against `flag_type.raise_permission_enum`.

**`NO_PERM` is truthy (confirmed, latent).** `FlagPermissionLevel.NO_PERM` is `'0'`, a non-empty string, so the
`if not permission_level` guard in `flags/models/models.py:FlagCollection.flags` never fires and a NO_PERM user gets
every flag, private ones too (only USERS is filtered). No caller passes a user at present. Compare with
`== FlagPermissionLevel.NO_PERM` or `.level`.

**A `watch` POST 500s (confirmed, dead path).** `FlagsView` POST with `watch` calls `fc.set_watcher`, removed in 2022;
`watch_toggle()` in `flags.js` still exists but nothing calls it. Delete both rather than revive them.

**`classification_not_public` is disabled by accident of JavaScript.** Without `CLINVAR_EXPORT` mode,
`FlagHelper.to_json` sets that type's `permission` / `raise_permission` to `FlagPermissionLevel.SYSTEM` (the string `"A"`)
instead of its `.level` (4); `flags.js` compares numbers with `>=`, and `n >= "A"` is always false, so it works.

**FlagType edits need a restart.** `ObjectManagerCachingImmutable` caches `get()` per process (not under unit tests), so
changing a type in admin or a shell is invisible to running gunicorn/celery workers until they restart. Add types by
data migration.

**Every flag needs its opening comment.** A comment's `resolution` is null unless it changed the resolution, so the
state at a date is "the latest comment with a resolution before it" (`flag_chanced_since`); creating a Flag without
going through `add_flag` leaves no opening comment and breaks that. `get_or_create_open_flag_of_type` with
`reopen_if_bot_closed` would also fail on such a flag (`FlagComment.last` returns None).

**`source_object` can be None.** An orphaned collection (owner deleted - the FK is on the owner, CASCADE runs the other
way) logs a warning and yields NO_PERM for everyone but superusers and admin_bot.
