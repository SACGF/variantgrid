# Integration status — "when did each external system last check in?"

Bring Mocha's dashboard "Sync" cards to VariantGrid, on the Server Status page, with a
registration mechanism so a deployment-specific app (SA Path) contributes its own cards without
core ever importing it.

## What Mocha does today

`core/templates/home.html` renders a row of cards under a **Sync** heading, one per external
system, each showing a few named timestamps and (optionally) a manual re-sync button:

| Card | Detail lines | Manual sync |
|---|---|---|
| Helix | Last Run, Last Modified | `import_ngs_helix.delay()` |
| Omico | Last Action, Last Success, Last Fail | — |
| MO Spreadsheet | Last Import | `import_mo_legacy_spreadsheet.delay()` |

`core/views.py::HomeView` builds them by hand — one `_get_*_sync()` method per system, each
returning `{"manual_sync": <post key or False>, "sync_details": {label: datetime}}`, all
hard-coded into one view. That shape is right; the hard-coding is what we replace.

Every card's data comes from a "when did this run" record that already exists somewhere —
`HelixNGSOrderSync` rows, `Action` rows against `ExternalSystem.OMICO_API`, `SpreadsheetImport`
rows. Mocha's Omico card is the interesting one: it shows last attempt / last success / last
failure separately, which is what tells you an integration is *running but broken* rather than
*not running*.

## What VariantGrid has to build on

`library/health_check.py` already has a signal-driven, app-pluggable status system —
`health_check_signal` / `health_check_overall_stats_signal`, receivers registered from each app's
`ready()` (`sync/apps.py`, `snpdb/apps.py`, `eventlog/apps.py`), collected with `send_robust` so
one broken receiver reports itself instead of taking the page down. `HealthCheckAge`
(`sync/signals/sync_health_check.py`) is already reporting "SYNC (Shariant) — 0 days old" into
the nightly Slack digest, and `variantopedia/views.py::health_check_details` already renders
those same receivers to HTML on a Server Status tab.

So the registration mechanism the user asked for exists and is proven — this plan follows its
shape rather than inventing a second one.

Health checks stay what they are: a nightly Slack digest, day-granularity, mixing record counts
with ages. Integrations get their own signal because the questions differ — an integration card
wants attempt/success/failure as three separate timestamps at minute granularity, a link to the
records it wrote, and a button to run it now. §6 wires the two together so a stale integration
still reaches Slack.

## The design

Five pieces. §1–§4 are core and land together; §5–§8 build on them.

---

### 1. The registry — `library/integration_status.py`

New module alongside `library/health_check.py`, same house style (dataclasses + a
`django.dispatch.Signal` + a collector that uses `send_robust`).

```python
@dataclass
class IntegrationDetail:
    """ One "Last Run: 4 minutes ago" item inside an integration's box """
    label: str
    timestamp: Optional[datetime]
    status: str = "info"          # bootstrap contextual class - info/success/warning/danger
    extra: Optional[str] = None   # eg the error message against "Last Failure"


@dataclass
class IntegrationTrigger:
    """ A "Sync now" button. The callable is held here rather than named in the POST, so the
        page can only ever run something an app registered """
    action_id: str                # unique across the deployment, eg "sapath-helix-sync"
    run: Callable[[], Any]        # usually some_task.delay
    label: str = "Sync now"


@dataclass
class IntegrationStatus:
    name: str                                # "Helix NGS Database"
    details: list[IntegrationDetail]
    direction: IntegrationDirection = IntegrationDirection.INBOUND
    description: Optional[str] = None        # one line under the title
    url: Optional[str] = None                # link to the records this integration writes
    record_count: Optional[int] = None
    trigger: Optional[IntegrationTrigger] = None
    warning_age: Optional[timedelta] = None  # feeds §6
    enabled: bool = True                     # renders greyed with a "disabled" badge
    sort_order: int = 0


class IntegrationDirection(models.TextChoices):
    INBOUND = "I", "Inbound"     # something outside pushes/we pull into VG
    OUTBOUND = "O", "Outbound"   # VG calls out


integration_status_signal = django.dispatch.Signal()


def get_integration_statuses() -> tuple[list[IntegrationStatus], list[str]]:
    """ Collected with send_robust: returns the statuses plus one message per receiver that
        raised, so a broken provider shows as a message on the page instead of a 500 """


def run_integration_trigger(action_id: str) -> Optional[IntegrationStatus]:
    """ Re-collects the statuses, finds the one whose trigger matches, calls run() """
```

`IntegrationStatus.last_activity` — a property returning the newest non-null timestamp across
`details` — gives §6 and the default sort something to work from without every provider
repeating itself.

Receivers return an `IntegrationStatus`, a list of them, or `None` (nothing configured on this
deployment — a Shariant box has no Helix row). `get_integration_statuses` flattens and drops
the `None`s, exactly as `populate_health_check` does.

---

### 2. Recording model — `eventlog.models.IntegrationActivity`

Several integrations have no natural "when did this run" record. Mocha's outbound API import
(`sapath/tasks/import_mocha.py`) writes `MochaSampleExtraction` rows but keeps no run log, so
"last run" is currently indistinguishable from "last row that happened to change". Incoming
VariantGrid API calls have no record at all — `eventlog`'s `PageViewsMiddleware` explicitly skips
them (`IGNORE_SEGMENTS = {"api", ...}`).

Rather than every integration growing its own `*Sync` model, core gets one:

```python
class IntegrationActivityStatus(models.TextChoices):
    SUCCESS = "S", "Success"
    NO_CHANGE = "N", "No change"
    ERROR = "E", "Error"


class IntegrationActivity(TimeStampedModel):
    """ One row per integration, updated in place - a rolling "last seen" for external systems
        that keep no run log of their own. Bounded however chatty the integration is """
    key = models.TextField(unique=True)          # "sapath-mocha-api", "api:specimen_measure_bulk_create"
    name = models.TextField()                    # "Mocha API" - the label on the Server Status row
    direction = models.CharField(max_length=1, choices=IntegrationDirection.choices)
    last_attempt = models.DateTimeField(null=True)
    last_success = models.DateTimeField(null=True)
    last_change = models.DateTimeField(null=True)   # last attempt that actually wrote something
    last_error = models.DateTimeField(null=True)
    last_error_message = models.TextField(null=True, blank=True)
    attempt_count = models.IntegerField(default=0)
    success_count = models.IntegerField(default=0)
    error_count = models.IntegerField(default=0)
```

`eventlog` is the home: it already owns audit/activity models, has migrations, and is
unconditionally installed. New migration `eventlog/migrations/0004_integrationactivity.py` —
create-table only, so no `ManualOperation` is needed.

Two ways to write:

```python
@classmethod
def record(cls, key, name=None, direction=..., status=..., changed=False, message=None):
    """ Single upsert. Counters go through F() so concurrent workers add up """

@classmethod
@contextmanager
def track(cls, key, name=None, direction=..., ):
    """
    with IntegrationActivity.track("sapath-mocha-api", name="Mocha API") as activity:
        activity.changed = bool(import_mocha_sample_extractions())
    """
```

`track` stamps `last_attempt` on entry, `last_success` (and `last_change` when the body set
`changed`) on clean exit, and on exception stamps `last_error` + `last_error_message` and
re-raises — the caller's error handling and Rollbar reporting are untouched.

---

### 3. Rows from recorded activity — `eventlog/signals/integration_activity_status.py`

One shipped receiver turns every `IntegrationActivity` row into an `IntegrationStatus`:

```python
@receiver(signal=integration_status_signal)
def integration_activity_status(sender, **kwargs) -> list[IntegrationStatus]:
    ...  # Last Run / Last Success / Last Change / Last Failure(+message) per row
```

`Last Failure` renders `status="danger"` when it is newer than `last_success`, `"warning"`
otherwise — so a system that failed once yesterday and has worked since reads differently from
one that is failing now. This is Mocha's Omico card, generalised.

Registered from `EventlogConfig.ready()` next to the existing `active_users_health_check` import.
The payoff: **an integration that calls `IntegrationActivity.track(...)` appears on Server Status
with no registration at all**, and apps only write a receiver when they want something richer.

---

### 4. The page — an "Integrations" section on Server Status

The section sits **directly below Celery Workers**, inside the existing embedded Server Status
tab in `variantopedia/templates/variantopedia/server_status.html`, and is built out of the same
parts that section is:

```
{% if integrations %}
    <h4>Integrations</h4>
    {% for integration in integrations %}
        {% labelled id=integration.slug label=integration.name %}
            <div class="alert {% if integration.ok %}alert-primary{% else %}alert-danger{% endif %}">
                ...
            </div>
        {% endlabelled %}
    {% endfor %}
{% endif %}
```

So each integration is one `{% labelled %}` row — label on the left, a full-width contextual
`alert` box on the right — matching the Celery worker boxes in size, spacing and colour. The
alert reads `alert-primary` when the integration's most recent outcome was a success and
`alert-danger` when it was a failure, the same `worker_info.ok` convention.

Inside the box, the detail lines run on one line the way a worker's `Jobs: / Tasks: / Scheduled:`
run, each timestamp through the existing `{% timestamp ts time_ago=True %}` tag
(`uicore/templatetags/js_tags.py:227`) so it reads "4 minutes ago" with the absolute time in the
user's timezone on hover:

> **Helix NGS Database**   Last Run: 4 minutes ago   Last Change: 6 hours ago   Orders: 18,204

A failing integration adds `Last Failure: 20 minutes ago` and its message. A disabled one renders
`alert-secondary` with a `disabled` badge. An integration that has never run says so in place of
its timestamps, rather than being hidden.

`IntegrationStatus.slug` — `slugify(name)` — gives `{% labelled %}` the stable `id` it wants,
the way the worker name does today.

**Wiring**: `variantopedia/views.py::server_status` calls `get_integration_statuses()` and puts
`integrations` (and any collector messages) into the context it already builds — no new view, no
new URL, no ajax tab. The whole section is inside `{% if integrations %}`, so deployments with
nothing registered see the page exactly as it is now.

Building the section out of the shared `IntegrationStatus` dataclasses (rather than shaping the
data in the view) is what lets it move to a landing page or its own URL later — which is where
the "put it on the server status page **to start**" is heading — without the providers changing.

---

### 5. Manual "Sync now" buttons

Each integration with a `trigger` renders a small `btn-warning` submit inside its alert box, in a
`<form method="POST">` carrying a hidden `action_id` and `action=run-integration` — the same
plain post-to-the-current-page the `Test Slack` / `Health Check` buttons use
(`server_status_settings_detail.html:11`). The existing action dispatch in
`variantopedia/views.py::server_status` gains one branch that calls
`run_integration_trigger(action_id)` and adds a `messages.info` naming what it started, matching
how `Test Slack` / `Health Check` already behave. The page is `@require_superuser`, and the
`action_id` only resolves against triggers an app registered, so nothing arbitrary is reachable.

---

### 6. Stale integrations reach Slack

One bridging receiver in `library/health_check.py`'s consumers — a new
`variantopedia/signals/integration_health_check.py` — connects to
`health_check_overall_stats_signal`, walks `get_integration_statuses()`, and emits a
`HealthCheckAge(name=..., last_performed=status.last_activity, warning_age=status.warning_age)`
for every enabled integration that declared a `warning_age`.

Every integration is then visible on the page live *and* in the nightly digest when it goes
quiet, from a single registration. `sync/signals/sync_health_check.py` keeps its own
`HealthCheckAge` — §7 gives sync a card that reads the same `SyncRun` rows, and the two stay
independent.

---

### 7. Core providers shipped with the framework

Each is a receiver module registered from its app's `ready()`, following
`sync/signals/sync_health_check.py` line for line.

| Provider | File | Reads | Direction | Shows |
|---|---|---|---|---|
| Sync destinations | `sync/signals/sync_integration_status.py` | `SyncDestination` + latest `SyncRun` per destination | per destination | Last Run / Last Success / Last Failure, `enabled=destination.enabled`, links to the sync destination page |
| Classification imports | `classification/signals/classification_import_integration_status.py` | `ClassificationImportRun` | Inbound | Last Import, ongoing-import count, `warning_age` left unset (labs import on their own cadence) |
| ClinVar export | `classification/signals/clinvar_export_integration_status.py` | `ClinVarExportBatch` (latest `SUBMITTED`, latest `REJECTED`) per `ClinVarKey` | Outbound | Last Submitted / Last Rejected, links to the batch |
| Recorded activity | `eventlog/signals/integration_activity_status.py` | `IntegrationActivity` (§3) | per row | §3 |

Sequencing (`SequencingRun` last seen, behind `settings.SEQAUTO_ENABLED`) and annotation
(`VariantAnnotationVersion` / ClinVar annotation dates) are the obvious next two and slot in the
same way once the first four are on a page and the shape has been eyeballed.

---

### 8. SA Path registration (`../variantgrid_sapath`, separate PR)

This is the test of whether the mechanism holds up — core stays clean and SA Path gets three
rows. `SapathConfig.ready()` (`sapath/apps.py`) already imports its signal modules; it gains
`from sapath.signals import integration_status  # noqa: F401`.

`sapath/signals/integration_status.py`:

- **Helix NGS Database** — from `HelixNGSOrderSync.get_last_run()` / `get_last_update()`, plus
  `HelixNGSOrder.objects.count()` as `record_count` and a link to the orders grid. Trigger
  `sapath_helix_load_if_changed.delay()`. `warning_age=timedelta(days=1)` — the CSV drop is daily,
  so a day of silence is worth a Slack line.
- **Mocha API** — outbound. `sapath/tasks/import_mocha.py::import_mocha_sample_extractions` wraps
  its body in `IntegrationActivity.track("sapath-mocha-api", name="Mocha API",
  direction=IntegrationDirection.OUTBOUND)` and sets `changed` from the returned list, which is
  enough for §3 to render the row unaided. A thin receiver adds `record_count` and the grid link
  on top.
- **VariantGrid API (lab client)** — §9's middleware, switched on in the SA Path env settings for
  the specimen-measure endpoints the pipeline posts to.

---

### 9. Incoming VariantGrid API calls

`eventlog/middleware.py` gains `IntegrationApiMiddleware`, appended to `MIDDLEWARE` in
`default_settings.py`:

```python
INTEGRATION_API_TRACKING = {}   # {url path prefix: display name}, eg {"/patients/api/": "Lab pipeline API"}
```

On response, when `request.path_info` starts with a configured prefix, it calls
`IntegrationActivity.record()` against `api:<resolved view_name>` — success for a 2xx/3xx,
error carrying the status code for 4xx/5xx. One indexed row-level `UPDATE` per API request,
against endpoints that are already doing real work. Empty dict (the default) means the middleware
returns immediately, so deployments that do not want it pay nothing.

§3 then renders those rows with no further wiring, which answers "when did the lab client last
send us anything?" — the VariantGrid equivalent of Mocha's *API Reports* card.

---

## Delivery order

1. **Framework** — §1 `library/integration_status.py`, §2 `IntegrationActivity` + migration,
   §3 auto-provider, §4 the Server Status section.
2. **Core providers** — §7 sync, classification imports, ClinVar export. First real rows.
3. **Triggers** — §5 the `Sync now` button path.
4. **Slack bridge** — §6.
5. **Incoming API** — §9 middleware.
6. **SA Path** — §8, in `../variantgrid_sapath`, once 1–5 are merged.

Each step is independently useful and independently reviewable; 1 and 2 are the bulk of it.

## Testing

Worth keeping (all cover logic written here, none restate Django):

- `library/tests/test_integration_status.py` — `get_integration_statuses()` collects from a
  temporarily-connected receiver, drops `None` returns, flattens lists, and reports a raising
  receiver as a message while still returning the other receivers' statuses.
- `run_integration_trigger` runs a registered `action_id` and leaves an unregistered one alone.
- `eventlog/tests/test_integration_activity.py` — `track()` on the clean path stamps
  `last_success`, with `changed` set also stamps `last_change`; on exception stamps `last_error` +
  message and re-raises; counters accumulate over repeated calls.
- §3's receiver: `Last Failure` is `danger` when newer than `last_success` and `warning` when
  older.
- §9 middleware: records against a configured prefix, ignores an unconfigured one, records error
  status for a 4xx.
- The Server Status section renders an integration with no run history, and hides itself entirely
  when nothing is registered.

## Settings summary

| Setting | Default | Effect |
|---|---|---|
| `INTEGRATION_STATUS_ENABLED` | `True` | Renders the Integrations section on Server Status |
| `INTEGRATION_API_TRACKING` | `{}` | Path prefix → card name for incoming API recording (§9) |

Everything else is registration in code, so a deployment's integrations arrive with its apps.

## Decisions taken

- **`eventlog` is the home for `IntegrationActivity`** — it is the audit/activity app, is always
  installed, and already has migrations (latest `0003_viewevent_referer`). `variantopedia` owns the
  page but has no models or migrations of its own.
- **One list, ordered by `sort_order` then name.** `direction` stays on the dataclass and in the
  model so a later page can group by it; the Server Status section renders a single flat run of
  rows, which is what keeps it looking like the Celery Workers section above it.
- **`warning_age` left unset keeps an integration off the nightly Slack digest** (§6). Each
  provider opts in with a period that suits its cadence — daily for the Helix CSV drop, unset for
  classification imports, which arrive whenever a lab sends them.
- **Receivers are imported inside each app's `ready()`**, with the `# pylint:
  disable=import-outside-toplevel,unused-import` / `# noqa: F401` comments the existing configs
  carry (`sync/apps.py`, `eventlog/apps.py`) — that is the established registration point, and the
  `noqa` is what stops an autofix from silently unregistering a receiver.
