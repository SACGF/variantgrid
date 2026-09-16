# eventlog — research notes

Verified against 7c4408c62 on 2026-09-06

Three kinds of record live here (`claude/maps/models.md#eventlog`). `eventlog/models.py:Event` is an explicit application
event - a login, a search, an import, an error - written through `create_event` or the logging handler. `ViewEvent` is a page
view captured by middleware on deployments that ask for it, and feeds the classification view-metrics pages.
`IntegrationActivity` is a one-row-per-integration "last seen" for external systems, updated in place, which the server
status page and nightly digest read. `vg status` lists the last ten `Event` rows at ERROR (`library/vg/status.py`). The
readme is `eventlog/__eventlog_readme.md`; the pages are `claude/maps/urls.md#eventlog`.

## Flows

`eventlog/models.py:create_event` writes an `Event` with `date=now`, inferring `app_name` from the caller's module via
`inspect.stack()` when none is given, and by default also emits the row through Python `logging` at the matching level.
About 25 call sites use it (VCF archive, partition archive, zygosity count repair, classification import, graph generation).
`eventlog/models.py:create_login_event` is wired to `user_logged_in` at import time, which is where the "logins" filter on
the log page gets its rows. The browser can also post one: `eventlog/views.py:create_event` is `@require_POST`, validates the
severity against `LogLevel.CHOICES`, and trims `app_name`, `name` and `details` to fixed lengths with a `[trimmed …]` footer
(#1519); `analysis.js` uses it for analysis-editor events. `eventlog/views.py:eventlog_view` is a decorator that records a
view hit as an Event - its own docstring says not to use it, because it costs a write per request.

The second writer is the logging handler. `variantgrid/settings/components/default_settings.py` attaches
`eventlog/loggers.py:EventLogHandler` as the `db` handler of the `django` logger, so every `django.request` warning or error
raised while serving a request becomes an `Event` named `django_exception` with the message and traceback in `details`. The
handler imports the model by string at emit time because logging is configured before apps load, returns without writing
when the record carries no authenticated `request.user`, and drops the "Not Found: …js.map" noise.

Page views: `eventlog/middleware.py:PageViewsMiddleware.process_view` builds a `ViewEvent` for non-AJAX requests whose first
path segment is in `settings.LOG_ACTIVITY_APPS`, skipping `IGNORE_SEGMENTS` (api, datatable, citations_json), any path
containing `detail` or `metrics`, and view names ending `_detail` / `_autocomplete`. The view kwargs plus GET and POST are
flattened into `args`, booleans and ints coerced, `csrfmiddlewaretoken` dropped, ontology terms un-mangled back to `MONDO:…`,
and a `classification_id` of the form `id.timestamp` split into `classification_id` and `modification_timestamp`. The row is
only saved in `PageViewsMiddleware.__call__` after the view returns, and redirects are discarded except for
`variantopedia:search`, so a search that jumped straight to a result still counts. Only `shariantcommon.py` installs the
middleware and defines `LOG_ACTIVITY_APPS`. Readers are `classification/views/classification_view_metrics.py:ViewEventCounts`
and `classification/views/search_view_metrics.py`, plus `eventlog/admin.py:ViewEventAdmin` with its lab, organisation and
"exclude admin/test users" filters.

Integrations: `eventlog/models.py:IntegrationActivity.record` and the `IntegrationActivity.track` context manager stamp
`last_attempt` / `last_success` / `last_change` / `last_error` on one row keyed by a string, with counters through `F()` so
concurrent workers add up and `_upsert` doing a single UPDATE once the row exists. `eventlog/middleware.py:IntegrationApiMiddleware`
does the same for inbound API prefixes named in `settings.INTEGRATION_API_TRACKING`, recording a 4xx/5xx as an error and any
successful non-GET as a change. `eventlog/signals/integration_activity_status.py:integration_activity_status` turns every row
into an `IntegrationStatus` for the server status page with no further registration, including a "Dismiss error" trigger
that sets `last_error_acknowledged` so a failure that has been seen goes quiet until a newer one arrives.

Reading: `eventlog/grids.py:EventColumns` shows own events only unless the user is a superuser, supports the page's
`filter` (logins, errors, warnings_and_errors, events, searches) and `exclude_admin` (drops superusers, bots and user-less
rows), and expands rows through `eventlog/views.py:eventlog_detail`, which checks superuser-or-owner.
`eventlog/signals/active_users_health_check.py:active_users_health_check` answers the nightly health check by counting
distinct users with an `Event` or `ViewEvent` in the window, admins on their own line.

## Why it is shaped this way

`Event` predates Rollbar and Slack and is now mostly an in-database audit trail; the Slack and email paths in
`library/log_utils.py` cover alerting. `ViewEvent` exists for Shariant's usage reporting, which is why it is opt-in per
deployment and skips the noisy paths. `IntegrationActivity` was added because outbound and inbound integrations (Alissa,
the SA Pathology Mocha API) left nothing behind but the records they wrote, so "when did this last run" had no answer; one
mutable row per integration keeps the table bounded however chatty the integration is.

## History

Security fixes (#3819, March 2026) added the owner check to `eventlog_detail`, restricted the grid to superusers rather than
staff, and made `create_event` validate its inputs. #1519 (May 2026) added the length caps. `IntegrationActivity`,
`IntegrationApiMiddleware` and the dismiss trigger arrived with the server status integration panel in August 2026. Emoji in
event details render since #1098 (`eventlog/grids.py:EventColumns.render_data` runs `emoji_to_unicode`).

## Traps

`settings.LOG_ACTIVITY_APPS` has no default: `PageViewsMiddleware` raises `AttributeError` on a deployment that installs it
without defining the set. The middleware infers the app from the first URL segment, which is not the Django app name.
`EventLogHandler` swallows anonymous-request errors entirely - a 500 on a public path is in Rollbar, not the event log.
`Event.can_write` is superuser-or-owner and an `Event` with `user=None` is writable by admins only. Tests:
`eventlog/tests/test_integration_activity.py`, `eventlog/tests/test_integration_api_middleware.py`, `eventlog/tests/test_models.py`.
