# VariantGrid Eventlog App — Reference Document

## Purpose and Overview

The `eventlog` app provides comprehensive activity logging, recording both automatic view access events and explicit business logic events. Creates an audit trail of user actions for compliance, debugging, and analytics.

**Key Design:** Middleware-based automatic view logging + explicit API for business events; superuser sees all, non-staff sees own only.

---

## Models

### ViewEvent
Automatically recorded HTTP request events.
- Fields: user FK (nullable), view_name (text), args (JSON), path (text), method (text), referer (text, nullable)
- Properties: `is_get` — True if method is GET or empty
- Ignored paths: api, datatable, citations_json; text: detail, metrics; AJAX requests; URL suffixes: _detail, _autocomplete

### Event
Explicit business logic events created via API.
- Fields: user FK (nullable), date (DateTimeField), app_name (text), name (text), details (text, nullable), severity (LogLevel: I/W/E/D), filename (text, nullable)
- Methods: `can_write(user_or_group)` — True if superuser or owns event

---

## Views and URL Patterns

```
''                 → eventlog()                — Main event log page
'detail/<int:pk>'  → eventlog_detail()         — Event detail view
'create_event'     → create_event()            — POST endpoint to create events
'datatable'        → DatabaseTableView(EventColumns) — AJAX datatable
```

### Key Views
- `eventlog(request)` — Renders main event log template with datatable
- `eventlog_detail(request, pk)` — Full details, filename, metadata for single event
- `create_event(request)` — POST only; accepts app_name, event_name, details, severity; creates Event record
- `@eventlog_view` decorator — Optional decorator to record view access; adds ~100ms per request

---

## Middleware and Logging

### PageViewsMiddleware
Automatically records HTTP requests as ViewEvent entries.

**process_view() — Before view:**
1. Skip AJAX requests
2. Check if app in LOG_ACTIVITY_APPS setting
3. Resolve URL to view name
4. Extract and normalize parameters (bool/int conversion, strip csrfmiddlewaretoken, URL-decode ontology terms)
5. Create ViewEvent and attach to request

**__call__() — After view:**
1. Save ViewEvent if attached to request
2. Skip redirect responses (except search)

### EventLogHandler
Custom Python logging handler that persists log records to Event model.
- Maps log level → LogLevel enum
- Extracts request object and user
- Filters out missing .js.map file errors
- Creates Event record with exception traceback if present

---

## Datatable Configuration (EventColumns)

Columns: date (desc default), severity (custom renderer), user, app_name, name, data (truncated 75 chars), id (hidden).

Filtering:
- Non-staff: own events only
- Query param `filter`: logins, errors, warnings_and_errors, events, searches
- `exclude_admin` JSON param removes superusers and bots

---

## Health Check Integration

**active_users_health_check()** — Queries Event and ViewEvent for activity in time window; separates superusers from regular users; returns HealthCheckRecentActivity with counts and usernames.

**email_health_check()** — Queries EmailLog for recent emails; counts by subject; shows summary.

---

## Integration Points

| App | Integration |
|-----|-------------|
| All apps | Any view can call create_event() to log explicit events |
| Email manager | email_health_check signal receiver queries EmailLog |
| Library | health_check signal integration |

---

## Settings

- `LOG_ACTIVITY_APPS` — List of app names to automatically log view access for
