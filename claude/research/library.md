# VariantGrid Library Module — Reference Document

## Purpose and Overview

`library/` is a comprehensive shared utility collection for VariantGrid. It is not a traditional Django app but a reusable library module that other apps depend on. It provides foundational services: caching, permissions, logging/notifications, health monitoring, graphs, and diverse utility functions.

---

## Caching and Performance

### timed_cache.py
- `timed_cache(size_limit, ttl, quick_key_access)` — Decorator for function result caching with size limits and TTL
- `clear_cached_property()` — Utility to clear cached_property instances from objects

---

## Permissions and Authorization

### guardian_utils.py
Django-guardian based permission system utilities.

**DjangoPermission** (static permission formatter):
- Constants: `READ = 'view'`, `WRITE = 'change'`
- `perm(obj, permission)` — Generate permission strings

**Utility functions:**
- `all_users_group()`, `public_group()`, `bot_group()` — Standard group getters
- `admin_bot()` — Get or create admin bot user (unit-test aware)
- `assign_permission_to_user_and_groups()` — Assign permissions to user and their groups
- `clear_permissions()` — Remove permissions from all users/groups
- `check_can_write()`, `check_can_delete()` — Permission checking with exceptions
- `groups_map()` — Fetch groups for multiple objects efficiently
- `highest_group()` — Determine highest permission level group

### django_utils/guardian_permissions_mixin.py

**GuardianPermissionsMixin**:
- `can_view(user)`, `can_write(user)` — Check permissions
- `filter_for_user(qs)` — Filter querysets by permissions
- `get_for_user(pk)` — Get single object with permission check

**GuardianPermissionsAutoInitialSaveMixin**:
- Automatically assigns permissions on first save

---

## Logging, Notifications, and Health Checks

### log_utils.py

**NotificationBuilder** — Main notification system supporting multiple output formats:
- Block classes: `HeaderBlock`, `FieldsBlock`, `MarkdownBlock`, `DividerBlock`
- Methods: `add_header()`, `add_field()`, `add_markdown()`, `add_divider()`, `merge()`
- Output: `as_text()`, `as_html()`, `as_slack()`

**AdminNotificationBuilder** — Extends NotificationBuilder for admin notifications via email

**Key functions:**
- `send_notification()` — Post to Slack webhooks with character limit handling
- `report_event()` — Log events to Rollbar and EventLog
- `report_message()` — Report non-fatal messages
- `report_exc_info()` — Report exceptions with tracebacks
- `log_admin_change()` — Log Django admin changes
- `log_saved_form()` — Log form saves with field changes
- `get_current_logged_in_user()` — Get authenticated user from thread locals

### health_check.py
System health monitoring with Slack notifications.
- `HealthCheckRequest` — Context wrapper with time window
- `HealthCheckStat` (abstract) — Base class for health statistics
- `HealthCheckRecentActivity` — Activity in time windows (created/modified/deleted); `simple_report()`
- `HealthCheckTotalAmount` — Total count statistics
- `HealthCheckCapacity` — Disk space and resource metrics
- `HealthCheckAge` — Age-based checks for periodic tasks
- `populate_health_check()` — Aggregate all health checks into notification
- `health_check_signal` — Django signal for apps to report health

### uptime_check.py
Public-facing uptime/status monitoring.
- `UptimeCheckStatus` — Enum: OKAY, NON_CRITICAL_FAILURE, CRITICAL_FAILURE, EXCEPTION
- `UptimeCheckResponse` — Single check result (name, status, note)
- `UptimeCheckOverall` — Aggregated results
- `retrieve_uptime_response()` — Collect all uptime checks
- `uptime_check_signal` — Signal for apps to register checks

---

## Data Handling and Serialization

### preview_request.py
System for generating rich previews/summaries of data objects.

**PreviewData** (dataclass):
- Fields: category, identifier, title, icon, summary, URLs
- `for_object()` — Factory from Django models
- `as_json()` — Serialize to JSON

**PreviewModelMixin** — Mixin for models to support previews:
- `preview_category()` — Human-readable category
- `preview_icon()` — Font awesome icon class
- `preview_enabled()` — Check if previews active

**PreviewKeyValue** (dataclass):
- `count()` — Create count previews

**Signals:**
- `preview_request_signal` — Apps register their preview handlers
- `preview_extra_signal` — Additional preview data

---

## Graphs and Visualization

### graphs/graph_base.py

**GraphBase** (abstract):
- `plot()` (abstract) — Draw the graph
- `save()` — Save to PNG via matplotlib
- `decorations()` — Set title/labels
- `post_plot()` — Add legends

---

## jQGrid Support

### jqgrid/jqgrid.py

**JqGrid** (base class, legacy):
- Fields: queryset, model, fields
- `json_encode()` — DjangoJSONEncoder wrapper
- Various operation formatting utilities

---

## Authentication and OAuth

### oauth.py

**ServerAuth** — Multi-auth support (Basic, OAuth2):
- `for_sync_details()` — Factory from settings dict
- `keycloak_connector()` — Get Keycloak connector
- `auth` property — Returns auth handler (OAuth2Session or HTTPBasicAuth)
- HTTP methods: `get()`, `post()`

### keycloak.py

**KeycloakNewUser** — New user data container

**Keycloak** — Keycloak admin API client:
- `ping()` — Health check
- `change_password()` — Send password reset email
- `check_groups()` — Ensure groups exist, create if needed
- `existing_user()` — Look up user by email
- `add_user()` — Create new user with groups

---

## Utility Functions (utils/)

20+ utility modules for common operations:

| Module | Key Contents |
|--------|-------------|
| `collection_utils.py` | invert_dict, sorted_nicely, batches |
| `django_utils.py` | is_ajax, render_ajax_view, get_cached_project_git_hash |
| `model_utils.py` | ModelUtilsMixin, ArrayLength, AsciiValue, model_has_field() |
| `diff_utils.py` | DiffBuilder, MultiDiff |
| `json_utils.py` | JSON safety and path manipulation |
| `hash_utils.py` | MD5/SHA256 functions |
| `date_utils.py` | Date/time helpers |
| `text_utils.py` | String manipulation |
| `html_utils.py` | html_link, html_id_safe |
| `file_utils.py`, `export_utils.py`, `misc_utils.py` | Various helpers |

---

## Integration Points

| Service/App | Integration |
|-------------|-------------|
| Guardian | Django-guardian backend for object-level permissions |
| Rollbar | Exception and message reporting |
| Slack | Webhook notifications for events and health checks |
| Keycloak | OIDC provider admin API |
| Django Signals | Health check and preview discovery mechanisms |
| All apps | Guardian utilities, logging, caching, preview system |
