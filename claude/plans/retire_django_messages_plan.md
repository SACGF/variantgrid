# Retire the `django_messages` fork — in-repo markdown messaging

## Background

VariantGrid depends on `django-messages` (upstream abandoned; we maintain a fork to track Django
releases). Its functional footprint here is tiny, and we already own every template. This plan
replaces it with a small in-repo app so the fork can be deleted.

### How it is used today
- **Model / manager**: `Message` (subject, body, sender, recipient, parent_msg, sent_at, read_at,
  replied_at, sender_deleted_at, recipient_deleted_at) + `inbox_for` / `outbox_for` / `trash_for`
  + `inbox_count_for`.
- **Views** (all login-gated by global middleware): `inbox`, `outbox`, `trash`, `compose`,
  `reply`, `view`, `delete`, `undelete`.
- **URLs**: included at `messages/` behind `settings.INBOX_ENABLED` (`True` by default, `False` on
  Shariant).
- **Programmatic sends — 2 sites, both system→user notices**:
  - `upload/vcf/vcf_import.py:433` — "VCF imported as vcf #N" (currently raw `<a href>` HTML).
  - `snpdb/signals/signal_handlers.py:41` — initial org/lab welcome (from `USER_CREATE_ORG_MESSAGE`,
    currently plain text).
- **UI**: nav inbox icon + unread badge via `{% inbox_count %}` (`uicore/.../base.html:190`);
  "Write internal message" link on the user profile (`snpdb/.../view_user.html:18`).
- **Email**: package `post_save` signal `new_message_email` emails the recipient on every message
  via raw `django.core.mail.send_mail` (currently active).
- **Templates**: all 6 already live in `variantgrid/templates/default_templates/django_messages/`.

## Goal

A self-contained `user_messages` app that preserves current behaviour (inbox, unread badge,
compose/reply, delete/undelete, system notices), stores bodies as **markdown**, renders them
**safely** (markdown → sanitized HTML), routes new-message email through our logged
`EmailLog.send_mail` infrastructure, preserves existing message rows, and lets us remove the
`django-messages` dependency and fork.

## New app: `user_messages`

### Model (`user_messages/models.py`)
Port the `Message` model verbatim in fields (so existing rows map 1:1), plus:
- Manager methods `inbox_for` / `outbox_for` / `trash_for` and module-level `inbox_count_for`.
- `sender` FK `on_delete=PROTECT`, `recipient` FK `on_delete=SET_NULL null=True` (as today).
- `get_absolute_url()` → `reverse('messages_detail', ...)`.
- `save()` stamps `sent_at` on first save.
- A `body_html` property that renders markdown → sanitized HTML (see below).
- `body` is documented as **markdown**.

### Safe markdown rendering
Add a helper (e.g. `user_messages/rendering.py` or a method on the model) that pipelines:
`markdown.markdown(body)` → `library.utils.sanitize_html(...)` → `SafeString`.
Reuse the existing `sanitize_html()` (bleach whitelist) added for the XSS fix so raw HTML embedded
in markdown is still neutralized. Templates render `{{ message.body_html }}` (already safe), so the
`{% autoescape off %}` block is retired.

### Views (`user_messages/views.py`)
Port `inbox`, `outbox`, `trash`, `compose`, `reply`, `view`, `delete`, `undelete` with the same
behaviour and the same sender/recipient authorization checks (`get_object_or_404` + 404 when the
user is neither party; mark `read_at` on first view by recipient). Login and CSRF come from global
middleware, so views carry no per-view decorators. The compose success flow keeps using the Django
contrib `messages.info(...)` confirmation.

### Form (`user_messages/forms.py`)
Port `ComposeForm` (recipient / subject / body) and the `CommaSeparatedUserField` recipient widget.
The body field label notes markdown is supported; keep the `Textarea`.

### URLs (`user_messages/urls.py`)
Reproduce the **same URL names** so nav, the compose link, and templates keep working unchanged:
`messages_redirect`, `messages_inbox`, `messages_outbox`, `messages_compose`,
`messages_compose_to`, `messages_reply`, `messages_detail`, `messages_delete`,
`messages_undelete`, `messages_trash`. Wire into `variantgrid/urls.py` behind `INBOX_ENABLED`,
replacing the `django_messages.urls` include.

### Templatetag (`user_messages/templatetags/user_messages_tags.py`)
Provide `{% inbox_count %}` (and `{% inbox_count as var %}`) returning the unread count for the
current user. Update `uicore/.../base.html` to `{% load %}` the new tag library.

### Templates
Move the 6 existing templates from `default_templates/django_messages/` to the `user_messages`
namespace, adjust `{% extends %}` / `{% load %}` and the message-body block to render
`{{ message.body_html }}`. `compose.html` gains a short "markdown supported" hint.

## Data migration (adopt the existing table)

Preserve existing rows with zero copying by taking over the `django_messages_message` table:
- `user_messages/migrations/0001_initial.py` sets the model's `Meta.db_table = 'django_messages_message'`
  and wraps `CreateModel` in `migrations.SeparateDatabaseAndState(state_operations=[...])` so Django
  records the model state without issuing `CREATE TABLE` (the table already exists).
- `0002_rename_table.py` runs `AlterModelTable` to rename the physical table to
  `user_messages_message` (a cheap `ALTER TABLE ... RENAME`), and drops the interim `db_table`
  override from the model `Meta`.

## Send-site updates (markdown)

- `upload/vcf/vcf_import.py` — build the body as markdown: `f"VCF {vcf.name} imported as [vcf #{vcf.pk}]({url})"`,
  importing `Message` from `user_messages.models`.
- `snpdb/signals/signal_handlers.py` — import `Message` from `user_messages.models`; the
  `USER_CREATE_ORG_MESSAGE` text is already plain and renders fine as markdown.

## Email notification via our infrastructure

Replace the package's raw `send_mail` signal with an explicit send in the message-create path (form
`save()` and the two system sites, or a single helper both call) that emails the recipient through
`email_manager.models.EmailLog.send_mail(...)` — the same logged, `SEND_EMAILS`-gated path that the
`NotificationBuilder` family uses internally. It emails `recipient.email` when present, records an
`EmailLog`, and renders the markdown body for the email content.

## Settings & dependency removal

- Remove `'django_messages'` from `INSTALLED_APPS`; add `'user_messages'`.
- Point `variantgrid/urls.py` at `user_messages.urls`.
- Fix the two spurious imports to use Django directly:
  `from django.contrib.auth.models import User` in `annotation/vcf_files/import_clinvar_vcf.py` and
  `classification/signals/classification_hooks_significant_change.py`.
- Remove `django-messages` from `requirements.txt` / `requirements-django5.txt` and delete the fork
  reference once deployed.

## Tests (`user_messages/tests/`)
- Model/manager: inbox/outbox/trash partitioning, `inbox_count_for`, `sent_at` stamping.
- Rendering: markdown links render; embedded `<script>` / `onerror` / `javascript:` hrefs are
  neutralized by `sanitize_html`; newlines preserved.
- Views: compose creates a row + emails recipient; `view` marks `read_at`; non-party access 404s;
  delete/undelete toggle the correct `*_deleted_at`.
- URL smoke test via `URLTestCase._test_urls()` for the inbox/compose/view/delete routes.
- Migration check: `makemigrations --check` clean; the `0001` state-only + `0002` rename apply on a
  DB seeded with legacy `django_messages_message` rows and preserve them.

## Rollout ordering
1. New app + model + rendering + views + forms + urls + templatetag + templates.
2. State-only `0001` (db_table override) then `0002` table rename.
3. Repoint send sites, nav `{% load %}`, and `variantgrid/urls.py`; fix spurious `User` imports.
4. Remove `django_messages` from settings + requirements.
5. Tests green (`--keepdb`), `makemigrations --check` clean.
