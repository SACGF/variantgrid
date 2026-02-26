# VariantGrid Email Manager App — Reference Document

## Purpose and Overview

The `email_manager` app provides a centralized email sending and logging system. It records all email operations for audit purposes, supports per-recipient emails to prevent address leakage, and provides test modes for non-production environments.

**Key Design:** All emails logged regardless of send success; privacy mode sends individual emails per recipient; SEND_EMAILS setting for safe testing.

---

## Models

### EmailLog
Records all email send operations.
- Fields: subject, html, text, from_email (nullable), recipient_list (semicolon-separated), probably_sent (bool), single_email (bool — True if batch mode)
- TimeStampedModel (created, modified)
- Implements PreviewModelMixin with fa-solid fa-envelope icon and "Emails" category
- `get_absolute_url()` → `/email_manager/detail/<pk>/`

#### `EmailLog.send_mail(subject, html, text, from_email, recipient_list, allow_users_to_see_others=False)` → bool
Main email sending interface.

**Batch mode** (`allow_users_to_see_others=True`): Single email to all recipients (all addresses visible to each other). Uses Django's `send_mail()`.

**Privacy mode** (`allow_users_to_see_others=False`, default): Separate `EmailMultiAlternatives` per recipient, sharing a single SMTP connection for efficiency. Prevents recipient list leakage.

- Always creates EmailLog record regardless of send success
- HTML attached as text/html alternative to plain text
- Checks `settings.SEND_EMAILS` flag — if False, logs but does not send

---

## Views and URL Patterns

```
''                        → email_manager_view()  — Main email log page (superuser only)
'detail/<int:email_id>'   → email_detail()         — Email detail view (superuser only)
'pure/<int:email_id>'     → email_pure()           — Raw HTML email view
'datatable'               → DatabaseTableView(EmailColumns) — AJAX datatable
```

### Key Views
- `email_manager_view(request)` — `@require_superuser`; renders email log with datatable
- `email_detail(request, email_id)` — `@require_superuser`; parses recipients, matches to User objects, shows unrecognized addresses
- `email_pure(request, email_id)` — Renders raw HTML email body (for template debugging)

---

## Admin Interface (EmailLogAdmin)
- `has_add_permission` → False (emails created programmatically)
- `has_change_permission` → False (immutable)
- `has_delete_permission` → False (preserve audit trail)
- List per page: 500; default ordering: -created
- Search fields: recipient_list, subject, text

---

## Datatable Configuration (EmailColumns)
Columns: id, created (desc default), subject, recipient_list, probably_sent (bool renderer, error highlight if false).
- Row expansion via AJAX (email_detail view)
- Superuser-only access
- power_search on recipient_list and subject

---

## Integration Points

| App | Integration |
|-----|-------------|
| Classification | Sends classification notifications and updates |
| snpdb | Sends system notifications (lab notifications) |
| Health check (eventlog) | email_health_check signal receiver counts recent emails |

---

## Settings

- `SEND_EMAILS` (bool) — Controls whether emails are actually sent. False in dev/test (emails logged but not sent). True in production.

---

## Privacy Considerations

- **Per-recipient mode (default)**: Each recipient gets individual email; doesn't see other recipients. Use for sensitive notifications.
- **Batch mode**: All recipients listed in "To:" field; visible to each other. Use for public announcements.
