# eventlog — agent notes
Owns: Event, ViewEvent, IntegrationActivity and ApiRequestLog (one row per call to an `INTEGRATION_API_TRACKING`
endpoint), the middleware that writes them, the API Log page, and notifying humans (Slack, admin email).
Start with: `eventlog/models.py` · `eventlog/middleware.py:IntegrationApiMiddleware` · `eventlog/grids.py:ApiRequestLogColumns` ·
`eventlog/notifications.py:NotificationBuilder`

- Notify humans with `eventlog/notifications.py:AdminNotificationBuilder` (Slack, plus superuser email when
  `is_communication=True`) or snpdb's `LabNotificationBuilder`; add blocks with `add_header`/`add_field`/`add_markdown`
  (Markdown, converted to HTML and Slack) and call `send()`. `eventlog/notifications.py:send_notification` is the raw
  Slack hook, falls back to `report_event` (Rollbar info + an Event) when Slack is unconfigured, and truncates blocks at
  `SLACK_CHARACTER_LIMIT`.
- `ApiRequestLog.changes` / `num_changes` are what the call wrote, collected while the view runs by
  `library/django_utils/change_tracking.py:collect_changes`: creates and updates from `post_save`, deletes for the models
  passed to `track_deletes` (seqauto's, from `seqauto/apps.py`), unchanged rows from `update_or_create_if_changed` /
  `get_or_create_recorded`. A tracked endpoint that writes with `bulk_create`, `QuerySet.update` or a fast delete
  records those rows itself with `active_collector().record(...)`, since no signal fires. Rows from before #1413 hold
  `{}`, meaning "not recorded" rather than "no changes".
- Deletes are opt-in per model because a `post_delete` receiver stops Django fast-deleting that model; never connect
  one without a sender.
- `ApiRequestLog.notes` are free-text lines from `eventlog/models.py:api_request_note`, kept on the call's
  `ChangeCollector` so code with no request in hand (a signal handler) can add one; outside a tracked call it only logs.
