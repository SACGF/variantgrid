# eventlog — agent notes
Owns: Event, ViewEvent, IntegrationActivity and ApiRequestLog (one row per call to an `INTEGRATION_API_TRACKING`
endpoint), the middleware that writes them, and the API Log page.
Start with: `eventlog/models.py` · `eventlog/middleware.py:IntegrationApiMiddleware` · `eventlog/grids.py:ApiRequestLogColumns`

- `ApiRequestLog.changes` / `num_changes` are what the call wrote, collected while the view runs by
  `library/django_utils/change_tracking.py:collect_changes`: creates and updates from `post_save`, deletes for the models
  passed to `track_deletes` (seqauto's, from `seqauto/apps.py`), unchanged rows from `update_or_create_if_changed` /
  `get_or_create_recorded`. A tracked endpoint that writes with `bulk_create`, `QuerySet.update` or a fast delete
  records those rows itself with `active_collector().record(...)`, since no signal fires. Rows from before #1413 hold
  `{}`, meaning "not recorded" rather than "no changes".
- Deletes are opt-in per model because a `post_delete` receiver stops Django fast-deleting that model; never connect
  one without a sender.
