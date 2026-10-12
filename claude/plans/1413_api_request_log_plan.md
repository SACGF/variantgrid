# #1413 — Log API requests

Written by Claude Opus 5.5 (claude-opus-5-5), 2026-10-12
Status: in progress

[#1413](https://github.com/SACGF/variantgrid/issues/1413). External pipelines post to the seqauto and specimen
measure APIs. When a call goes wrong, all we have afterwards is whatever records it wrote. `eventlog/middleware.py:IntegrationApiMiddleware`
already watches the prefixes in `INTEGRATION_API_TRACKING`, but it only keeps one rolling "last seen" row per endpoint
(`eventlog/models.py:IntegrationActivity`). This plan keeps one row per call for those same prefixes, so an admin
can see what was sent and what came back.

## Model

```python
class ApiRequestLog(models.Model):
    """ One row per call to an INTEGRATION_API_TRACKING endpoint, kept forever unless API_REQUEST_LOG_DAYS is set """
    created = models.DateTimeField(auto_now_add=True, db_index=True)
    user = models.ForeignKey(User, null=True, on_delete=SET_NULL)
    method = models.CharField(max_length=10)
    path = models.TextField()
    query_string = models.TextField(blank=True)
    view_name = models.TextField()
    view_kwargs = models.JSONField(default=dict, encoder=DjangoJSONEncoder)
    status_code = models.IntegerField()
    duration_ms = models.IntegerField()
    request_body = models.TextField(null=True)   # trimmed, secrets redacted
    response_body = models.TextField(null=True)  # trimmed, only kept for status >= 400
    content_type = models.ForeignKey(ContentType, null=True, on_delete=SET_NULL)  # the view's api_log_model
    object_pks = ArrayField(models.TextField(), default=list)  # empty when the call failed before finding one

    class Meta:
        indexes = [models.Index(fields=["view_name", "created"])]
```

### Linking a call to the objects it was about

- A view names the model its calls are about with a class attribute, `api_log_model`. The middleware reads it from
  `request.resolver_match.func.cls` (DRF) or `.view_class` (Django), so a call that fails before the view finds
  anything still records `content_type`.
- As the view resolves objects it calls `eventlog/models.py:log_api_request_objects` with their pks, which collects
  them on the Django request. The middleware writes them into `object_pks` on the same INSERT.
- `ApiRequestLog.for_object(content_type, pk)` returns the calls whose `object_pks` contain pk, plus the calls for that
  model with empty `object_pks` (failed before finding one) whose request body has pk as a JSON string. The body
  search runs only over the failures.
- Seqauto: the run-scoped viewsets subclass `seqauto/views_rest.py:SequencingRunModelViewSet`, which logs the run of
  each record from `get_object`, `perform_create` and `perform_update`. The bulk-create and link-extraction views
  subclass `seqauto/views_rest.py:SequencingRunWriteAPIView` and log the runs of what they saved. Catalogue endpoints
  (enrichment kit, sequencer, experiment, variant caller) declare no model and get no link.

## Middleware

`IntegrationApiMiddleware` stays the only entry point, and `INTEGRATION_API_TRACKING` stays the opt-in, so a
deployment with an empty setting writes nothing.

- **Before the view:** it captures the request body, because DRF consumes the stream.
  - JSON or form bodies up to `API_REQUEST_LOG_MAX_READ_BYTES` (1 MB) are read through `request.body`. Django caches
    it, so the view parses it as normal.
  - Bigger bodies and `multipart/*` are recorded as a placeholder (`[multipart/form-data, 12345 bytes]`), without
    reading them. This keeps uploads out of memory and avoids `DATA_UPLOAD_MAX_MEMORY_SIZE` errors.
  - In a JSON body, any key containing `password`, `token`, `secret` or `api_key` is replaced with `"[redacted]"`, at
    any depth.
  - The text is trimmed to `API_REQUEST_LOG_MAX_BODY` (10,000 characters).
- **Around the view:** `time.perf_counter()` gives `duration_ms`.
- **After the view:** it records the `IntegrationActivity` update as before, plus one `ApiRequestLog` row.
  - The user is `request.user` once the view has run. DRF copies its authenticated user onto the Django request.
  - The response body is kept only for status ≥ 400 on non-streaming responses, trimmed the same way.

Writing the row is one INSERT in the request. That costs nothing at seqauto volumes (fewer than 100 calls a week).
Measure it before adding a busy prefix such as the Shariant classification API.

## Retention

Rows are kept forever by default: the expected volume is a few thousand calls a year, and old calls help debugging.
A deployment that sets `API_REQUEST_LOG_DAYS` gets `eventlog/tasks.py:prune_api_request_log`, scheduled daily on
`db_workers` in `variantgrid/celery.py`, which deletes rows older than that many days.

## UI

- `eventlog/grids.py:ApiRequestLogColumns` has these columns: created, user, method, view name, status, duration and
  path.
  - The search box matches path, view name and request body.
  - A `view_name` param narrows the grid to one endpoint, and an "errors only" checkbox keeps status ≥ 400.
  - Expanding a row loads `api_request_log_detail`, which shows the request and response bodies.
- Page `api_request_log` (Settings menu → "API Log", admin only). The page, the grid and the detail view are all
  `require_superuser`.
- The sequencing run page has an "API Log" tab for superusers, shown once any call has been linked to a
  SequencingRun. It loads `api_request_log_object_tab` (app label, model, pk), whose grid is
  `ApiRequestLog.for_object`.
- `eventlog/signals/integration_activity_status.py:integration_status_for_activity` sets `url` on rows whose key
  starts with `api:` to the API log filtered by that view name. The Server Status integration box then links
  straight to that endpoint's calls.

## Tests

`eventlog/tests/test_integration_api_middleware.py` covers:
- a row is written with status, user and body
- JSON secrets are redacted at depth
- a multipart body is not read
- a response body is kept only for an error
- nothing is written for an unconfigured prefix

`seqauto/tests/test_api_request_log.py` covers a call linked to its run, a failed call found by the run name in its
body (and not by a run whose name it only starts), and the run page tab.

`eventlog/tests/test_urls.py` covers the grid, with a superuser.
