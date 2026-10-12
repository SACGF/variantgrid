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

    class Meta:
        indexes = [models.Index(fields=["view_name", "created"])]
```

The model has no link to domain objects. Searching `request_body` for a sequencing run name answers "what was sent
for run X" while the table is seqauto-sized. A generic object link can come later if people use the log.

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

`eventlog/tests/test_urls.py` covers the grid, with a superuser.
