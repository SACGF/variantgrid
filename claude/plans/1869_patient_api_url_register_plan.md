# #1869 — URLS_NAME_REGISTER does not gate DRF router URLs

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-14
Status: draft

[#1869](https://github.com/SACGF/variantgrid/issues/1869), split out of the TSO 500 plan
(SACGF/variantgrid_sapath#431). No models change; the data section is empty.

## The problem, verified

`variantgrid/perm_path.py:_perm_path` wraps a view in `library/django_utils/__init__.py:require_superuser`
when `settings.URLS_NAME_REGISTER[name]` is False. It only sees the names passed to `perm_path.path`, so a
DRF router's patterns, built inside `rest_framework.routers.SimpleRouter.get_urls` with Django's own
`re_path` and handed over as `path('', include(router.urls), name='patients_apis')` in `patients/urls.py`,
are never consulted. The nine `api_patient-list` / `-detail` style names that
`variantgrid/settings/env/shariantcommon.py` sets False are dead entries; only
`api_specimen_measure_bulk_create`, registered through `perm_path.path`, is enforced. Confirmed by reading
both call paths; no Shariant box was needed.

Drift found on the way:

- `name='patients_apis'` (and `seqauto_apis` in `seqauto/urls.py`) is silently dropped by Django - `_path`
  builds a `URLResolver` for an `include()` and ignores `name` - so it is not a resolvable name and
  `get_visible_url_names` never lists it. `perm_path` does look it up, though, and a deployment that ever set
  it False would get `require_superuser(<include tuple>)`: `functools.wraps` tolerates the tuple, the route
  becomes a `URLPattern` at `''`, and every API call would die with `TypeError` when the wrapper calls the
  tuple. Never exercised; removed below.
- `seqauto/urls.py` has the same shape but is app-gated on Shariant (`URLS_APP_REGISTER["seqauto"] = False`),
  so nothing leaks there today; it gets the same fix for consistency.
- `DefaultRouter` also mounts `api-root` at `/patients/` and `/seqauto/`, listing the registered endpoints.
  The register defaults it True and this plan leaves it alone.

## Decision: make router registration honour the register

The register is the one mechanism for "this deployment does not serve that URL": `perm_path` enforces it,
`variantgrid/perm_path.py:get_visible_url_names` drives menus and tips from it, and the Shariant names are
already written. Applying it to router patterns makes those entries true and fixes the class of drift for
any future router. The alternative - skipping `router.register` on those deployments - needs a second
setting or a register lookup at registration time, leaves the existing names dead, and gives 404 where every
other unregistered URL gives 403 to a non-superuser.

Semantics follow `perm_path`: a name set False is superuser-only, not absent. `require_superuser` raises
`PermissionDenied` before DRF's dispatch, so a non-superuser API client gets Django's 403 (anonymous too,
rather than DRF's 401) - fine for an endpoint the deployment does not offer.

## Code changes

`variantgrid/perm_path.py` - one new function next to `path`:

```python
def router_urls(router) -> list[URLPattern]:
    """ A DRF router's patterns with URLS_NAME_REGISTER applied, as path() does for a named view """
```

For each pattern in `router.urls` whose `name` the register sets False, return
`URLPattern(url.pattern, require_superuser(url.callback), url.default_args, url.name)`; others pass through.
`require_superuser` uses `functools.wraps`, so the attributes DRF's `as_view()` leaves on the callback
(`cls`, `actions`, `csrf_exempt`) survive. Format-suffix variants (`.json`) share their name and are covered.
Extend the module docstring: named views go through `path`, routers through `router_urls`.

`patients/urls.py` and `seqauto/urls.py` - replace `path('', include(router.urls), name='..._apis')` with
`urlpatterns += router_urls(router)`, dropping the `include` import where nothing else uses it. The
`api_specimen_measure_bulk_create` path stays ahead of the router patterns.

`variantgrid/settings/env/shariantcommon.py` - unchanged; its nine API names now take effect.

## Tests

New module variantgrid/tests/test_perm_path.py, a plain `TestCase` (no DB rows needed; `RequestFactory` and users
from `User.objects.create_user` / `create_superuser`):

- Register `patients/views_rest.py:PatientViewSet` on a fresh `DefaultRouter` as `api_patient`, call
  `router_urls` under `override_settings(URLS_NAME_REGISTER=defaultdict(lambda: True, {"api_patient-list": False}))`
  (the `variantgrid/tests/test_tips.py` pattern), pick the pattern named `api_patient-list`: a non-superuser
  GET raises `PermissionDenied`, a superuser GET returns 200.
- Same build: `api_patient-detail` (True) returns the untouched viewset callback.

The register is read when a URLconf is imported, so the test builds its own router rather than reversing
through the global one (`snpdb/tests/test_lab_members.py` records that trap).
`patients/tests/test_urls.py:Test.testApiUrls` keeps covering the default-True path through the real URLconf.

## Manual verification

On this box `vg page /patients/api/v1/patient/` still returns 200 as `claude_agent` (register True). On a
Shariant test deployment after deploy: `curl -u <non-superuser>` against `/patients/api/v1/patient/` is 403,
the same as `/patients/api/v1/specimen_measure/bulk_create` is today, and the Patients menu is still absent.

## Definition of done

- `scripts/vg tests --explain` names `variantgrid.tests.test_perm_path`, `patients.tests.test_urls` and
  `seqauto.tests.test_urls`, and they pass.
- One line in `library/CLAUDE.md` (or the `perm_path` docstring) next to the `perm_path` note: router
  patterns go through `router_urls` or the register does not apply to them.
- `scripts/vg map` refreshed (URL map) and `scripts/vg docs check` passes.
