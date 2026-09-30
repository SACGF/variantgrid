# vcauth / oidc_auth — research notes

Verified against a96540a68 on 2026-09-25

Two small apps around Django's `User`. `vcauth` is installed everywhere and only swaps the admin's User page for one
with an "Email weekly summary" action (`vcauth/user_admin.py:CustomUserAdmin`, calling the classification app's
`send_summary_email_to_user`). `oidc_auth` is installed only by the Shariant settings (`variantgrid/settings/env/shariantcommon.py`
prepends it to `INSTALLED_APPS`, adds the backend and middleware, and sets `USE_OIDC = True`): it hands login and user
management to Keycloak through mozilla-django-oidc, and Keycloak becomes the source of truth for who a user is, which
labs they are in and whether they are an admin. Login and CSRF enforcement is global middleware, covered in
`claude/guides/operations.md#authentication-surface`; URLs, the uptime signal and settings are in
[urls](../maps/urls.md), [signals](../maps/signals.md) and [settings](../maps/settings.md).

## Flows

### Browser login

`LOGIN_URL = '/oidc_login/'`, so the global login middleware sends anonymous users to `variantgrid/views.py:oidc_login`,
which calls mozilla's `oidc_authentication_init` view directly (PKCE S256, RS256 ID tokens checked against the realm's JWKS).
The `/oidc/` URLs are mounted `login_not_required`, and only when `USE_OIDC` (`variantgrid/urls.py`). On the callback,
mozilla fetches userinfo and calls `oidc_auth/backend.py:VariantGridOIDCAuthenticationBackend.filter_users_by_claims`:
match on `preferred_username` first, then case-insensitive email (more than one email match is refused by mozilla).
Found or newly created, the user goes through `create_or_update`, which rewrites the local row from the claims every login:

1. Username, email (validated, `SuspiciousOperation` if invalid), first/last name, with control characters stripped.
2. `OIDC_REQUIRED_GROUP` (eg `/variantgrid/shariant_production`) must be in `groups`, otherwise the user is saved
   inactive and shown which Shariant environments their groups *do* allow. One Keycloak realm (`agha`) serves prod,
   test, demo and security; this group is what keeps a test-only account out of production.
3. Groups are `/associations/<org>/<lab>` and `/variantgrid/<role>`. No group under either prefix: inactive, "doesn't
   belong to any labs". Each association adds every prefix level as a Django `Group` (`org`, `org/lab`), matching
   `Lab.group_name` / `Organization.group_name` (`snpdb/models/models.py:Lab`), created on the fly if missing.
   `/variantgrid/admin` sets `is_superuser` and `is_staff`; `bot` and `tester` become the Django groups
   "variantgrid/bot" and "variantgrid/tester". Everyone also gets `PUBLIC_GROUP_NAME` and `LOGGED_IN_USERS_GROUP_NAME`.
4. Under `MAINTENANCE_MODE` only testers and non-bot admins get through (refused, but not saved inactive).
5. Django groups not in that computed set are removed, the Keycloak `sub` is stored on `UserSettingsOverride.oauth_sub`,
   and `snpdb/models/models_user_settings.py:UserSettingsOverride.auto_set_default_lab` drops or fills the default lab.

Refusals go through `oidc_auth/backend.py:VariantGridOIDCAuthenticationBackend._refuse_login`, which returns the user
inactive; mozilla's callback view then refuses it (`login_failure`), so steps 2-4 end at the login page with the message.
Returning `None` would instead fire `user_login_failed`, which django-axes records against the IP with no username. A
first login builds the user unsaved (`create_user`) and saves only once it is accepted, so a refused first login leaves
no row; an existing user refused at step 2 or 3 is saved inactive. A missing `groups` claim counts as no groups; a
missing `preferred_username` is a `SuspiciousOperation`.

### Staying logged in, logging out, stale callbacks

`oidc_auth/session_refresh.py:VariantGridSessionRefresh` re-checks the ID token's expiry with Keycloak on GET page loads,
but not for any path containing `/api/` or for `X-Requested-With: XMLHttpRequest`, so a long-open page's AJAX calls ride
the Django session rather than being bounced into a redirect. Logout posts to `oidc_logout` (the settings menu picks it
when `USE_OIDC`); `oidc_auth/backend.py:provider_logout` builds the Keycloak end-session URL with `id_token_hint`, which
is why `OIDC_STORE_ID_TOKEN = True`. `oidc_auth/oidc_error_handler.py:HandleOIDC400Middleware` turns a 400 from the
token endpoint on `/oidc/` (a browser replaying a used authorization code) into a redirect to `/` for a fresh login.

### API clients and user creation

DRF on Shariant authenticates with `mozilla_django_oidc.contrib.drf.OIDCAuthentication` (bearer access token) then
`SessionAuthentication`. With no `OIDC_DRF_AUTH_BACKEND` set, mozilla uses the same `VariantGridOIDCAuthenticationBackend`,
so every bearer-token request calls Keycloak's userinfo endpoint and re-runs `create_or_update`, group sync included (read
from mozilla-django-oidc 5.0.2, not measured).

Userinfo accepts a token issued to any client in the realm, and one realm serves every Shariant environment, so
`oidc_auth/backend.py:VariantGridOIDCAuthenticationBackend._verify_bearer_token_client` (called from `get_userinfo` when
there is no ID token, i.e. the DRF path) also requires the token's `azp` to be `OIDC_RP_CLIENT_ID` or one of
`OIDC_API_EXTRA_CLIENT_IDS`, or its `aud` to include `OIDC_RP_CLIENT_ID` (a Keycloak audience mapper); otherwise 401, with
a `report_message` naming the client. The claims are read without a signature check because userinfo has just validated
that exact token. The DRF path has no request, so `_refuse_login` raises `PermissionDenied` (403) rather than returning an
inactive user, which DRF's `IsAuthenticated` would accept.

Users are created in Keycloak, not Django: superusers use `variantgrid/views.py:keycloak_admin` (`/system/keycloak_admin`),
which calls `snpdb/keycloak.py:Keycloak.add_user` with username = email, the lab's two `/associations/` groups plus
`OIDC_REQUIRED_GROUP`, and a password-reset email. The local row appears on first login. `snpdb/keycloak.py:Keycloak`
(admin REST API, credentials in `KEYCLOAK_SYNC_DETAILS`) also sends password resets from user settings and answers
`oidc_auth/signals/keycloak_uptime_check.py:keycloak_uptime_check` for the uptime page.

## Why it is shaped this way

- **Keycloak is authoritative, and trusted fully.** Admin rights, group membership and `is_active` are recomputed on
  every login, and groups are created from whatever paths Keycloak sends, because lab groups are often set up in
  Keycloak before the lab exists in the app (comment in `create_or_update`). Anyone who can edit groups in the realm can
  make a Shariant superuser.
- **`ModelBackend` stays in `AUTHENTICATION_BACKENDS`**, and `/accounts/` and `/admin/` login still take passwords, so
  a local account with a usable password is a way in when Keycloak is down. OIDC-created users have no usable password.
  django-axes only covers these password logins; OIDC brute-force protection is Keycloak's.
- **Deactivating rather than refusing** is how an environment says no: the settings comment ("login failure is generally
  user is inactive, which is how prod distinguishes between prod and test logins") relies on mozilla denying inactive
  users, and the attempt is reported with `report_message`.

## History

`oidc_auth` began as an app named `auth` (the "blank slate" 2020 import; a commented `OIDC_DRF_AUTH_BACKEND` in
`shariantcommon.py` still names `auth.backend`, and `apps.py` has a TODO to rename the config). June 2026 moved to the
Keycloak end-session logout and hardened the backend (claim sanitising, email validation, PKCE); September 2026
(`18a86999e`) narrowed `HandleOIDC400Middleware` from swallowing any error on `/oidc/` to the provider's 400 only, and
built the wrong-environment message with `format_html` so the email is escaped. Both have tests in
`oidc_auth/tests/test_oidc_hardening.py`.

## Traps

- **Refusals must not return `None` or an inactive user to DRF.** Any new refusal in `create_or_update` goes through
  `_refuse_login`: `self.request` is unset on the bearer path, and an inactive user returned there is authenticated.
- **Groups granted in Django admin and a locally granted `is_superuser` are silently undone** at that user's next login
  (deliberate: Keycloak is authoritative).
- **Adding an API client in Keycloak needs a setting.** A new client that fetches tokens for the API must be listed in
  `OIDC_API_EXTRA_CLIENT_IDS` (or given an audience mapper to this deployment's client), or every call is a 401.
- **The required group is also a lab-less pass.** `OIDC_REQUIRED_GROUP` is itself under `/variantgrid/`, so a user with
  only that group is not refused as "doesn't belong to any labs"; they log in with no lab.
- `oidc_auth/backend.py:provider_logout` sends a password login (`/accounts/`, `/admin/`), which has no `oidc_id_token`,
  straight to `LOGOUT_REDIRECT_URL`; only OIDC sessions go through Keycloak's end-session endpoint.
- `oidc_auth/session_refresh.py:VariantGridSessionRefresh` skips `/api/` and `/<app>/api/` paths (DRF), not any path
  containing `/api/`.
