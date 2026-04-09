# VariantGrid VCAuth and OIDC Auth Apps — Reference Document

---

## VCAuth App — User Administration

### Purpose and Overview

`vcauth/` is a minimal Django app that customizes Django's authentication and user administration. It extends the default Django User admin interface with custom actions and integrations with the classification system.

### Files

- `admin.py` — Django admin registration
- `user_admin.py` — Custom User admin class

### CustomUserAdmin

Extends Django's `UserAdmin`. Registered by unregistering the default User admin.

**Features:**
- `email_discordance()` action — Send weekly summary emails to selected users
  - Calls `send_summary_email_to_user()` from classification app
  - Provides user feedback messages on success/failure
  - Shows warning if email server issues occur

**Actions:** "Email weekly summary" — bulk action on user list

### Integration Points

- **Django Admin**: Unregisters default User admin, registers CustomUserAdmin
- **Classification App**: Uses `send_summary_email_to_user()` for summaries
- **Email Manager**: Relies on email sending infrastructure

---

## OIDC Auth App — OpenID Connect Authentication

### Purpose and Overview

`oidc_auth/` handles OpenID Connect (OIDC) authentication and session management. Extends mozilla-django-oidc with custom logic for Keycloak integration, group management, and environment-specific access control.

### Backend Authentication (backend.py)

**VariantGridOIDCAuthenticationBackend** — Extends `mozilla_django_oidc.auth.OIDCAuthenticationBackend`

**Key Methods:**
- `filter_users_by_claims()` — Find user by username first, fall back to email (case-insensitive)
- `create_user()`, `update_user()` — Delegate to `create_or_update()`
- `create_or_update()` — Main sync logic:
  - Updates user profile from OIDC claims (username, email, names)
  - Validates `OIDC_REQUIRED_GROUP` if configured
  - Parses nested Keycloak group structure:
    - `/associations/org/lab` → Association groups
    - `/variantgrid/admin`, `bot`, `tester` → Role-based permissions
  - Syncs group membership (add/remove from Django groups)
  - Sets superuser/staff flags based on admin role
  - Stores `oauth_sub` in UserSettingsOverride
  - Handles maintenance mode (only admins/testers allowed)

**Group Processing:**
- Extracts nested path components as individual groups
- Example: `/associations/pathology/lab1` → groups `pathology/lab1` and `pathology`

### Session Management (session_refresh.py)

**VariantGridSessionRefresh** — Extends `mozilla_django_oidc.middleware.SessionRefresh`
- Detects AJAX requests via `HTTP_X_REQUESTED_WITH` header
- Skips session refresh for API endpoints (`/api/*`)
- Prevents unnecessary OIDC refresh for API calls

### Error Handling (oidc_error_handler.py)

**HandleOIDC400Middleware** — Catches OIDC-related errors
- Catches 400 errors on `/oidc/*` paths
- Redirects to home page to force fresh token fetch

### App Configuration (apps.py)

`OIDCAuthConfig` — imports `keycloak_uptime_check` signal handler in `ready()`

---

### User Login Flow
1. User initiates OIDC authentication
2. `VariantGridOIDCAuthenticationBackend` receives OIDC claims
3. Filter by username/email to find or create user
4. Validate required group membership
5. Parse and sync group memberships from Keycloak
6. Set superuser/staff flags
7. Store OAuth sub ID
8. Return authenticated user

---

### Required Settings

| Setting | Purpose |
|---------|---------|
| `OIDC_RP_CLIENT_ID` | OIDC client ID |
| `OIDC_RP_CLIENT_SECRET` | OIDC client secret |
| `OIDC_OP_AUTHORIZATION_ENDPOINT` | Keycloak auth endpoint |
| `OIDC_OP_TOKEN_ENDPOINT` | Keycloak token endpoint |
| `OIDC_OP_USER_ENDPOINT` | Keycloak user info endpoint |
| `OIDC_REQUIRED_GROUP` | Optional group path to enforce |
| `KEYCLOAK_SYNC_DETAILS` | Keycloak API credentials dict |
| `MAINTENANCE_MODE` | Boolean for read-only mode |
| `PUBLIC_GROUP_NAME` | Public group constant |
| `LOGGED_IN_USERS_GROUP_NAME` | Logged-in group constant |

---

### Integration Points

| App | Integration |
|-----|-------------|
| Keycloak | OIDC provider and admin API |
| Django Auth | User model and permissions |
| Library | Keycloak uptime check signal receiver, library.keycloak.Keycloak |
| snpdb | UserSettingsOverride stores OAuth metadata |
| Django Groups | Syncs Keycloak groups to Django groups |

---

## Expression App (Legacy/Unused)

`expression/` is a legacy/dormant Django app related to gene expression data analysis.

- `models.py` — Empty (single newline)
- `migrations/` — 2 migration files; models were created then removed (CuffDiffRecord, cuff_diff_file)
- No active views, URLs, admin configuration, or business logic
- Registered as a Django app (migrations tracked) but no active functionality
- Historical use: possibly Cufflinks gene expression analysis / RNA-seq differential expression
