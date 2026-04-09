# VariantGrid Flags App — Reference Document

## Purpose and Overview

The `flags` app is a sophisticated flagging system that allows users to raise, resolve, and comment on issues ("flags") related to different types of objects in VariantGrid (Alleles, Classifications, Clinical Contexts, etc.). It provides a flexible framework for tracking problems, tracking their lifecycle, and managing permissions around who can view, raise, and resolve different flag types.

**Key Design:** Signal-based architecture for decoupled integration; four-tier permission hierarchy; audit trail via comments.

---

## Models

### FlagResolution
Possible resolution states for flags.
- Fields: id (PK, text, e.g. 'open', 'closed', 'rejected'), label, description, status (FlagStatus: OPEN='O', CLOSED='C', REJECTED='R')
- TimeStampedModel

### FlagTypeContext
Groups flag types by domain (e.g., 'classification', 'allele', 'clinical_context').
- Fields: id (PK, text), label

### FlagType
Defines a specific type of flag that can be raised.
- Fields: id (PK, text, e.g. 'classification_discordant'), context FK, label, description, help_text (only shown to users who can edit), raise_permission, permission (FlagPermissionLevel), only_one (bool, one open flag per collection), attributes (JSON), comments_enabled (bool), importance (int)
- Manager: ObjectManagerCachingImmutable (cached, immutable)
- Methods: resolution_for_status(status), default_resolution()

### FlagTypeResolution
Junction table: FlagType ↔ FlagResolution (not all resolutions valid for all types).
- Unique together: (flag_type, resolution)

### Flag
An actual flagged issue instance attached to a FlagCollection.
- Fields: collection FK, flag_type FK, user FK, resolution FK (PROTECT), user_private (bool), data (JSON)
- Key method: `flag_action(resolution=None, user=None, comment=None, permission_check=True, first_comment=False)` — atomic, creates FlagComment, updates resolution, emits flag_comment_action signal

### FlagComment
Audit trail and discussion on a flag.
- Fields: flag FK, user FK, text, resolution FK (nullable — if this comment changed the resolution)
- Static method: `FlagComment.last(flag)` — most recent comment

### FlagCollection
Container for all flags related to a specific object.
- Base: GuardianPermissionsMixin
- Fields: context FK (FlagTypeContext, PROTECT)
- Custom permission: admin_flagcollection
- Key methods:
  - `permission_level(user)` → FlagPermissionLevel (returns SYSTEM for admin_bot, ADMIN for superuser, else delegates to source_object.flag_user_permission())
  - `flags(user=None, only_open=False)` → QuerySet[Flag] (filters private flags)
  - `add_flag(flag_type, user=None, comment=None, ...)` → Flag
  - `get_or_create_open_flag_of_type(flag_type, ...)` → (Flag, bool) — handles reopen, old_data, close_other_data
  - `close_open_flags_of_type(flag_type, ...)` → int
  - `ensure_resolution(flag_type, resolution, comment=None)` — idempotent
  - `filter_for_flags(qs, flag_types=None, open_only=True)` — static, filter objects with flags
  - `source_object` property — lazy-loads the related object
  - `extra_info` property — populated via flag_collection_extra_info_signal

### FlagsMixin (Abstract)
Mixin for models that can have flags.
- Fields: flag_collection FK (FlagCollection, nullable)
- Abstract methods: flag_type_context() → FlagTypeContext
- Methods:
  - `flag_user_permission(user)` → FlagPermissionLevel (default: OWNER if can_write, USERS if can_view)
  - `flags_of_type(flag_type)` → QuerySet
  - `has_open_flags` property
  - `has_open_flag_with_attribute(attribute, value)` → bool
  - `close_open_flags_of_type(flag_type, comment=None, ...)`
  - `has_flag_activity(since=None)` → bool
  - `flag_collection_safe` property — lazy-creates FlagCollection if needed
- Known implementations: **Allele**, **Classification**, **ClinicalContext**

### FlagInfos
Batch helper for loading and manipulating flags with related data.
- Methods: set_extra_info(pk, extra_info, source_object), extra_flag_info(flag)
- Attributes: flag_collections, flags_for_collection, ids

---

## Permission System

### FlagPermissionLevel Enum (totally ordered)
| Level | Value | Meaning |
|-------|-------|---------|
| NO_PERM | 0 | Cannot see flags |
| USERS | 1 ('U') | Regular user — can see public flags, raise user-level flags |
| OWNER | 2 ('O') | Object owner — can manage object-related flags |
| ADMIN | 3 ('X') | Admin — can view all flags including private ones |
| SYSTEM | 4 ('A') | Bot — can perform any action without checks |

Each FlagType has `raise_permission` (minimum to create) and `permission` (minimum to resolve).
Private flags (`user_private=True`) visible only to creator and users with permission >= OWNER.

---

## Signals

### flag_collection_extra_info_signal
- Sender: FlagCollection class
- Args: flag_infos (FlagInfos), user (User)
- Purpose: Allow other apps to attach additional metadata to FlagCollections
- Usage: Classification app adds classification label

### flag_comment_action
- Sender: Flag class
- Args: flag_comment (FlagComment), old_resolution (FlagResolution)
- Purpose: Notify subscribers when flag is changed or commented on
- Usage: Can trigger business logic like updating related objects

---

## REST API

### FlagsView — `/api/flags/<flag_collection_id>`
- Supports comma-separated IDs
- **GET params**: history (include full history), since (changes since timestamp)
- **POST params**: flag_type, comment, user_private, resolution, watch
- Returns: flags, comments, flag types, resolutions, user info

### FlagView — `/api/flag/<flag_id>`
- **GET**: Flag details with full comment history
- **POST params**: resolution (new resolution ID), comment
- Returns: Updated flag data

### FlagHelper
Internal helper that serializes flags and related data to JSON, handles permission filtering, lazy-loads and caches flag types/resolutions/users.

---

## URL Patterns
```python
path('api/flags/<flag_collection_id>', FlagsView.as_view(), name='flags_api'),
path('api/flag/<int:flag_id>', FlagView.as_view(), name='flag_api'),
```

---

## Admin Interface

- **FlagAdmin**: List + search by ID, filter by flag_type/resolution/user; FlagCommentAdminTabular inline
- **FlagCommentAdmin**: id, flag, user, text, resolution, created, modified
- **FlagTypeAdmin**: id, context, label, description, help_text, raise_permission
- **FlagCollectionAdmin**: id, text; FlagInline (read-only); no creation/modification

---

## Business Logic

### Flag Lifecycle
1. FlagsMixin object creates FlagCollection (if needed) via `flag_collection_safe`
2. Raise flag: `flag_collection.add_flag()` or `get_or_create_open_flag_of_type()`
3. Comment: `flag.flag_action(comment="...")`
4. Resolve: `flag.flag_action(resolution=new_resolution)`
5. Closure: Resolution changes to CLOSED or REJECTED status

### only_one=True Flags
Only one open instance of this flag type per collection at a time. `get_or_create_open_flag_of_type()` automatically handles reuse/reopen.

### Data Tracking
Flags store arbitrary JSON `data` — useful with `close_other_data` to close flags with different data while keeping the current one.

---

## Key Integration Points

| App | Integration |
|-----|-------------|
| Classification | Flags for discordance, significance change, transcript version change, unshared state |
| Allele | Flags for variant matching issues, GRCh37/38 mismatches |
| ClinicalContext | Flags for discordance detection |
| Any app | Any model implementing FlagsMixin can have flags |
