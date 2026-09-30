# #2007 / #695 — Menus as data, user chrome out of the page, cacheable public pages

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-30; revised by Claude Opus 5.5 (claude-opus-5-5), 2026-09-30
Status: in progress - steps 1-2 (menus as data, #2007) done; steps 3-4 (#695) not started

[#2007](https://github.com/SACGF/variantgrid/issues/2007) (sub-menu rework: menus as configuration, not templates) and
[#695](https://github.com/SACGF/variantgrid/issues/695) (anonymous browsing of gene and variant pages). One plan because
the second depends on the first: a page can only be served from a shared cache once nothing in it depends on who is
asking, and the menu is the first of those things.

## Where it stands

Menus are data in `variantgrid/menus.py:MENUS` (`settings.MENUS`; the classes are `uicore/menus.py`), rendered server-side by `uicore/templatetags/ui_menus.py:menu_bar_main` and
`menu_bar_sub` from the request's url name; the per-area menu bars, the wrapper templates and the context-processor
`menu_*_base` variables are gone, and no page template chooses a menu. The how-to is in `uicore/AGENTS.md` and
`claude/research/uicore.md`.

The registry is in-house rather than django-sitetree: sitetree ships models and migrations even for trees declared in
code, matches the current item by exact path (we match by url name), and our admin / POST / external items would have
needed its templates and access hooks overridden. What we kept from it is the shape: detail pages belong to a menu
(`pages`) without being listed, which is what fixes highlighting and makes breadcrumbs a small addition later.

For #695, the base page (`uicore/templates/uicore/page/base.html`) still reads the user in five places: the Rollbar
config block (id, username, email), the inbox link and count, the username and avatar, site messages and Django
messages, and the menu (admin-only items). Any of those forces `Vary: Cookie` (reading `request.session`, which
`request.user` does lazily, is what sets it), and `cache_page` then keys per session, so a shared cache of an expensive
page is impossible while they are in the body. `analysis/views/views_grid.py` already pairs `cache_page` with
`vary_on_cookie` for exactly this reason.

## Data

The chrome payload (JSON, one request per page load):

```python
@dataclass
class Chrome:
    menu_html: str          # top bar + sub-menu for the requested path, already rendered
    username: str           # '' for anonymous
    avatar_html: str
    inbox_unread: int
    site_messages_html: str
    messages_html: str      # Django messages, consumed by this request
    rollbar_person: dict    # {} for anonymous
```

## Design

### Chrome endpoint

`GET /uicore/chrome/?path=<page path>` in a new uicore urls module, `@login_not_required` (it decides what to show from
`request.user` itself). It resolves `path` to a url name, renders `menu_bar_main` / `menu_bar_sub` for it and fills the
rest of `Chrome`. Cache key: `(role bucket, url_name, CACHE_VERSION)`, where role bucket is anonymous / user /
superuser (the only thing `uicore/menus.py:MenuItem` visibility depends on besides settings). Anonymous users see only
items marked for guests - a `guest` flag on `MenuItem` and `Menu`, added with the first public page. Django messages
and the inbox count are per user and excluded from the cached part. `global.js` fetches it on `DOMContentLoaded` and
fills the top bar, side bar and the navbar right-hand side; the navbar reserves its height so nothing shifts.

`base.html` then loses `{% menu_bar_main %}`, `{% block submenu %}`, the user block and `{% site_messages %}`.

### Public pages

- The view takes `@login_not_required`. `LoginRequiredMiddleware` then skips the user check, so the session is never
  read for anonymous requests.
- A `public_page_cache(timeout)` decorator in `library/django_utils/` caches on `(path, PUBLIC_PAGE_VERSION)`, ignoring
  `Vary`. Django's own `cache_page` is deliberately not used: any stray session read would silently turn one entry into
  one per session, and a fixed key makes that a visible bug (a logged-in user seeing anonymous content) rather than a
  silent cache miss. The decorator refuses to store a response that carries `Vary: Cookie`, and logs it.
- `PUBLIC_PAGE_VERSION` is a Redis counter bumped by the annotation-version switch and by classification publishing
  (or by a management command); the popular-pages warmer (a celery task on the `annotation` queue, top-N genes and
  alleles by hit count, run after each bump) refills them.
- The cached body must contain no `{% csrf_token %}` (the tag also sets `Vary: Cookie`); forms on public pages POST via
  JS, which already sends the token from the cookie (`variantgrid/static_files/default_static/js/global.js`).
- Lab-scoped content on the variant and allele pages (classifications visible to the user's labs, the user's tags)
  becomes AJAX fragments with their own per-lab cache keys, or the public page is a reduced reference template. Decide
  per section when porting those two pages; everything else on them is user-independent.

### Order of work

1. ~~Menus as data: `variantgrid/menus.py`, rendered server-side.~~ Done.
2. ~~Port every page, delete the menu bars, wrappers and context-processor variables.~~ Done.
3. Chrome endpoint: move the menu and the rest of the user chrome into it; strip `base.html`.
4. `public_page_cache`, the version counter and the warmer; open view_gene_symbol first (no lab-scoped content), then
   view_allele and view_variant once their lab-scoped sections are fragments.

### Tests worth keeping

- The chrome endpoint's cached part is the same for two users in the same role bucket and differs between anonymous
  and superuser.
- `public_page_cache` refuses a response with `Vary: Cookie`, and two anonymous requests for the same path hit the
  cache once.

## Open questions

- Whether the public variant page is the full template with per-lab fragments or a reduced reference page. Decide when
  step 4 starts, after step 3 has shown how much of the page is already user-independent.
