# #2007 / #695 — Menus as data (django-sitetree), user chrome out of the page, cacheable public pages

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-30
Status: draft

[#2007](https://github.com/SACGF/variantgrid/issues/2007) (sub-menu rework: menus as configuration, not templates) and
[#695](https://github.com/SACGF/variantgrid/issues/695) (anonymous browsing of gene and variant pages). One plan because
the second depends on the first: a page can only be served from a shared cache once nothing in it depends on who is
asking, and the menu is the first of those things.

## The problem

Menus are templates, and every page picks its own sub-menu:

- 11 sub-menu bars (`uicore/templates/uicore/menus/menu_bar_data.html` and siblings), each rendered by a one-line
  inclusion tag in `uicore/templatetags/ui_menu_bars.py`.
- 10 wrapper templates (`snpdb/templates/snpdb/menu/menu_data_base.html` and siblings, plus
  `seqauto/templates/seqauto/menu_sequencing_data_base.html`) whose only job is to fill `{% block submenu %}`.
- 163 page templates choose a sub-menu: 95 inline (`{% block submenu %}{% menu_bar_x %}{% endblock %}`) and 68 by
  `{% extends menu_x_base %}`, through ten context variables that `snpdb/processors.py:settings_context_processor`
  injects for a clinician-restricted mode that no longer exists.
- Two "active item" heuristics that disagree. `uicore/templatetags/ui_menus.py:menu_top` matches URL path prefixes
  hard-coded in `uicore/templates/uicore/menus/menu_bar_main.html` (`app_name='patients|snpdb/cohorts|pedigree|...'`);
  `uicore/templatetags/ui_menus.py:menu_item` matches url_name, then route prefix, then an `other_urls` list. A page
  whose URL is not itself a menu item (`variantopedia/views_allele.py:view_allele`, view_vcf, view_classification) is
  highlighted only if someone remembered the override in its template.
- Per-deployment visibility is already data: `variantgrid/perm_path.py:get_visible_url_names` reads
  `URLS_NAME_REGISTER`. Only the structure is not.

For #695, the base page (`uicore/templates/uicore/page/base.html`) reads the user in five places: the Rollbar config
block (id, username, email), the inbox link and count, the username and avatar, site messages and Django messages, and
the menu. Any of those forces `Vary: Cookie` (reading `request.session`, which `request.user` does lazily, is what sets
it), and `cache_page` then keys per session, so a shared cache of an expensive page is impossible while they are in the
body. `analysis/views/views_grid.py` already pairs `cache_page` with `vary_on_cookie` for exactly this reason.

## Decision

1. **django-sitetree**, trees declared in code, `SITETREE_DYNAMIC_ONLY = True`. Read in full at 1.18.1 (March 2026,
   Python 3.10+). It gives: `tree()`/`item()` per app in a sitetrees module; items by URL name, with args
   (`'view_allele allele.pk'`); `in_menu=False` items so a detail page has a place in the tree without being a menu
   entry (the proper fix for the highlight problem); `access_guest`, `access_loggedin`, `access_by_perms` and a per-item
   `access_check(context)` hook (the anonymous story of #695, and where `URLS_NAME_REGISTER` and `admin_only` go);
   breadcrumbs and page-title tags; `register_dynamic_trees(..., target_tree_alias, parent_tree_item_alias)` so the
   sapath and shariant apps graft items into our menus without template overrides; tree structure cached in Django
   cache. django-simple-menu (Jazzband, 330 lines, classifiers to Django 4.2) was the alternative: it selects the
   current item by regex on `request.path`, the heuristic we are trying to leave, and offers little we would not write
   ourselves.
2. **One chrome endpoint.** Everything in the base page that depends on the user (menu, username, avatar, inbox, site
   and Django messages, Rollbar person) is fetched by one AJAX call after load and cached per role bucket and url name.
   The page body renders without touching `request.user` or the session.
3. **Public pages** opt in with `@login_not_required` and a page cache keyed on path with a version prefix, warmed for
   the popular genes and variants.

## Data

No database models. The menu is a Python data structure:

```python
# <app>/sitetrees.py, autodiscovered by sitetree
sitetrees = (
    tree('main', items=[
        item('Variants', 'variant_tags', alias='variants', access_check=visible, children=[
            item('Liftover', 'liftover_runs', access_check=superuser_only),
            item('Tagged Variants', 'variant_tags', access_check=visible),
            item('Enter Variants', 'manual_variant_entry', access_check=visible),
            item('Variants', 'variants', access_check=visible),
            item('Variant Wiki', 'variant_wiki', access_check=visible),
            # detail pages: in the tree for highlighting and breadcrumbs, never listed
            item('Allele', 'view_allele allele.pk', in_menu=False, access_guest=True, access_loggedin=True),
            item('Variant', 'view_variant variant.pk', in_menu=False, access_guest=True, access_loggedin=True),
        ]),
        ...
    ]),
)
```

`visible` is `lambda context: get_visible_url_names()[item.url_name]` (sitetree passes the item as `context['item']`
- confirm when implementing); `superuser_only` is `visible` plus `request.user.is_superuser`. Per-deployment
differences stay in `URLS_NAME_REGISTER`; the tree is the same on every deployment.

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

### Current-item matching

sitetree resolves the current item by comparing `request.path` with each item's resolved URL, which works for
argument-less items and, for `'view_allele allele.pk'`, only when `allele` is in the template context. The chrome
endpoint has no such context, so a `SITETREE_CLS` subclass overrides `get_tree_current_item` to resolve the requested
path with `django.urls.resolve` and pick the item whose URL name matches. One matcher for everything, no prefix lists.

### Chrome endpoint

`GET /uicore/chrome/?path=<page path>` in a new uicore urls module, `@login_not_required` (it decides what to show from
`request.user` itself). It renders the sitetree menu for `path` as the current page and fills the rest of `Chrome`.
Cache key: `(role bucket, url_name, CACHE_VERSION)`, where role bucket is anonymous / user / superuser plus whatever
the tree's `access_check`s depend on (today: nothing else). Django messages and the inbox count are per user and
excluded from the cached part. `global.js` fetches it on `DOMContentLoaded` and fills `#menu-top`, `#submenu` and the
navbar right-hand side; the navbar reserves its height so nothing shifts.

`base.html` loses `{% menu_bar_main %}`, `{% block submenu %}`, the user block and `{% site_messages %}`. The two
`menu_top` / `menu_item` tags, `ui_menu_bars.py`, the 11 menu bars, the 10 wrapper templates and the
`MENU_BASE_TEMPLATES` loop in the context processor are deleted once every page is ported.

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

1. Add django-sitetree to `requirements.in`, compile, confirm the suite passes against Django 6.1 (its classifiers stop
   before 6; tox tests Django main). Settings: `sitetree` in `INSTALLED_APPS`, `SITETREE_DYNAMIC_ONLY = True`,
   `SITETREE_CLS` pointing at the url_name matcher.
2. Variants tree only: a variantopedia sitetrees module, render it from the chrome endpoint, port the 12 pages that use
   the Variants sub-menu, delete `menu_bar_variants.html` and `menu_variants_base.html`. `vg page` before and after on
   view_allele and variant_tags.
3. Remaining trees app by app (largest first: classifications 41 pages, patients 28, settings 23, data 24,
   sequencing 20, genes 12, analysis 12, annotation 8, tests 6, help 1). Each step deletes its bar and wrapper.
4. Move the rest of the user chrome into the endpoint; strip `base.html`; delete the tags, `ui_menu_bars.py` and the
   context-processor loop.
5. `public_page_cache`, the version counter and the warmer; open view_gene_symbol first (no lab-scoped content), then
   view_allele and view_variant once their lab-scoped sections are fragments.

### Tests worth keeping

- Tree items whose URL name is not in the resolver, or is hidden by `URLS_NAME_REGISTER`, are absent for a user and a
  superuser alike; `admin_only` items are absent for a non-superuser; anonymous sees only `access_guest` items.
- The url_name matcher highlights the Variants tree for `/variantopedia/view_allele/<id>` and nothing for an unknown
  path.
- `public_page_cache` refuses a response with `Vary: Cookie`, and two anonymous requests for the same path hit the
  cache once.

## Open questions

- Whether to keep the `Help` link and bug-report icon in the tree (as `item(... url='https://...')`, sitetree accepts
  absolute URLs) or leave them as static HTML in the navbar. Static is simpler and they never vary.
- Whether the public variant page is the full template with per-lab fragments or a reduced reference page. Decide when
  step 5 starts, after steps 1-4 have shown how much of the page is already user-independent.
