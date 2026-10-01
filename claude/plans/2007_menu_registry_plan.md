# #2007 / #695 — Menus as data, per-user page frame by AJAX, cacheable public pages

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-30; revised by Claude Opus 5.5 (claude-opus-5-5), 2026-10-01
Status: in progress - steps 1-3 (menus as data and the page frame endpoint, #2007) done; step 4 (#695) not started

[#2007](https://github.com/SACGF/variantgrid/issues/2007) (sub-menu rework: menus as configuration, not templates) and
[#695](https://github.com/SACGF/variantgrid/issues/695) (anonymous browsing of gene and variant pages). One plan because
the second depends on the first: a page can only be served from a shared cache once nothing in it depends on who is
asking, and the menu is the first of those things.

## Where it stands

Menus are data in `variantgrid/menus.py:MENUS` (`settings.MENUS`; the classes are `uicore/menus.py`); the per-area
menu bars, the wrapper templates and the context-processor `menu_*_base` variables are gone, and no page template
chooses a menu. The menus and the rest of the per-user page frame come from `uicore/views/page_frame_view.py:page_frame`
after the page loads, so `base.html` no longer reads the user. The how-to is in `uicore/AGENTS.md` and
`claude/research/uicore.md`.

The registry is in-house rather than django-sitetree: sitetree ships models and migrations even for trees declared in
code, matches the current item by exact path (we match by url name), and our admin / POST / external items would have
needed its templates and access hooks overridden. What we kept from it is the shape: detail pages belong to a menu
(`pages`) without being listed, which is what fixes highlighting and makes breadcrumbs a small addition later.

For #695: anything in a response that reads the user forces `Vary: Cookie` (reading `request.session`, which
`request.user` does lazily, is what sets it), and `cache_page` then keys per session, so a shared cache of an expensive
page needs a body that never reads the user. `base.html` is now there; the page templates themselves are the remaining
step, along with `snpdb/processors.py:settings_context_processor`, which reads `request.user` for `somalier_enabled`
when Somalier is on. `analysis/views/views_grid.py` already pairs `cache_page` with `vary_on_cookie` for exactly this reason.

## Data

The frame payload, `uicore/page_frame.py:PageFrame` (JSON, one request per page load):

```python
@dataclass
class PageFrame:
    menu_main_html: str
    menu_sub_html: str
    user_html: str          # inbox link and username / avatar title, '' for anonymous
    username: str           # '' for anonymous
    site_messages_html: str
    messages_html: str      # Django messages, consumed by this request
    rollbar_person: dict    # {} for anonymous
```

## Design

### Page frame endpoint (done)

`GET /uicore/page_frame?url_name=<page url name>` (`uicore/views/page_frame_view.py:page_frame`, `@login_not_required`,
`never_cache`). The menus are `uicore/page_frame.py:menu_html(url_name, MenuRole)`, memoised per process because they
depend only on code and settings; the rest is built per request. Anonymous users get empty menus and no site
messages - a `guest` flag on `MenuItem` and `Menu` comes with the first public page. `loadPageFrame` in `global.js`
starts the request from the head and fills the placeholders on document ready.

### Public pages

- The view takes `@login_not_required`. `LoginRequiredMiddleware` then skips the user check, so the session is never
  read for anonymous requests.
- A `public_page_cache(timeout)` decorator in `library/django_utils/` caches on
  `(path, genome build, PUBLIC_PAGE_VERSION)`, ignoring `Vary`. The build is part of the key because gene and variant
  pages without a build in the URL resolve it through `snpdb/genome_build_manager.py:GenomeBuildManager` (URL, then the
  user's `UserSettings.default_genome_build`, then the first annotated build): anonymous users get the default build's
  entry and a logged-in user the entry for their own build, and the same body serves everyone on that build. Django's own `cache_page` is deliberately not used: any stray session read would silently turn one entry into
  one per session, and a fixed key makes that a visible bug (a logged-in user seeing anonymous content) rather than a
  silent cache miss. The decorator refuses to store a response that carries `Vary: Cookie`, and logs it.
- `PUBLIC_PAGE_VERSION` is a Redis counter bumped by the annotation-version switch and by classification publishing
  (or by a management command); the popular-pages warmer (a celery task on the `annotation` queue, top-N genes and
  alleles by hit count, run after each bump) refills them.
- The cached body must contain no `{% csrf_token %}` (the tag also sets `Vary: Cookie`); forms on public pages POST via
  JS, which already sends the token from the cookie (`variantgrid/static_files/default_static/js/global.js`).
- The cached body is the same for every viewer on that build; anything that depends on the user is an AJAX fragment
  loaded after the page, with its own per-user or per-lab cache key (or none). Anonymous users get the fragment's
  empty or public form. Lab-scoped content on the variant and allele pages (classifications visible to the user's
  labs, the user's tags) goes the same way, or the public page is a reduced reference template - decide per section
  when porting those two pages.
- No `messages.add_message` in a public view: it writes the session, which sets `Vary: Cookie`, and the decorator then
  refuses to cache. A warning derived from the page's own inputs (path and build) is the same for everyone who gets
  that cache entry, so it is rendered into the body. Django messages queued by an earlier request (after a POST) still
  reach the user through the page frame endpoint.

#### view_gene_symbol

`genes/views/views.py:view_gene_symbol` and `genes/gene_symbol_view_info.py:GeneSymbolViewInfo` are mostly
user-independent, but not entirely:

| Part | Depends on | Public page |
|---|---|---|
| gene version, `genome_build` | the build (URL or user default) | in the cache key |
| `warnings()` (symbol not in the requested build, no genes) | symbol and build only | rendered in the body instead of `messages.add_message` |
| `annotation_description` | `UserSettings.tool_tips` | always rendered; tool tips shown or hidden client-side from the user's setting |
| `gene_in_gene_lists` | `GeneList.filter_for_user` | fragment (anonymous: public gene lists only, or the section is hidden) |
| `unmatched_classifications` | `ClassificationModification.filter_for_user` | fragment |
| `classifications` | `filter_for_user` | already AJAX; the endpoint needs an anonymous path (shared classifications only) |
| `has_variants`, `has_samples_in_other_builds`, `has_gene_coverage` | global | in the body - but see Open questions |
| grid columns ("from your User Settings") | the user's column settings | the grids already load by AJAX; anonymous gets the default columns |

### Order of work

1. ~~Menus as data: `variantgrid/menus.py`, rendered server-side.~~ Done.
2. ~~Port every page, delete the menu bars, wrappers and context-processor variables.~~ Done.
3. ~~Page frame endpoint: move the menu and the rest of the per-user page frame into it; strip `base.html`.~~ Done.
4. `public_page_cache`, the version counter and the warmer; open view_gene_symbol first (its user-dependent parts are
   the few listed above), then view_allele and view_variant once their lab-scoped sections are fragments.

### Tests worth keeping

- `public_page_cache` refuses a response with `Vary: Cookie`, two anonymous requests for the same path hit the cache
  once, and the same path on two genome builds is two entries.

## Open questions

- Whether the public variant page is the full template with per-lab fragments or a reduced reference page. Decide when
  step 4 starts, after step 3 has shown how much of the page is already user-independent.
- Whether anonymous users may see the instance-wide flags on the gene page (`has_variants`: observed, tagged,
  classified and ClinVar variants exist for the gene here; `has_samples_in_other_builds`). They are the same for every
  user, so they cache, but they say what this deployment holds - hide the observed/tagged ones and their graphs from
  anonymous users, or accept it for deployments that open gene pages.
