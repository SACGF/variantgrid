# Headless browser harness

Screenshots and layout measurements of real pages, for an agent (or a developer) checking template and CSS
changes: `vg page` gives the HTML and the query count, this gives the rendered boxes and a PNG to look at.
Everything runs from a working-tree `runserver`, logged in as the `claude_agent` user, so it never touches
gunicorn or edits a file.

The harness lives outside the repo in `$VG_BROWSER_HOME` (default `~/.cache/claude-browser`): a uv venv with
playwright and, on boxes without root, the chromium shared libraries side-loaded from Ubuntu packages.
Nothing is added to the project venv or `requirements*.in`.

## One-off setup

```bash
scripts/browser/setup.sh            # venv + playwright + chromium; side-loads missing libs via apt-get download
```

## Each session

```bash
python3 manage.py runserver localhost:8765 --insecure --noreload --skip-checks &   # working tree, static from STATIC_ROOT
python3 manage.py shell -c "$(cat scripts/browser/agent_session.py)"               # writes $VG_BROWSER_HOME/sessionid
PY=~/.cache/claude-browser/pw/bin/python
$PY scripts/browser/shot.py --measure ".hotspot-graph" /genes/view_gene_symbol/BRCA1/GRCh37
$PY scripts/browser/jsq.py /snpdb/data "() => document.querySelector('#samples-datatable').offsetHeight"
kill %1                                                                            # stop the runserver when done
```

`shot.py` renders each URL in the modes asked for (`--modes quirks,standards`; standards prepends
`<!DOCTYPE html>` to every document response through a playwright route, so both modes come from the same
server), saves `shots/<slug>__<mode>_<vh>.png`, and prints status, `document.compatMode`, page height and the
boxes of any `--measure` selectors. Options worth knowing:

- `--patch 'OLD=>NEW'` rewrites served HTML (documents and ajax fragments) before the browser sees it - try a
  CSS fix without editing a template.
- `--click SELECTOR` opens a tab or node before measuring.
- `--unroll` opens the `.main-content` scroller so a full-page screenshot shows the whole page; `--clip SELECTOR`
  screenshots just the region around one element.
- `--dump` records every element's box; `compare.py shots/results.json` then lists what collapsed, appeared,
  resized or changed font size between the two modes.

`jsq.py URL 'JS'` evaluates one expression (an `async () =>` is fine) in each mode and prints both results.

## Gotchas

- `--skip-checks`: the project's `variantgrid.E001` check refuses `runserver` while `SESSION_COOKIE_SECURE` is
  on. Chromium treats `http://localhost` as a secure context and sends Secure cookies anyway, so the check is
  the only thing in the way.
- Sessions are in the cache backend (`SESSION_ENGINE`), not `django_session`: `agent_session.py` uses the
  configured engine, a `SessionStore` from the db backend is never found.
- `playwright install` prunes browser builds other versions used; pin the playwright version in `setup.sh`
  when changing it.
- Pages that render as bare fragments (no navbar, `$ is not defined` in `errors`) are ajax tab bodies, not
  pages: find the page that includes them.
