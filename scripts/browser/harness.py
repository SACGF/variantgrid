"""
Shared plumbing for shot.py and jsq.py: where the harness lives, the agent session cookie, a chromium launch
that can see the side-loaded libraries, and the route that turns a served page into standards mode.
"""
import os
from pathlib import Path

HOME = Path(os.environ.get("VG_BROWSER_HOME", Path.home() / ".cache" / "claude-browser"))
BASE = os.environ.get("VG_BROWSER_BASE", "http://localhost:8765")
LIBS = HOME / "libs/usr/lib/x86_64-linux-gnu"
VIEWPORT = {"width": 1280, "height": 900}


def session_key() -> str:
    if key := os.environ.get("VG_BROWSER_SESSION"):  # vg page --js-writes-to-database passes its own session
        return key
    path = HOME / "sessionid"
    if not path.exists():
        raise SystemExit(f"{path} missing: run  python3 manage.py shell -c \"$(cat scripts/browser/agent_session.py)\"")
    return path.read_text().strip()


def launch(playwright):
    return playwright.chromium.launch(env={**os.environ, "LD_LIBRARY_PATH": str(LIBS)})


def new_context(browser, mode: str, patches=(), viewport=None):
    """ A logged-in context. mode 'standards' prepends <!DOCTYPE html> to every document; 'quirks' serves pages
        as they are. patches are (old, new) text replacements applied to every HTML response, documents and
        ajax fragments alike, so a CSS or markup fix can be tried without editing a template. """
    ctx = browser.new_context(viewport=viewport or VIEWPORT)
    ctx.add_cookies([{"name": "sessionid", "value": session_key(), "url": BASE}])

    def rewrite(route, request):
        resp = route.fetch()
        if "text/html" not in resp.headers.get("content-type", ""):
            return route.fulfill(response=resp)
        body = resp.text()
        for old, new in patches:
            body = body.replace(old, new)
        if request.resource_type == "document" and mode == "standards" \
                and not body.lstrip().lower().startswith("<!doctype"):
            body = "<!DOCTYPE html>\n" + body
        route.fulfill(response=resp, body=body)

    ctx.route(lambda url: True, lambda route, request: rewrite(route, request)
              if request.resource_type in ("document", "xhr", "fetch") else route.continue_())
    return ctx
