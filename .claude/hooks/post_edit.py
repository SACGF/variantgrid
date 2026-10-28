"""
PostToolUse hook for Edit/Write: the fastest feedback loop.

  *.py    -> ruff check on that file; findings go to stderr with exit 2 so the agent sees them
  *.js    -> eslint on that file when node is installed (never on the vendored js/lib/)
  *.scss  -> remind that the compiled .css needs the same change by hand (AGENTS.md "SCSS / CSS")
  *.md    -> vg docs check on that file when it is one of the agent docs, so a dead citation is caught
             as it is written rather than by the Agent maps job on push
  *.html, *.scss, our *.js -> vg css unused, as a warning: a class this edit stopped naming (or a rule it added
             that nothing names) is reported while the agent can still act on it. Deliberately not a CI check -
             the scan cannot see classes built from data or added by JS, so a human or agent decides
  *.html with a <script>, our *.js -> once per file per session, a reminder to load a page that runs it with
             `vg page <url> --js-writes-to-database`, naming candidate URLs from claude/maps/urls.md
"""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = os.environ.get("CLAUDE_PROJECT_DIR") or os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)

from library.vg.css import find_unused, render_unused_report  # noqa: E402
from library.vg.docs import check_docs, doc_files, render_report  # noqa: E402

CSS_SOURCE_SUFFIXES = (".html", ".scss", ".js")
URL_MAP = Path(REPO) / "claude" / "maps" / "urls.md"
JS_REMINDED_DIR = Path(tempfile.gettempdir()) / "vg_js_check_reminders"
TEMPLATE_NAME_RE = re.compile(r"/templates/(.+\.html)$")
STATIC_NAME_RE = re.compile(r"/[a-z]+_static/(.+\.js)$")
MAX_SUGGESTED_URLS = 3
# More matching URLs than this means a base template: any page will do
EVERY_PAGE_URLS = 30


def _is_our_static_js(relative: str) -> bool:
    return relative.endswith(".js") and "static_files" in relative and "/js/lib/" not in relative


def _css_unused_warning(relative: str) -> str:
    if not relative.endswith(CSS_SOURCE_SUFFIXES) or (relative.endswith(".js") and not _is_our_static_js(relative)):
        return ""
    report = find_unused()
    if report.ok:
        return ""
    return (f"Warning (not blocking) - vg css unused after editing {relative}:\n{render_unused_report(report)}\n"
            "Delete a rule nothing uses, or if the class is built from data or added by JS, add it to "
            "DATA_DERIVED / LIBRARY_ADDED in library/vg/css.py.\n")


def _templates_naming(needle: str) -> set[str]:
    """ Templates whose source mentions needle (an {% include %}/{% extends %} name, or a static path) """
    names = set()
    for template in Path(REPO).glob("*/templates/**/*.html"):
        if needle in template.read_text(errors="replace"):
            if m := TEMPLATE_NAME_RE.search(template.as_posix()):
                names.add(m[1])
    return names


def _candidate_urls(relative: str) -> list[str]:
    """ URL patterns whose view renders the edited template, or a template one include/extends up from it """
    if not URL_MAP.exists():
        return []
    if m := TEMPLATE_NAME_RE.search("/" + relative):
        templates = {m[1]} | _templates_naming(m[1])
    elif m := STATIC_NAME_RE.search("/" + relative):
        templates = _templates_naming(m[1])
        templates |= {includer for t in templates for includer in _templates_naming(t)}
    else:
        return []
    urls = []
    for row in URL_MAP.read_text().splitlines():
        cells = [c.strip() for c in row.split("|")]
        if len(cells) > 4 and cells[2].startswith("`") and templates & set(cells[4].split(", ")):
            urls.append("/" + cells[2].strip("`"))
    if len(urls) > EVERY_PAGE_URLS:
        return ["/variantopedia/dashboard (it loads on most pages)"]
    return urls[:MAX_SUGGESTED_URLS]


def _js_check_reminder(relative: str, path: str, session_id: str) -> str:
    """ Once per file per session, so repeated edits to one template don't repeat it """
    if not (_is_our_static_js(relative) or (relative.endswith(".html") and "<script" in Path(path).read_text(errors="replace"))):
        return ""
    marker = JS_REMINDED_DIR / session_id / re.sub(r"[^A-Za-z0-9]+", "_", relative)
    if marker.exists():
        return ""
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.touch()
    urls = _candidate_urls(relative)
    where = f" - pages that render it: {', '.join(urls)}" if urls else ""
    static = " A static JS edit needs collectstatic first." if relative.endswith(".js") else ""
    return (f"Reminder (not blocking): {relative} has page JS. Before reporting back, load a page that runs it with "
            f"`python3 manage.py vg page <url> --js-writes-to-database`{where}. Pick a URL whose data reaches the code "
            f"you changed. It WRITES to the database (not rolled back): fine on a developer box, ask first on a shared "
            f"test box, refused elsewhere.{static}\n")


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except ValueError:
        return 0
    path = (payload.get("tool_input") or {}).get("file_path") or ""
    if not path or not os.path.exists(path):
        return 0
    relative = os.path.relpath(path, REPO)
    if relative.startswith("..") or "/migrations/" in relative:
        return 0

    if relative.endswith(".py"):
        ruff = os.path.join(REPO, ".venv", "bin", "ruff")
        if not os.path.exists(ruff):
            ruff = shutil.which("ruff")
        if ruff:
            result = subprocess.run([ruff, "check", "--quiet", "--output-format", "concise", path],
                                    capture_output=True, text=True, check=False)
            if result.returncode:
                sys.stderr.write(f"ruff on {relative}:\n{result.stdout}{result.stderr}")
                return 2
    elif _is_our_static_js(relative):
        npx = shutil.which("npx")
        if npx:
            result = subprocess.run([npx, "--no-install", "eslint", path], capture_output=True, text=True, check=False, cwd=REPO)
            if result.returncode and result.stdout.strip():
                sys.stderr.write(f"eslint on {relative}:\n{result.stdout}")
                return 2
    elif relative.endswith(".md") and Path(path).resolve() in set(doc_files()):
        report = check_docs([path])
        if not report.ok:
            sys.stderr.write(f"vg docs check on {relative}:\n{render_report(report)}\n")
            return 2
    elif relative.endswith(".scss"):
        css = relative[:-5] + ".css"
        sys.stderr.write(f"{relative} is compiled by a PyCharm watcher; hand-apply the same minimal change to {css} "
                         f"(match its formatting, leave the .css.map alone).\n{_css_unused_warning(relative)}")
        return 2
    warning = _css_unused_warning(relative) + _js_check_reminder(relative, path, payload.get("session_id") or "no-session")
    if warning:
        sys.stderr.write(warning)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
