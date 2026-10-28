"""
`vg page --js-writes-to-database`: load one page in headless chromium as `claude_agent` and report what its JavaScript complained
about - uncaught exceptions, console errors and warnings, failed requests and 4xx/5xx responses - with each
inline-script frame mapped back to the template line it came from.

The page is served from this process: a threaded WSGI server on a free localhost port running the working tree,
with static files from STATIC_ROOT (templates name the manifest's hashed files, so a JS or CSS edit shows up
after collectstatic). The browser is scripts/browser/console.py in the harness venv ($VG_BROWSER_HOME, made by
scripts/browser/setup.sh), logged in with a short-lived session. Unlike `vg page`, NOTHING IS ROLLED BACK, hence
the flag's name: the page and the ajax it fires on load write to the database as they would for a user. Celery's send_task, Rollbar,
outgoing POSTs and email are stubbed for the duration. It refuses to run outside a developer box or a listed
shared test deployment (`check_box_allows_writes`), so production is never written to.

Entry points: `run_page_js` (returns a JsReport) and `render_js_report`.
"""
import json
import os
import re
import subprocess
import threading
from collections import Counter
from contextlib import contextmanager
from dataclasses import dataclass, field
from importlib import import_module
from pathlib import Path
from urllib.parse import urlsplit

from django.conf import settings
from django.contrib.staticfiles.handlers import StaticFilesHandler
from django.core.handlers.wsgi import WSGIHandler
from django.core.servers.basehttp import ThreadedWSGIServer, WSGIRequestHandler
from django.template import TemplateDoesNotExist
from django.template.loader import get_template
from django.views.static import serve as static_serve

from library.vg.page import AGENT_USERNAME, agent_session_key, get_agent_user, quiet_side_effects
from library.vg.repo import REPO_ROOT
from library.vg.settings_chain import resolved_settings_module

# An allow-list, so production - and a new host nobody has classified yet - is refused rather than missed
SHARED_TEST_SETTINGS = {"vgtest2", "sharianttest", "shariantdemo", "shariantsecurity"}
HARNESS_HOME = Path(os.environ.get("VG_BROWSER_HOME", Path.home() / ".cache" / "claude-browser"))
HARNESS_PYTHON = HARNESS_HOME / "pw" / "bin" / "python"
CONSOLE_SCRIPT = REPO_ROOT / "scripts" / "browser" / "console.py"
SESSION_SECONDS = 600
MAX_FRAMES = 6
# Short lines ("});", "}") appear in every template, so only a line this long identifies one
MIN_MATCH_LENGTH = 15
_FRAME_RE = re.compile(r"(?P<url>(?:https?|file)://[^\s()]+?):(?P<line>\d+):(?P<col>\d+)")
_MANIFEST_HASH_RE = re.compile(r"\.[0-9a-f]{12}(\.\w+)$")


class HarnessMissing(Exception):
    pass


class HarnessFailed(Exception):
    pass


class WritesRefused(Exception):
    pass


def check_box_allows_writes(settings_module: str | None = None):
    """ Developer boxes (env_developers) and the shared test deployments only """
    module = settings_module or resolved_settings_module()
    folder, stem = module.split(" ")[0].split(".")[-2:]
    if folder == "env_developers" or stem in SHARED_TEST_SETTINGS:
        return
    raise WritesRefused(f"Refusing --js-writes-to-database on {module}: the page's writes are not rolled back, and this "
                        f"is not a developer box or a test deployment ({', '.join(sorted(SHARED_TEST_SETTINGS))}). "
                        f"Use plain `vg page`, which rolls back.")


@dataclass
class JsReport:
    url: str
    status: int | None
    errors: list[str] = field(default_factory=list)
    console: list[str] = field(default_factory=list)
    requests: list[str] = field(default_factory=list)
    click_failures: list[str] = field(default_factory=list)


class _StaticRootHandler(StaticFilesHandler):
    """ Serves STATIC_URL from STATIC_ROOT: the finders behind StaticFilesHandler don't know hashed names """

    def serve(self, request):
        return static_serve(request, self.file_path(request.path), document_root=settings.STATIC_ROOT)


class _QuietRequestHandler(WSGIRequestHandler):
    def log_message(self, format, *args):  # pylint: disable=redefined-builtin
        pass


@contextmanager
def _serve():
    server = ThreadedWSGIServer(("localhost", 0), _QuietRequestHandler, allow_reuse_address=False)
    server.set_app(_StaticRootHandler(WSGIHandler()))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://localhost:{server.server_address[1]}"
    finally:
        server.shutdown()
        server.server_close()


class TemplateLines:
    """ Finds the template line an inline script line came from, by its text """

    def __init__(self, template_names: list[str]):
        self._index: dict[str, str] = {}
        for name in template_names:
            try:
                origin = get_template(name).origin.name
            except TemplateDoesNotExist:
                continue
            if not origin or not Path(origin).exists():
                continue
            path = Path(origin)
            shown = path.relative_to(REPO_ROOT) if path.is_relative_to(REPO_ROOT) else path
            for line_number, line in enumerate(path.read_text(errors="replace").splitlines(), start=1):
                stripped = line.strip()
                if len(stripped) >= MIN_MATCH_LENGTH:
                    self._index.setdefault(stripped, f"{shown}:{line_number}")

    def find(self, line: str) -> str | None:
        return self._index.get(line.strip())


class FrameMapper:
    """ Rewrites a browser location (url:line:col) as somewhere to look: a template line for the page's own
        inline scripts, a source path for static files """

    def __init__(self, document_url: str, document: str, template_lines: TemplateLines):
        self.document_path = urlsplit(document_url).path
        self.document_lines = document.splitlines()
        self.template_lines = template_lines

    def location(self, url: str, line: int) -> str:
        parts = urlsplit(url)
        if parts.path == self.document_path:
            if 0 < line <= len(self.document_lines):
                source = self.document_lines[line - 1]
                if found := self.template_lines.find(source):
                    return found
                return f"page line {line}: {source.strip()[:120]}"
            return f"page line {line}"
        return f"{_MANIFEST_HASH_RE.sub(r'\1', parts.path)}:{line}"

    def text(self, text: str) -> str:
        """ A message that carries its own stack (jQuery's Deferred warning) """
        return _FRAME_RE.sub(lambda m: self.location(m["url"], int(m["line"])), text)

    def stack(self, stack: str) -> list[str]:
        frames = []
        for frame in stack.splitlines():
            frame = frame.strip()
            if not frame.startswith("at "):
                continue
            if m := _FRAME_RE.search(frame):
                where = self.location(m["url"], int(m["line"]))
                function = frame[3:m.start()].strip(" (")
                frame = f"at {function} ({where})" if function else f"at {where}"
            frames.append(frame)
        if len(frames) > MAX_FRAMES:
            frames = [*frames[:MAX_FRAMES], f"... {len(frames) - MAX_FRAMES} more"]
        return frames


def _counted(lines: list[str]) -> list[str]:
    return [f"{line} (x{n})" if n > 1 else line for line, n in Counter(lines).items()]


def build_report(data: dict, template_names: list[str]) -> JsReport:
    mapper = FrameMapper(data["url"], data["document"], TemplateLines(template_names))
    report = JsReport(url=urlsplit(data["url"]).path, status=data["status"], click_failures=data["click_failures"])
    for error in data["errors"]:
        report.errors.append("\n".join([error["message"], *(f"    {f}" for f in mapper.stack(error["stack"]))]))
    report.console = _counted([
        # console locations are 0-based lines, stack frames 1-based
        f"{m['type']}: {mapper.text(m['text'])}" + (f"  [{mapper.location(m['url'], m['lineNumber'] + 1)}]" if m.get("url") else "")
        for m in data["console"]
    ])
    report.requests = _counted([
        f"{r['method']} {urlsplit(r['url']).path} -> {r.get('status') or r.get('error')}" for r in data["requests"]
    ])
    return report


def run_page_js(url: str, username: str | None = None, template_names: list[str] | None = None,
                wait_ms: int = 1500, clicks: list[str] | None = None) -> JsReport:
    check_box_allows_writes()
    if not HARNESS_PYTHON.exists():
        raise HarnessMissing(f"No browser harness at {HARNESS_HOME}: run scripts/browser/setup.sh once on this box")
    user = get_agent_user(username or AGENT_USERNAME)
    session_key = agent_session_key(user, SESSION_SECONDS)
    command = [str(HARNESS_PYTHON), str(CONSOLE_SCRIPT), "--wait", str(wait_ms),
               *(arg for selector in clicks or [] for arg in ("--click", selector)), url]
    try:
        with quiet_side_effects(block_subprocesses=False), _serve() as base:
            env = {**os.environ, "VG_BROWSER_BASE": base, "VG_BROWSER_SESSION": session_key}
            proc = subprocess.run(command, env=env, capture_output=True, text=True, timeout=300, check=False)
    finally:
        import_module(settings.SESSION_ENGINE).SessionStore(session_key=session_key).delete()
    if proc.returncode:
        raise HarnessFailed(f"{CONSOLE_SCRIPT.name} exited {proc.returncode}:\n{proc.stderr[-2000:]}")
    return build_report(json.loads(proc.stdout), template_names or [])


def render_js_report(report: JsReport) -> str:
    lines = [f"js: {len(report.errors)} uncaught error(s), {len(report.console)} console error/warning(s), "
             f"{len(report.requests)} failed request(s)"]
    lines += [f"  error: {e}" for e in report.errors]
    lines += ["  " + c.replace("\n", "\n    ") for c in report.console]
    lines += [f"  request: {r}" for r in report.requests]
    lines += [f"  click failed (its JS never ran): {c}" for c in report.click_failures]
    return "\n".join(lines)
