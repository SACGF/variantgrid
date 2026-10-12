"""
`vg deprecations`: first-party code that uses something deprecated, found without running the test suite.

Static: basedpyright's `reportDeprecated` over every first-party package reads the PEP 702 `@deprecated` markers in
typeshed, the installed stubs (django-stubs) and our own code. It runs as `LOOKAHEAD_PYTHON`, so a call deprecated in a
Python we don't deploy yet already shows. A hit is "ours" when the deprecated name is one a first-party `@deprecated`
declares (a migration still in progress); otherwise it is the standard library or a package, to fix before an upgrade
removes it.

Runtime: `manage.py check` under `-W always::DeprecationWarning` (and pending / future) reports the `warnings.warn`
calls made while every app, model and URLconf is imported, which is how Django and most packages deprecate. Warnings
raised only when a function is called show up in test output instead.

Entry points are `find_deprecations` and `render_deprecations`. No Django here: both halves are subprocesses.
"""
import ast
import json
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from library.vg.repo import REPO_ROOT, first_party_packages, iter_python_files

# The newest Python typeshed has stubs for: its version-gated @deprecated markers are what we'll hit on upgrade
LOOKAHEAD_PYTHON = "3.15"

_PYRIGHT_SYMBOL_RE = re.compile(r'"(?P<name>[^"]+)"(?: in class "(?P<cls>[^"]+)")?')
_RUNTIME_WARNING_RE = re.compile(r"^(?P<path>\S+\.py):(?P<line>\d+): (?P<category>\w+): (?P<message>.*)$")
_DEPRECATION_CATEGORY_RE = re.compile(r"Deprecation|RemovedIn|FutureWarning")
_RUNTIME_WARNING_FLAGS = ["-W", "always::DeprecationWarning", "-W", "always::PendingDeprecationWarning",
                          "-W", "always::FutureWarning"]


@dataclass
class StaticDeprecation:
    path: str
    line: int
    symbol: str
    reason: str
    ours: bool


@dataclass
class RuntimeDeprecation:
    path: str
    line: int
    category: str
    message: str
    ours: bool


@dataclass
class DeprecationReport:
    static: list[StaticDeprecation] = field(default_factory=list)
    runtime: list[RuntimeDeprecation] = field(default_factory=list)
    unavailable: dict[str, str] = field(default_factory=dict)

    @property
    def external(self) -> list[StaticDeprecation]:
        return [d for d in self.static if not d.ours]

    @property
    def ours(self) -> list[StaticDeprecation]:
        return [d for d in self.static if d.ours]

    @property
    def runtime_in_repo(self) -> list[RuntimeDeprecation]:
        return [d for d in self.runtime if d.ours]

    @property
    def ok(self) -> bool:
        return not self.external and not self.runtime_in_repo


def _is_deprecated_decorator(node: ast.expr) -> bool:
    target = node.func if isinstance(node, ast.Call) else node
    return (isinstance(target, ast.Name) and target.id == "deprecated") or \
        (isinstance(target, ast.Attribute) and target.attr == "deprecated")


def deprecated_names_in_source(source: str) -> set[str]:
    """ Names a module marks `@deprecated`: 'func', 'Class', or 'Class.method' """
    names = set()

    def visit(body, prefix: str):
        for node in body:
            if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
                if any(map(_is_deprecated_decorator, node.decorator_list)):
                    names.add(prefix + node.name)
                if isinstance(node, ast.ClassDef):
                    visit(node.body, f"{prefix}{node.name}.")

    visit(ast.parse(source).body, "")
    return names


def our_deprecated_names() -> set[str]:
    names = set()
    for _, path in iter_python_files():
        source = path.read_text(errors="replace")
        if "deprecated" in source:
            try:
                names |= deprecated_names_in_source(source)
            except SyntaxError:
                continue
    return names


def pyright_symbol(message: str) -> str:
    """ 'The method "utcnow" in class "datetime" is deprecated' -> 'datetime.utcnow' """
    if m := _PYRIGHT_SYMBOL_RE.search(message):
        return f"{m['cls']}.{m['name']}" if m["cls"] else m["name"]
    return message.splitlines()[0]


def _pyright_config() -> dict:
    return {
        "include": [str(REPO_ROOT / package) for package in first_party_packages()],
        "exclude": ["**/node_modules", "**/__pycache__", "**/.*"],
        "executionEnvironments": [{"root": str(REPO_ROOT)}],
        "venvPath": str(REPO_ROOT),
        "venv": ".venv",
        "pythonVersion": LOOKAHEAD_PYTHON,
        "typeCheckingMode": "off",
        "reportDeprecated": "warning",
    }


def _basedpyright() -> str | None:
    candidate = REPO_ROOT / ".venv" / "bin" / "basedpyright"
    return str(candidate) if candidate.exists() else shutil.which("basedpyright")


def static_deprecations(pyright_json: dict, ours: set[str]) -> list[StaticDeprecation]:
    found = []
    for diagnostic in pyright_json.get("generalDiagnostics", []):
        if diagnostic.get("rule") != "reportDeprecated":
            continue
        headline, _, reason = diagnostic["message"].partition("\n")
        symbol = pyright_symbol(headline)
        found.append(StaticDeprecation(path=str(Path(diagnostic["file"]).relative_to(REPO_ROOT)),
                                       line=diagnostic["range"]["start"]["line"] + 1,
                                       symbol=symbol, reason=reason.strip(), ours=symbol in ours))
    return found


def runtime_deprecations(stderr: str) -> list[RuntimeDeprecation]:
    """ Deprecation warnings in Python's 'path:line: Category: message' form, de-duplicated """
    found = {}
    for line in stderr.splitlines():
        if (m := _RUNTIME_WARNING_RE.match(line)) and _DEPRECATION_CATEGORY_RE.search(m["category"]):
            path = Path(m["path"])
            in_repo = path.is_relative_to(REPO_ROOT) and not path.is_relative_to(REPO_ROOT / ".venv")
            if in_repo:
                shown = str(path.relative_to(REPO_ROOT))
            else:
                shown = m["path"].split("site-packages/")[-1]
            key = (shown, int(m["line"]), m["category"], m["message"])
            found[key] = RuntimeDeprecation(*key, ours=in_repo)
    return list(found.values())


def find_deprecations(runtime: bool = True) -> DeprecationReport:
    report = DeprecationReport()
    check = None
    if runtime:
        python = REPO_ROOT / ".venv" / "bin" / "python"
        check = subprocess.Popen([str(python), *_RUNTIME_WARNING_FLAGS, str(REPO_ROOT / "manage.py"), "check"],
                                 cwd=REPO_ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)

    if pyright := _basedpyright():
        with tempfile.TemporaryDirectory() as config_dir:
            config = Path(config_dir) / "pyrightconfig.json"
            config.write_text(json.dumps(_pyright_config()))
            result = subprocess.run([pyright, "-p", str(config), "--outputjson"], cwd=REPO_ROOT,
                                    capture_output=True, text=True, timeout=900)
        try:
            report.static = static_deprecations(json.loads(result.stdout), our_deprecated_names())
        except json.JSONDecodeError:
            report.unavailable["static"] = f"basedpyright exited {result.returncode}: {result.stderr.strip()[:200]}"
    else:
        report.unavailable["static"] = "basedpyright not installed - scripts/install_requirements.sh --dev"

    if check:
        _, stderr = check.communicate(timeout=600)
        report.runtime = runtime_deprecations(stderr)
        if check.returncode and not report.runtime:
            report.unavailable["runtime"] = f"manage.py check exited {check.returncode}: {stderr.strip()[-200:]}"
    return report


def _by_symbol(hits: list[StaticDeprecation]) -> dict[str, list[StaticDeprecation]]:
    grouped = defaultdict(list)
    for hit in sorted(hits, key=lambda h: (h.symbol, h.path, h.line)):
        grouped[hit.symbol].append(hit)
    return grouped


def render_deprecations(report: DeprecationReport) -> str:
    lines = [f"Standard library / packages, as Python {LOOKAHEAD_PYTHON} (fix these): {len(report.external)}"]
    for symbol, hits in _by_symbol(report.external).items():
        lines.append(f"  {symbol} - {hits[0].reason}")
        lines.extend(f"    {h.path}:{h.line}" for h in hits)
    lines.append(f"Our own @deprecated, still in use: {len(report.ours)}")
    for symbol, hits in _by_symbol(report.ours).items():
        files = len({h.path for h in hits})
        lines.append(f"  {symbol} - {len(hits)} use(s) in {files} file(s) - {hits[0].reason}")
    lines.append(f"Warned at import (manage.py check): {len(report.runtime)}")
    for hit in sorted(report.runtime, key=lambda h: (not h.ours, h.path, h.line)):
        where = "" if hit.ours else " [installed package]"
        lines.append(f"  {hit.path}:{hit.line}{where} {hit.category}: {hit.message}")
    lines.extend(f"{name} unavailable: {why}" for name, why in report.unavailable.items())
    return "\n".join(lines)
