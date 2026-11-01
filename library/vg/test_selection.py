"""
`vg tests --changed`: which test modules does a set of changed files put at risk?

Rules:
- a changed test module selects itself
- a changed Python module selects every test module that transitively imports it
- a changed web module (anything under a views/, grids/, forms/ or templatetags/ package, or a module whose name
  has a view / grid / form word, e.g. views_rest.py, grid_export.py, node_view.py) also selects the owning app's
  tests.test_urls: tests reach views through URLs, not imports. One that still selects nothing (its app has no
  test_urls) is reported in `warnings`
- a changed template tag library also selects tests.test_urls of every app whose templates `{% load %}` it
- a changed template, static JS/CSS or urls.py selects the owning app's tests.test_urls
- a changed migration selects the owning app's whole tests package
- a changed settings module selects nothing, with a warning: no import edge says which code reads a setting
- everything else (docs, data, scripts) selects nothing

Entry points: `select_tests(changed=None, base=None)`, and `paths_as_changed` for paths named on a command line.
"""
import re
from dataclasses import dataclass, field
from functools import cache
from pathlib import Path

from library.vg.import_graph import ImportGraph, build_import_graph
from library.vg.repo import (
    REPO_ROOT,
    app_for_path,
    changed_files,
    first_party_packages,
    module_name_for,
)

_STATIC_APP_RE = re.compile(r"static_files/[^/]+/(?:js|css|scss)/([^/]+)/")
_LOAD_TAG_RE = re.compile(r"{%\s*load\s+([^%]+?)\s*%}")
_WEB_PACKAGES = {"views", "grids", "forms", "templatetags"}
_WEB_WORDS = {"view", "views", "grid", "grids", "form", "forms"}
_SETTINGS_DIR = Path("variantgrid", "settings")


@dataclass
class TestSelection:
    changed: list[str]
    labels: list[str]                                    # manage.py test labels, sorted
    reasons: dict[str, list[str]] = field(default_factory=dict)  # label -> changed paths that selected it
    warnings: list[str] = field(default_factory=list)    # changed paths whose risk the rules cannot see


def paths_as_changed(paths: list[str]) -> list[str]:
    """ Command-line paths (relative to the shell's directory, or to the repo root) as repo-relative paths """
    changed = []
    for path in paths:
        resolved = Path(path).resolve()
        if resolved.exists() and resolved.is_relative_to(REPO_ROOT):
            changed.append(str(resolved.relative_to(REPO_ROOT)))
        else:
            changed.append(path)
    return changed


def _is_test_module(module: str) -> bool:
    parts = module.split(".")
    return "tests" in parts and (parts[-1].startswith("test") or parts[-1] == "tests")


def _is_web_module(path: Path) -> bool:
    return bool(_WEB_PACKAGES.intersection(path.parts[:-1])) or bool(_WEB_WORDS.intersection(path.stem.split("_")))


def _test_urls_label(app: str) -> str | None:
    return f"{app}.tests.test_urls" if (REPO_ROOT / app / "tests" / "test_urls.py").exists() else None


def _app_tests_label(app: str) -> str | None:
    return f"{app}.tests" if (REPO_ROOT / app / "tests").is_dir() else None


@cache
def _template_tag_loaders() -> dict[str, set[str]]:
    """ Template tag library name -> apps with a template that loads it """
    loaders: dict[str, set[str]] = {}
    for package in first_party_packages():
        for template in (REPO_ROOT / package).rglob("templates/**/*.html"):
            for match in _LOAD_TAG_RE.finditer(template.read_text(errors="replace")):
                words = match.group(1).split()
                libraries = words[words.index("from") + 1:] if "from" in words else words
                for library in libraries:
                    loaders.setdefault(library, set()).add(package)
    return loaders


def _select(selection: dict[str, set[str]], label: str | None, reason: str) -> bool:
    if label:
        selection.setdefault(label, set()).add(reason)
    return label is not None


def select_tests(changed: list[str] | None = None, base: str | None = None,
                 graph: ImportGraph | None = None) -> TestSelection:
    changed = changed if changed is not None else changed_files(base)
    graph = graph or build_import_graph()
    selection: dict[str, set[str]] = {}
    warnings: list[str] = []

    for path in changed:
        app = app_for_path(path)
        if not app:
            continue
        parts = Path(path).parts

        if Path(path).is_relative_to(_SETTINGS_DIR):
            warnings.append(f"{path}: settings changes select no tests; run the suite or the tests of the app "
                            f"the setting belongs to")
            continue

        if "migrations" in parts:
            _select(selection, _app_tests_label(app), path)
            continue

        if path.endswith(".py"):
            module = module_name_for(REPO_ROOT / path)
            if not module:
                continue
            if _is_test_module(module):
                _select(selection, module, path)
                continue
            selected = False
            for dependent in graph.dependents(module):
                if _is_test_module(dependent):
                    selected |= _select(selection, dependent, path)
            if Path(path).name == "urls.py" or _is_web_module(Path(path)):
                selected |= _select(selection, _test_urls_label(app), path)
            if "templatetags" in parts:
                for loader in _template_tag_loaders().get(Path(path).stem, ()):
                    selected |= _select(selection, _test_urls_label(loader), path)
            if not selected and _is_web_module(Path(path)):
                warnings.append(f"{path}: {app} has no tests/test_urls.py and no test imports this, so it selects no tests")
            continue

        if "templates" in parts or path.endswith((".html", ".js", ".css", ".scss")):
            static_match = _STATIC_APP_RE.search(path)
            owner = static_match.group(1) if static_match and app_for_path(static_match.group(1)) else app
            _select(selection, _test_urls_label(owner), path)

    # A test package label (from a migration) makes its individual module labels redundant
    packages = {label for label in selection if label.endswith(".tests")}
    for label in list(selection):
        package = next((p for p in packages if label != p and label.startswith(p + ".")), None)
        if package:
            selection[package].update(selection.pop(label))

    labels = sorted(selection)
    return TestSelection(changed=changed, labels=labels,
                         reasons={label: sorted(selection[label]) for label in labels}, warnings=warnings)
