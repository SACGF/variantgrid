"""
`vg tests --changed`: which test modules does a set of changed files put at risk?

Rules (see claude/plans/1816_agent_system_phase1_plan.md):
- a changed test module selects itself
- a changed Python module selects every test module that transitively imports it
- a changed template, static JS/CSS or urls.py selects the owning app's tests.test_urls
- a changed migration selects the owning app's whole tests package
- everything else (docs, data, scripts) selects nothing
"""
import re
from dataclasses import dataclass, field
from pathlib import Path

from library.vg.import_graph import ImportGraph, build_import_graph
from library.vg.repo import REPO_ROOT, app_for_path, changed_files, module_name_for

_STATIC_APP_RE = re.compile(r"static_files/[^/]+/(?:js|css|scss)/([^/]+)/")


@dataclass
class TestSelection:
    changed: list[str]
    labels: list[str]                                    # manage.py test labels, sorted
    reasons: dict[str, list[str]] = field(default_factory=dict)  # label -> changed paths that selected it


def _is_test_module(module: str) -> bool:
    parts = module.split(".")
    return "tests" in parts and (parts[-1].startswith("test") or parts[-1] == "tests")


def _test_urls_label(app: str) -> str | None:
    return f"{app}.tests.test_urls" if (REPO_ROOT / app / "tests" / "test_urls.py").exists() else None


def _app_tests_label(app: str) -> str | None:
    return f"{app}.tests" if (REPO_ROOT / app / "tests").is_dir() else None


def _select(selection: dict[str, set[str]], label: str | None, reason: str):
    if label:
        selection.setdefault(label, set()).add(reason)


def select_tests(changed: list[str] | None = None, base: str | None = None,
                 graph: ImportGraph | None = None) -> TestSelection:
    changed = changed if changed is not None else changed_files(base)
    graph = graph or build_import_graph()
    selection: dict[str, set[str]] = {}

    for path in changed:
        app = app_for_path(path)
        if not app:
            continue
        parts = Path(path).parts

        if "migrations" in parts:
            _select(selection, _app_tests_label(app), path)
            continue

        if path.endswith(".py"):
            module = module_name_for(REPO_ROOT / path)
            if not module:
                continue
            if _is_test_module(module):
                _select(selection, module, path)
            for dependent in graph.dependents(module):
                if _is_test_module(dependent):
                    _select(selection, dependent, path)
            if Path(path).name == "urls.py":
                _select(selection, _test_urls_label(app), path)
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
                         reasons={label: sorted(selection[label]) for label in labels})
