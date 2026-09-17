"""
First-party import graph built from the AST, without importing anything. Used by test selection
(which test modules transitively import a changed module) and by `vg map` to rank fan-in.

`from a.b import c` resolves to the module `a.b.c` when that module exists, otherwise to `a.b` (c is a
symbol). Relative imports are resolved against the importing package. Third-party imports are dropped.

Parsing ~2,500 modules takes a few seconds, so each file's raw import statements are cached in
.vg_cache/ keyed on (mtime, size); a warm run resolves the graph in well under a second.
"""
import ast
import json
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from library.vg.repo import REPO_ROOT, first_party_packages, iter_python_files

CACHE_FILE = REPO_ROOT / ".vg_cache" / "imports.json"

_BODY_ATTRS = ("body", "orelse", "finalbody", "handlers")


@dataclass
class ImportGraph:
    modules: dict[str, Path] = field(default_factory=dict)
    imports: dict[str, set[str]] = field(default_factory=lambda: defaultdict(set))  # module -> modules it imports
    imported_by: dict[str, set[str]] = field(default_factory=lambda: defaultdict(set))

    def add_edge(self, importer: str, imported: str):
        self.imports[importer].add(imported)
        self.imported_by[imported].add(importer)

    def dependents(self, module: str) -> set[str]:
        """ Every module that transitively imports `module` (excluding itself) """
        seen: set[str] = set()
        frontier = [module]
        while frontier:
            current = frontier.pop()
            for importer in self.imported_by.get(current, ()):
                if importer not in seen:
                    seen.add(importer)
                    frontier.append(importer)
        seen.discard(module)
        return seen

    def fan_in(self) -> dict[str, int]:
        return {module: len(importers) for module, importers in self.imported_by.items()}


def iter_import_statements(node):
    """ Import statements anywhere in the file, walking only statement bodies (not expressions) """
    for attr in _BODY_ATTRS:
        for child in getattr(node, attr, ()) or ():
            if isinstance(child, (ast.Import, ast.ImportFrom)):
                yield child
            elif isinstance(child, (ast.stmt, ast.ExceptHandler)):
                yield from iter_import_statements(child)


def _raw_imports(path: Path) -> list[list]:
    """ [["import", "a.b"], ["from", level, base_or_None, [names]]] for one file """
    try:
        tree = ast.parse(path.read_bytes(), filename=str(path))
    except (SyntaxError, ValueError):
        return []
    raw = []
    for node in iter_import_statements(tree):
        if isinstance(node, ast.Import):
            raw.extend(["import", alias.name] for alias in node.names)
        else:
            raw.append(["from", node.level, node.module, [alias.name for alias in node.names]])
    return raw


def _load_cache() -> dict:
    try:
        return json.loads(CACHE_FILE.read_text())
    except (OSError, ValueError):
        return {}


def _save_cache(cache: dict):
    try:
        CACHE_FILE.parent.mkdir(exist_ok=True)
        CACHE_FILE.write_text(json.dumps(cache))
    except OSError:
        pass


def resolve_from_import(known: set[str], base: str, name: str) -> str | None:
    """ from <base> import <name>: the submodule if it is one, else the package/module itself """
    candidate = f"{base}.{name}" if base else name
    if candidate in known:
        return candidate
    if base in known:
        return base
    return None


def relative_base(importer: str, is_package: bool, level: int, module: str | None) -> str:
    parts = importer.split(".")
    if not is_package:
        parts = parts[:-1]
    drop = level - 1
    if drop:
        parts = parts[:-drop] if drop <= len(parts) else []
    if module:
        parts.append(module)
    return ".".join(parts)


def build_import_graph() -> ImportGraph:
    graph = ImportGraph()
    for module, path in iter_python_files():
        if module:
            graph.modules[module] = path
    known = set(graph.modules)
    first_party = set(first_party_packages())

    cache = _load_cache()
    fresh_cache = {}
    for module, path in graph.modules.items():
        stat = path.stat()
        key = str(path.relative_to(REPO_ROOT))
        stamp = [stat.st_mtime_ns, stat.st_size]
        entry = cache.get(key)
        if not entry or entry["stamp"] != stamp:
            entry = {"stamp": stamp, "imports": _raw_imports(path)}
        fresh_cache[key] = entry

        is_package = path.name == "__init__.py"
        for raw in entry["imports"]:
            if raw[0] == "import":
                target = raw[1]
                if target.split(".")[0] not in first_party:
                    continue
                while target and target not in known:
                    target = target.rpartition(".")[0]
                if target:
                    graph.add_edge(module, target)
            else:
                _, level, base_module, names = raw
                base = relative_base(module, is_package, level, base_module) if level else (base_module or "")
                if base.split(".")[0] not in first_party:
                    continue
                for name in names:
                    target = resolve_from_import(known, base, name)
                    if target:
                        graph.add_edge(module, target)

    if fresh_cache != cache:
        _save_cache(fresh_cache)
    return graph
