"""
First-party import graph built from the AST, without importing anything. Used by test selection
(which test modules transitively import a changed module), by `vg map` to rank fan-in, and by
`vg imports cycles`, which fails on any cycle among the imports that run at module load.

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


def iter_module_load_imports(node):
    """ Import statements that run when the module is imported: not those inside a function body
        (the sanctioned way to break a cycle) or under `if TYPE_CHECKING:` """
    for attr in _BODY_ATTRS:
        for child in getattr(node, attr, ()) or ():
            if isinstance(child, (ast.Import, ast.ImportFrom)):
                yield child
            elif isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            elif isinstance(child, ast.If) and "TYPE_CHECKING" in ast.unparse(child.test):
                for orelse_child in child.orelse:
                    yield from iter_module_load_imports(ast.Module(body=[orelse_child], type_ignores=[]))
            elif isinstance(child, (ast.stmt, ast.ExceptHandler)):
                yield from iter_module_load_imports(child)


def _raw_imports(path: Path, module_load_only: bool = False) -> list[list]:
    """ [["import", "a.b"], ["from", level, base_or_None, [names]]] for one file """
    try:
        tree = ast.parse(path.read_bytes(), filename=str(path))
    except (SyntaxError, ValueError):
        return []
    raw = []
    walker = iter_module_load_imports if module_load_only else iter_import_statements
    for node in walker(tree):
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


def build_import_graph(module_load_only: bool = False) -> ImportGraph:
    """ module_load_only: just the imports that run at import time (uncached - the cache holds every import) """
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
        if module_load_only:
            entry = {"imports": _raw_imports(path, module_load_only=True)}
        else:
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

    if not module_load_only and fresh_cache != cache:
        _save_cache(fresh_cache)
    return graph


def _is_production_module(module: str) -> bool:
    """ Tests are leaves nothing imports and migrations are frozen, so neither can be part of a cycle worth fixing """
    parts = module.split(".")
    return not any(part in ("tests", "migrations") or part.startswith("test_") for part in parts)


def import_cycles(graph: ImportGraph) -> list[list[str]]:
    """ Strongly connected components of more than one module (Tarjan, iterative), largest first """
    edges = {module: sorted(t for t in targets if _is_production_module(t))
             for module, targets in graph.imports.items() if _is_production_module(module)}
    index: dict[str, int] = {}
    lowlink: dict[str, int] = {}
    on_stack: set[str] = set()
    stack: list[str] = []
    cycles = []
    for root in sorted(edges):
        if root in index:
            continue
        work = [(root, iter(edges.get(root, ())))]
        index[root] = lowlink[root] = len(index)
        stack.append(root)
        on_stack.add(root)
        while work:
            node, children = work[-1]
            child = next(children, None)
            if child is not None:
                if child not in index:
                    index[child] = lowlink[child] = len(index)
                    stack.append(child)
                    on_stack.add(child)
                    work.append((child, iter(edges.get(child, ()))))
                elif child in on_stack:
                    lowlink[node] = min(lowlink[node], index[child])
                continue
            work.pop()
            if work:
                parent = work[-1][0]
                lowlink[parent] = min(lowlink[parent], lowlink[node])
            if lowlink[node] == index[node]:
                component = []
                while True:
                    member = stack.pop()
                    on_stack.discard(member)
                    component.append(member)
                    if member == node:
                        break
                if len(component) > 1:
                    cycles.append(sorted(component))
    return sorted(cycles, key=lambda c: (-len(c), c))


def render_cycles(graph: ImportGraph, cycles: list[list[str]]) -> str:
    if not cycles:
        return "No import cycles among module-load imports."
    lines = [f"{len(cycles)} import cycle(s) among module-load imports. Import each name from the module that "
             "defines it rather than through a package __init__, or move the import into the function that needs it:"]
    for cycle in cycles:
        members = set(cycle)
        lines.append("")
        for module in cycle:
            for target in sorted(graph.imports.get(module, ()) & members):
                lines.append(f"  {module} -> {target}")
    return "\n".join(lines)
