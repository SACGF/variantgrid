"""
`vg outline <file.py>`: classes, methods and functions with line numbers and the first line of each
docstring, so a 2,800-line module can be navigated without reading it. `vg outline --coverage` is the
module-docstring ratchet: how many non-test, non-migration modules state their contract, per package.
Pure AST, no Django.
"""
import ast
from dataclasses import dataclass
from pathlib import Path

from library.vg.repo import REPO_ROOT, iter_python_files


@dataclass
class OutlineEntry:
    kind: str          # "class" | "def" | "async def"
    name: str
    lineno: int
    end_lineno: int
    depth: int         # 0 for module level, 1 for methods / nested defs
    doc: str           # first docstring line, or ""
    bases: list[str]   # class bases, as written

    @property
    def length(self) -> int:
        return self.end_lineno - self.lineno + 1


def _source(node) -> str:
    return ast.unparse(node)


def _first_doc_line(node) -> str:
    doc = ast.get_docstring(node, clean=True)
    return doc.strip().splitlines()[0].strip() if doc else ""


def outline(path) -> list[OutlineEntry]:
    tree = ast.parse(Path(path).read_bytes(), filename=str(path))
    entries: list[OutlineEntry] = []

    def visit(body, depth):
        for node in body:
            if isinstance(node, ast.ClassDef):
                entries.append(OutlineEntry("class", node.name, node.lineno, node.end_lineno, depth,
                                            _first_doc_line(node), [_source(b) for b in node.bases]))
                visit(node.body, depth + 1)
            elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                kind = "async def" if isinstance(node, ast.AsyncFunctionDef) else "def"
                entries.append(OutlineEntry(kind, node.name, node.lineno, node.end_lineno, depth,
                                            _first_doc_line(node), []))
                if depth == 0:
                    visit(node.body, depth + 1)  # nested defs one level down are worth seeing; deeper are not

    visit(tree.body, 0)
    return entries


def module_doc(path) -> str:
    return _first_doc_line(ast.parse(Path(path).read_bytes()))


def render_outline(path, entries: list[OutlineEntry], min_lines: int = 0) -> str:
    """ One line per entry: `L123-456  class Variant(Model)  - first doc line`; the header names the module """
    lines = [f"{path}  ({_line_count(path)} lines)"]
    doc = module_doc(path)
    lines.append(f"  {doc}" if doc else "  (no module docstring)")
    for entry in entries:
        if entry.depth and entry.length < min_lines:
            continue
        indent = "    " * entry.depth
        signature = entry.name + (f"({', '.join(entry.bases)})" if entry.bases else "")
        doc = f"  - {entry.doc}" if entry.doc else ""
        lines.append(f"{entry.lineno:>5}-{entry.end_lineno:<5} {indent}{entry.kind} {signature}{doc}")
    return "\n".join(lines)


def _line_count(path) -> int:
    with open(path, "rb") as f:
        return sum(1 for _ in f)


def has_module_doc(path) -> bool:
    try:
        return ast.get_docstring(ast.parse(Path(path).read_bytes())) is not None
    except (SyntaxError, ValueError):
        return False


def _counts_toward_coverage(path: Path) -> bool:
    """ Tests and migrations are out; an __init__ counts only when it holds code (a re-exporting package is a contract too) """
    parts = path.relative_to(REPO_ROOT).parts
    if "tests" in parts or "migrations" in parts:
        return False
    return path.name != "__init__.py" or bool(path.read_text().strip())


def docstring_coverage(packages=None) -> dict[str, tuple[int, int]]:
    """ package -> (modules with a docstring, modules counted); tests, migrations and empty __init__ files are excluded """
    coverage: dict[str, list[int]] = {}
    for module, path in iter_python_files(packages):
        if not module or not _counts_toward_coverage(path):
            continue
        package = module.split(".")[0]
        counts = coverage.setdefault(package, [0, 0])
        counts[0] += has_module_doc(path)
        counts[1] += 1
    return {package: (documented, total) for package, (documented, total) in sorted(coverage.items())}


def render_coverage(packages=None) -> str:
    coverage = docstring_coverage(packages)
    documented = sum(d for d, _ in coverage.values())
    total = sum(t for _, t in coverage.values())
    lines = [f"module docstring coverage: {documented}/{total} ({documented / total:.0%})" if total else "no modules"]
    lines += [f"  {package:<16} {d:>4}/{t:<4} {d / t:.0%}" for package, (d, t) in coverage.items()]
    return "\n".join(lines)
