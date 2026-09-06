"""
`vg outline <file.py>`: classes, methods and functions with line numbers and the first line of each
docstring, so a 2,800-line module can be navigated without reading it. Pure AST, no Django.
"""
import ast
from dataclasses import dataclass
from pathlib import Path


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
