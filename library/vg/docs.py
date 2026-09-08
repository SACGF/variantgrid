"""
`vg docs check`: every path and `path:Symbol` citation in the agent docs, resolved against the working
tree, so a doc with dead citations announces itself instead of misleading the next reader.

Owns: which files count as docs (DOC_GLOBS), what counts as a citation (a backticked path, a bare
`path.py:Symbol`, a Markdown link to a repo file) and how a citation resolves (file exists; the dotted
symbol is a top-level class / function / assignment or a member of one, by AST). Entry points are
`check_docs` (returns a list of DeadCitation) and `render_report`. Pure AST and filesystem - no Django,
so scripts/vg runs it in well under a second and CI runs it without a database.

Paths resolve against the repo root first, then the citing doc's directory (so `snpdb/CLAUDE.md` may say
`models/models_variant.py:Variant`); a bare filename with no directory resolves to any file of that name under
the doc's directory tree, then anywhere in the repo. A `module.path:Symbol` spelling is accepted as
`module/path.py:Symbol`. Fenced code blocks are skipped: what is inside them is a command or a listing, not a
claim about the tree. `claude/maps/` is generated and verified by `vg map --check`, so it is not a doc here.
A citation of gitignored build output (BUILD_OUTPUTS - collectstatic's `variantgrid/sitestatic/`, the `lint.txt`
report) resolves without touching the filesystem, so a doc describing the build reads the same on CI as on a box
that has run it.

Plans are checked while they are live: a `Status:` of draft, approved or in progress. A landed or superseded
plan, or one with no Status line, describes code as it was and is reported as unchecked rather than failed.
"""
import ast
import re
from dataclasses import dataclass
from functools import cache
from pathlib import Path

from library.vg.repo import REPO_ROOT, first_party_packages

DOC_GLOBS = ("claude/**/*.md", "*/CLAUDE.md", "*/__*_readme.md", "CLAUDE.md")
GENERATED_DIRS = ("claude/maps",)
# Gitignored build output: present on a box that has run the build, absent from a fresh checkout
BUILD_OUTPUTS = ("variantgrid/sitestatic", "lint.txt")
LIVE_PLAN_STATUSES = ("draft", "approved", "in progress")
_STATUS_RE = re.compile(r"^Status:\s*(.+)$", re.MULTILINE)

# Extensions that make a bare token a file citation even without a directory separator
_FILE_SUFFIXES = (".py", ".md", ".html", ".js", ".scss", ".css", ".sh", ".yml", ".yaml", ".json", ".txt",
                  ".ini", ".toml", ".cfg", ".csv", ".tsv", ".sql", ".vcf", ".gz", ".png")
# Directories a slash-separated token may start with and still be a repo path (packages are added at runtime)
_REPO_DIRS = ("claude", "scripts", "config", ".claude", ".github", "variantgrid", "docs", "data", "test_data")

_FENCE_RE = re.compile(r"^(```|~~~)")
_INLINE_CODE_RE = re.compile(r"`([^`\n]+)`")
_MD_LINK_RE = re.compile(r"\[[^\]]*\]\(([^)\s]+)\)")
_DOTTED_NAME = r"[A-Za-z_]\w*(?:\.[A-Za-z_]\w*)*"
_BARE_CITATION_RE = re.compile(r"(?<![\w/`.\-])((?:[\w.\-]+/)*[\w\-]+\.py):(" + _DOTTED_NAME + ")")
_PATH_TOKEN_RE = re.compile(r"^[\w./\-]+(?:#[\w\-]+)?$")
_SYMBOL_TOKEN_RE = re.compile(r"^([\w./\-]+):(" + _DOTTED_NAME + r")\.?$")
_MODULE_SPELLING_RE = re.compile(r"^[A-Za-z_][\w]*(?:\.[A-Za-z_]\w*)+$")


@dataclass
class Citation:
    doc: Path
    lineno: int
    text: str            # the citation as written
    path: str            # the path part
    symbol: str = ""     # dotted symbol after the colon, if any
    anchor: str = ""     # #fragment on a markdown path, if any


@dataclass
class DeadCitation:
    citation: Citation
    reason: str

    @property
    def location(self) -> str:
        return f"{self.citation.doc.relative_to(REPO_ROOT)}:{self.citation.lineno}"


def doc_files(paths=None) -> list[Path]:
    """ The docs covered by the check (DOC_GLOBS), or the given files / directories narrowed to .md """
    if paths:
        files = []
        for given in paths:
            path = Path(given)
            if not path.is_absolute():
                path = REPO_ROOT / path
            if path.is_dir():
                files.extend(sorted(path.rglob("*.md")))
            else:
                files.append(path)
        return files
    files = set()
    for pattern in DOC_GLOBS:
        files.update(p for p in REPO_ROOT.glob(pattern) if ".venv" not in p.parts and "node_modules" not in p.parts)
    return sorted(p for p in files if not is_generated(p))


def is_generated(doc: Path) -> bool:
    relative = doc.relative_to(REPO_ROOT).as_posix()
    return any(relative.startswith(prefix + "/") for prefix in GENERATED_DIRS)


def is_build_output(path: str) -> bool:
    """ A citation of collectstatic output or the lint report: the doc is describing the build, and the file
        is gitignored, so it resolves the same on CI as on a box that has run it """
    return any(path == output or path.startswith(output + "/") for output in BUILD_OUTPUTS)


def is_plan(doc: Path) -> bool:
    relative = doc.relative_to(REPO_ROOT).as_posix()
    return relative.startswith("claude/") and (relative.startswith("claude/plans/") or relative.endswith("_plan.md"))


def plan_status(doc: Path) -> str | None:
    """ The plan's Status line (lowercased, first word or two), or None when it has none """
    match = _STATUS_RE.search(doc.read_text(errors="replace")[:4000])
    return match.group(1).strip().lower() if match else None


def is_live_plan(doc: Path) -> bool:
    status = plan_status(doc)
    return status is not None and status.startswith(LIVE_PLAN_STATUSES)


def _looks_like_path(token: str) -> bool:
    """ A backticked token is a file citation when it is a path (has a source suffix, or a slash from a repo dir) """
    if not _PATH_TOKEN_RE.match(token) or token.startswith(("/", "./", "../", "~")) or ".." in token.split("/"):
        return False
    base = token.split("#", 1)[0]
    if base.endswith("/"):
        base = base[:-1]
    if "*" in base or not base or (base.startswith(".") and "/" not in base):
        return False
    if any(base.endswith(suffix) for suffix in _FILE_SUFFIXES):
        return True
    if "/" in base:
        head = base.split("/", 1)[0]
        return head in _REPO_DIRS or head in first_party_packages()
    return False


def _module_to_path(spelling: str) -> str:
    return spelling.replace(".", "/") + ".py"


def _split_symbol_token(token: str) -> tuple[str, str] | None:
    """ 'snpdb/models/x.py:Variant.save' -> (path, symbol); 'snpdb.models.x:Variant' -> (path, symbol) """
    match = _SYMBOL_TOKEN_RE.match(token)
    if not match:
        return None
    path, symbol = match.groups()
    if _MODULE_SPELLING_RE.match(path) and path.split(".")[0] in first_party_packages():
        path = _module_to_path(path)
    if not _looks_like_path(path):
        return None
    return path, symbol


def _iter_prose_lines(text: str):
    """ (lineno, line) for lines outside fenced code blocks """
    in_fence = False
    for lineno, line in enumerate(text.splitlines(), start=1):
        if _FENCE_RE.match(line.strip()):
            in_fence = not in_fence
            continue
        if not in_fence:
            yield lineno, line


def citations_in(doc: Path) -> list[Citation]:
    citations: list[Citation] = []
    seen: set[tuple[int, str]] = set()

    def add(lineno, text, path, symbol="", anchor=""):
        if (lineno, text) not in seen:
            seen.add((lineno, text))
            citations.append(Citation(doc, lineno, text, path, symbol, anchor))

    for lineno, line in _iter_prose_lines(doc.read_text(errors="replace")):
        for match in _INLINE_CODE_RE.finditer(line):
            token = match.group(1).strip()
            split = _split_symbol_token(token)
            if split:
                add(lineno, token, *split)
            elif _looks_like_path(token):
                path, _, anchor = token.partition("#")
                add(lineno, token, path.rstrip("/"), anchor=anchor)
        for match in _MD_LINK_RE.finditer(_INLINE_CODE_RE.sub("", line)):
            target = match.group(1)
            if "://" in target or target.startswith(("#", "mailto:", "/")):
                continue
            path, _, anchor = target.partition("#")
            add(lineno, target, path, anchor=anchor)
        for match in _BARE_CITATION_RE.finditer(_INLINE_CODE_RE.sub("", line)):
            path, symbol = match.groups()
            if _looks_like_path(path):
                add(lineno, match.group(0), path, symbol)
    return citations


def resolve_path(citation: Citation) -> Path | None:
    candidates = [REPO_ROOT / citation.path, citation.doc.parent / citation.path]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    if "/" not in citation.path:
        return _resolve_bare_filename(citation)
    return None


_SKIP_DIRS = {".venv", "node_modules", "__pycache__", ".git", "sitestatic", ".vg_cache", "media_root"}


@cache
def _files_by_name() -> dict[str, list[Path]]:
    """ basename -> every file with that name under the repo's source and doc directories """
    index: dict[str, list[Path]] = {}
    roots = [REPO_ROOT / name for name in (*first_party_packages(), *_REPO_DIRS)]
    for root in roots:
        if not root.is_dir():
            continue
        for path in root.rglob("*"):
            if path.is_file() and not _SKIP_DIRS.intersection(path.parts):
                index.setdefault(path.name, []).append(path)
    return index


def _resolve_bare_filename(citation: Citation) -> Path | None:
    """ `default_settings.py` in operations.md means the one file of that name; the doc's own tree wins """
    matches = _files_by_name().get(citation.path, [])
    doc_dir = citation.doc.parent
    local = [m for m in matches if doc_dir in m.parents]
    return (local or matches or [None])[0]


@cache
def _module_symbols(path: Path) -> dict[str, set[str]]:
    """ {"": top-level names, "ClassName": member names, "Outer.Inner": ...} for one module """
    try:
        tree = ast.parse(path.read_bytes(), filename=str(path))
    except (SyntaxError, ValueError):
        return {}
    symbols: dict[str, set[str]] = {}

    def names_in(body) -> set[str]:
        names = set()
        for node in body:
            if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
                names.add(node.name)
            elif isinstance(node, ast.Assign):
                for target in node.targets:
                    names.update(_target_names(target))
            elif isinstance(node, (ast.AnnAssign, ast.AugAssign)):
                names.update(_target_names(node.target))
            elif isinstance(node, (ast.If, ast.Try, ast.With)):
                names.update(names_in(node.body))
                names.update(names_in(getattr(node, "orelse", [])))
                for handler in getattr(node, "handlers", []):
                    names.update(names_in(handler.body))
            elif isinstance(node, ast.ImportFrom):
                names.update(alias.asname or alias.name for alias in node.names)
            elif isinstance(node, ast.Import):
                names.update((alias.asname or alias.name).split(".")[0] for alias in node.names)
        return names

    def visit_class(class_node, prefix):
        members = names_in(class_node.body)
        for node in class_node.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                members.update(_self_attributes(node))
        symbols[prefix] = members
        for node in class_node.body:
            if isinstance(node, ast.ClassDef):
                visit_class(node, f"{prefix}.{node.name}")

    symbols[""] = names_in(tree.body)
    for node in tree.body:
        if isinstance(node, ast.ClassDef):
            visit_class(node, node.name)
    return symbols


def _target_names(target) -> set[str]:
    if isinstance(target, ast.Name):
        return {target.id}
    if isinstance(target, (ast.Tuple, ast.List)):
        return set().union(*(_target_names(t) for t in target.elts))
    return set()


def _self_attributes(func) -> set[str]:
    """ Attributes assigned as self.x = ... inside a method (instance attributes cited as Class.x) """
    names = set()
    for node in ast.walk(func):
        if isinstance(node, ast.Attribute) and isinstance(node.ctx, ast.Store) \
                and isinstance(node.value, ast.Name) and node.value.id == "self":
            names.add(node.attr)
    return names


def _package_symbols(package_dir: Path) -> set[str]:
    """ Names a package re-exports: anything defined or imported in its __init__ plus its submodule names """
    names = set()
    init = package_dir / "__init__.py"
    if init.exists():
        names.update(_module_symbols(init).get("", set()))
    for child in package_dir.iterdir():
        if child.suffix == ".py" and child.name != "__init__.py":
            names.add(child.stem)
        elif child.is_dir() and (child / "__init__.py").exists():
            names.add(child.name)
    return names


def resolve_symbol(path: Path, symbol: str) -> str | None:
    """ None when the dotted symbol resolves in the module at `path`; otherwise why not """
    if path.is_dir():
        head = symbol.split(".")[0]
        return None if head in _package_symbols(path) else f"no `{head}` in package {path.name}/"
    if path.suffix != ".py":
        return f"{path.name} is not a Python module"
    symbols = _module_symbols(path)
    parts = symbol.split(".")
    if parts[0] not in symbols.get("", set()):
        return f"no top-level `{parts[0]}`"
    # Walk Class.member[.member]: the member must be defined on that class (or nested class) in this file
    for depth in range(1, len(parts)):
        owner = ".".join(parts[:depth])
        if owner not in symbols:
            return f"`{owner}` is not a class in this module, so `{parts[depth]}` cannot be checked"
        if parts[depth] not in symbols[owner]:
            return f"`{owner}` has no member `{parts[depth]}` defined in this module"
    return None


@cache
def _markdown_anchors(path: Path) -> set[str]:
    anchors = set()
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("#"):
            heading = line.lstrip("#").strip()
            heading = re.sub(r"<a id=\"([\w\-]+)\"></a>", "", heading).strip()
            anchors.add(re.sub(r"[^\w\- ]", "", heading.lower()).strip().replace(" ", "-"))
        anchors.update(re.findall(r"<a id=\"([\w\-]+)\"></a>", line))
    return anchors


def check_citation(citation: Citation) -> str | None:
    if is_build_output(citation.path):
        return None
    path = resolve_path(citation)
    if path is None:
        return "no such file"
    if citation.symbol:
        return resolve_symbol(path, citation.symbol)
    if citation.anchor and path.suffix == ".md" and citation.anchor not in _markdown_anchors(path):
        return f"no heading or anchor `#{citation.anchor}` in {path.name}"
    return None


@dataclass
class DocsReport:
    dead: list[DeadCitation]
    docs_checked: int
    citations_checked: int
    unchecked_plans: dict[str, str]   # repo-relative path -> status (or "no Status line")

    @property
    def ok(self) -> bool:
        return not self.dead


def check_docs(paths=None, all_plans: bool = False) -> DocsReport:
    dead: list[DeadCitation] = []
    unchecked: dict[str, str] = {}
    checked = 0
    citations_checked = 0
    for doc in doc_files(paths):
        if not all_plans and is_plan(doc) and not is_live_plan(doc):
            unchecked[str(doc.relative_to(REPO_ROOT))] = plan_status(doc) or "no Status line"
            continue
        checked += 1
        for citation in citations_in(doc):
            citations_checked += 1
            reason = check_citation(citation)
            if reason:
                dead.append(DeadCitation(citation, reason))
    return DocsReport(dead, checked, citations_checked, unchecked)


def render_report(report: DocsReport) -> str:
    lines = []
    by_doc: dict[str, list[DeadCitation]] = {}
    for item in report.dead:
        by_doc.setdefault(str(item.citation.doc.relative_to(REPO_ROOT)), []).append(item)
    for doc, items in sorted(by_doc.items()):
        lines.append(f"{doc}: {len(items)} dead")
        lines += [f"  L{item.citation.lineno}  {item.citation.text}  - {item.reason}" for item in items]
    if report.unchecked_plans:
        lines.append(f"{len(report.unchecked_plans)} plan(s) not live (landed / superseded / no Status line), not checked; "
                     "--all-plans includes them")
    summary = f"{report.citations_checked} citations in {report.docs_checked} docs: {len(report.dead)} dead"
    lines.append(summary if report.ok else f"{summary} (fix the citation or the code it describes)")
    return "\n".join(lines)
