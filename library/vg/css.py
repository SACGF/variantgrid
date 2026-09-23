"""
`vg css unused`: the class and id selectors in our SCSS that nothing in the templates, JS or Python
mentions, so dead styling announces itself before a stylesheet is split or a rule is copied.

Owns: which stylesheets are ours (SCSS_GLOBS - the hand-written .scss under static_files), how a
selector's classes and ids are read out of a stylesheet (the text before each `{`, back to the previous
brace or semicolon, so nesting and one-line rules both count), and what counts as a mention. A mention is the whole name as a
token (`cs-r` is not mentioned by `cs-report`) in any .html / .js / .py / .jinja file under the first-party
packages, minified vendor JS included, so a class a library adds at runtime (`dataTables_filter`,
`select2-selection`) counts as used because its script names it. Built-up names are the trap: a template
writes `cs-{{ status }}` and the stylesheet lists `cs-vus_a`, `cs-withdrawn` and so on, so a name is
reported as *dynamic*, not unused, when any proper prefix of it that ends in `-` or `_` is itself a token
somewhere (`cs-` is a token of `"cs-" + status`, `cs-{{ x }}`, `f"cs-{x}"` and `cs-${x}` alike). Those
need a reader's eye, which is why they are listed (`--dynamic`) rather than dropped. A class added by code
that is not in the tree at all - Django's form `errorlist`, jQuery UI's `ui-*`, FontAwesome's SVG output -
is listed in LIBRARY_ADDED and reported as such, never as unused.

Entry points: `find_unused` (returns an UnusedReport) and `render_unused_report`. Pure filesystem - no Django.
"""
import re
from dataclasses import dataclass, field
from pathlib import Path

from library.vg.repo import REPO_ROOT, first_party_packages

SCSS_GLOBS = ("variantgrid/static_files/*_static/css/*.scss",)
SOURCE_SUFFIXES = (".html", ".js", ".py", ".jinja", ".jinja2")
# Directories under a package that are build output or vendored trees, never a source of class names
SKIPPED_DIRS = {"sitestatic", "node_modules", "__pycache__"}
# This module and its tests talk about selectors without using any, so they are not sources
SELF = {"library/vg/css.py", "library/tests/test_vg.py"}
# Class name prefixes added at runtime by code that lives outside the repo (CDN scripts, Django itself)
LIBRARY_ADDED = (
    "errorlist",  # Django form errors
    "ui-",  # jQuery UI widgets (CDN)
    "svg-inline--fa",  # FontAwesome's SVG rendering (CDN)
    "fileupload-", "template-download", "template-upload",  # jQuery File Upload
)

_TOKEN = re.compile(r"[\w-]+")
_CLASS = re.compile(r"\.(-?[A-Za-z_][\w-]*)")
_ID = re.compile(r"#([A-Za-z_][\w-]*)")


@dataclass
class Selector:
    name: str
    kind: str  # "class" or "id"
    locations: list[tuple[str, int]] = field(default_factory=list)  # (repo-relative stylesheet, line)

    @property
    def spelled(self) -> str:
        return ("." if self.kind == "class" else "#") + self.name


@dataclass
class UnusedReport:
    stylesheets: list[str]
    selectors: int
    unused: list[Selector]
    dynamic: list[Selector]  # a prefix is built up somewhere, so a reader must decide
    library: list[Selector]  # named only by code outside the tree (LIBRARY_ADDED)

    @property
    def ok(self) -> bool:
        return not self.unused


def stylesheets() -> list[Path]:
    paths = []
    for pattern in SCSS_GLOBS:
        paths.extend(REPO_ROOT.glob(pattern))
    return sorted(paths)


_COMMENT = re.compile(r"/\*.*?\*/|(?<![:\w])//[^\n]*", re.DOTALL)


def _selector_openers(text: str):
    """ Yields (line_number, selector text) for every `{`: the text since the previous `{`, `}` or `;` opens it """
    # Comments are blanked, not cut, so line numbers survive; `//` after a colon is a URL, not a comment
    text = _COMMENT.sub(lambda m: re.sub(r"[^\n]", " ", m.group()), text)
    for match in re.finditer(r"[^{};]*\{", text):
        selector = match.group()[:-1].strip()
        if not selector or selector.startswith(("@", "$")):
            continue
        yield text.count("\n", 0, match.end()) + 1, selector


def selectors_in(path: Path, into: dict[tuple[str, str], Selector] | None = None) -> dict[tuple[str, str], Selector]:
    """ Class and id selectors of one stylesheet, keyed (kind, name), merged into `into` when given """
    relative = str(path.relative_to(REPO_ROOT))
    found = into if into is not None else {}
    for number, selector in _selector_openers(path.read_text()):
        for pattern, kind in ((_CLASS, "class"), (_ID, "id")):
            for name in pattern.findall(selector):
                found.setdefault((kind, name), Selector(name, kind)).locations.append((relative, number))
    return found


def _source_files():
    for package in first_party_packages():
        for path in sorted((REPO_ROOT / package).rglob("*")):
            if path.suffix not in SOURCE_SUFFIXES or SKIPPED_DIRS & set(path.parts):
                continue
            if str(path.relative_to(REPO_ROOT)) in SELF:
                continue
            yield path


def source_tokens(files=None) -> set[str]:
    """ Every [\\w-]+ token in the source tree: class names, and the `cs-` stem of a built-up `cs-{{ x }}` """
    tokens: set[str] = set()
    for path in files if files is not None else _source_files():
        tokens.update(_TOKEN.findall(path.read_text(errors="replace")))
    return tokens


def _dynamic_prefixes(name: str):
    for i, char in enumerate(name[:-1], start=1):
        if char in "-_":
            yield name[:i]


def classify(name: str, tokens: set[str]) -> str:
    """ 'used', 'library', 'dynamic' (a `-`/`_` prefix of the name is a token) or 'unused' """
    if name.startswith(LIBRARY_ADDED):
        return "library"
    if name in tokens:
        return "used"
    if any(prefix in tokens for prefix in _dynamic_prefixes(name)):
        return "dynamic"
    return "unused"


def find_unused(paths=None) -> UnusedReport:
    sheets = [Path(p).resolve() for p in paths] if paths else stylesheets()
    selectors: dict[tuple[str, str], Selector] = {}
    for sheet in sheets:
        selectors_in(sheet, into=selectors)
    tokens = source_tokens()
    buckets: dict[str, list[Selector]] = {"unused": [], "dynamic": [], "library": []}
    for selector in sorted(selectors.values(), key=lambda s: s.locations[0]):
        verdict = classify(selector.name, tokens)
        if verdict in buckets:
            buckets[verdict].append(selector)
    return UnusedReport(
        stylesheets=[str(s.relative_to(REPO_ROOT)) for s in sheets],
        selectors=len(selectors),
        **buckets,
    )


def render_unused_report(report: UnusedReport, show_dynamic: bool = False) -> str:
    sheets = report.stylesheets[0] if len(report.stylesheets) == 1 else f"{len(report.stylesheets)} stylesheets"
    out = [f"{report.selectors} class/id selectors in {sheets}"]
    if report.unused:
        out.append(f"\n{len(report.unused)} unused - nothing in the templates, JS or Python names them:")
        out.extend(_render_selector(s) for s in report.unused)
    else:
        out.append("\nNo unused selectors.")
    if report.dynamic:
        out.append(f"\n{len(report.dynamic)} dynamic - a prefix is built up at runtime, so check by hand"
                   + (":" if show_dynamic else " (--dynamic lists them)"))
        if show_dynamic:
            out.extend(_render_selector(s) for s in report.dynamic)
    if report.library:
        out.append(f"\n{len(report.library)} named only by a library outside the tree, kept: "
                   + " ".join(s.spelled for s in report.library))
    return "\n".join(out)


def _render_selector(selector: Selector) -> str:
    where = " ".join(f"{sheet}:{line}" for sheet, line in selector.locations)
    return f"  {selector.spelled:40} {where}"
