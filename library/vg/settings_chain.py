"""
`vg settings`: which settings file the current box resolves to, the order its `from … import *` chain
is evaluated in, and every module-level assignment to a given NAME along that chain - so "where did
this value come from" is one command instead of a walk through variantgrid/settings/.

Everything here is AST over the settings package; nothing is imported, so secrets are never
evaluated. The resolved value itself needs django.conf.settings and is read by the command.
"""
import ast
import os
import socket
from dataclasses import dataclass
from pathlib import Path

from library.vg.repo import REPO_ROOT

SETTINGS_DIR = REPO_ROOT / "variantgrid" / "settings"
SETTINGS_PACKAGE = "variantgrid.settings"
HOSTNAME_SETTINGS_MODULE = "variantgrid.settings"  # the package's __init__ picks the env module from the hostname
_MUTATING_METHODS = ("update", "append", "extend", "setdefault", "pop", "insert", "remove")


@dataclass
class Assignment:
    path: Path
    lineno: int
    mutated: bool     # `X[...] = ` / `X.update(...)` rather than `X = `
    source: str       # the statement's first line, as written

    @property
    def label(self) -> str:
        return f"{self.path.relative_to(REPO_ROOT)}:{self.lineno}"


def flattened_hostname(hostname: str | None = None) -> str:
    """ The env module stem variantgrid/settings/__init__.py derives from a hostname: 'vg-test2' -> 'vgtest2' """
    stem = (hostname or socket.gethostname()).lower().split(".")[0].replace("-", "")
    return f"s{stem}" if stem[0].isnumeric() else stem


def resolved_settings_module(environ=None) -> str:
    """ The module DJANGO_SETTINGS_MODULE ends up meaning on this box, after the hostname lookup """
    environ = os.environ if environ is None else environ
    module = environ.get("DJANGO_SETTINGS_MODULE") or HOSTNAME_SETTINGS_MODULE
    if module != HOSTNAME_SETTINGS_MODULE:
        return module
    stem = flattened_hostname()
    for folder in ("env_developers", "env"):
        if (SETTINGS_DIR / folder / f"{stem}.py").exists():
            return f"{SETTINGS_PACKAGE}.{folder}.{stem}"
    return f"{SETTINGS_PACKAGE}.env.{stem} (missing - settings would fail to load)"


def module_path(module: str) -> Path | None:
    parts = module.split(".")
    if parts[:2] != ["variantgrid", "settings"]:
        return None
    path = REPO_ROOT.joinpath(*parts).with_suffix(".py")
    return path if path.exists() else None


def _star_imports(path: Path) -> list[str]:
    try:
        tree = ast.parse(path.read_bytes())
    except (OSError, SyntaxError):
        return []
    return [node.module for node in tree.body
            if isinstance(node, ast.ImportFrom) and node.module and node.level == 0
            and any(alias.name == "*" for alias in node.names)]


def settings_chain(module: str) -> list[Path]:
    """ Files whose module-level assignments land in the settings namespace, in evaluation order: a
        `from X import *` runs X (and X's own star-imports) to completion before the importer continues,
        and a module already executed is never re-run, so this is a post-order walk with a seen set. """
    chain: list[Path] = []
    seen: set[Path] = set()

    def visit(current: str):
        path = module_path(current)
        if path is None or path in seen:
            return
        seen.add(path)
        for imported in _star_imports(path):
            visit(imported)
        chain.append(path)

    visit(module)
    return chain


def _module_level_statements(body):
    """ Statements that run at import time: the module body plus the branches of its if / try / with blocks
        (settings files gate assignments on UNIT_TEST and on optional imports), never function or class bodies """
    for node in body:
        yield node
        if isinstance(node, (ast.If, ast.With)):
            yield from _module_level_statements(node.body)
            yield from _module_level_statements(getattr(node, "orelse", []))
        elif isinstance(node, ast.Try):
            yield from _module_level_statements(node.body)
            for handler in node.handlers:
                yield from _module_level_statements(handler.body)
            yield from _module_level_statements(node.orelse)
            yield from _module_level_statements(node.finalbody)


def assignments(path: Path) -> list[tuple[str, Assignment]]:
    """ (setting name, Assignment) for every module-level statement that assigns or mutates a name """
    try:
        source_lines = path.read_text().splitlines()
        tree = ast.parse("\n".join(source_lines))
    except (OSError, SyntaxError):
        return []
    found = []
    for node in _module_level_statements(tree.body):
        targets, mutated = [], False
        if isinstance(node, ast.Assign):
            targets = node.targets
        elif isinstance(node, (ast.AugAssign, ast.AnnAssign)):
            targets = [node.target]
        elif isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
            func = node.value.func
            if isinstance(func, ast.Attribute) and func.attr in _MUTATING_METHODS:
                targets, mutated = [func.value], True
        for target in targets:
            target_mutated = mutated
            while isinstance(target, (ast.Subscript, ast.Attribute)):
                target = target.value
                target_mutated = True
            names = [e for e in target.elts if isinstance(e, ast.Name)] if isinstance(target, ast.Tuple) else \
                [target] if isinstance(target, ast.Name) else []
            for name in names:
                found.append((name.id, Assignment(path, node.lineno, target_mutated, source_lines[node.lineno - 1].strip())))
    return found


def trail(name: str, module: str) -> list[Assignment]:
    """ Every assignment to `name` along the chain for `module`, in the order they take effect. A name
        only ever imported by name (settings_paths constants) is looked up across components/ as well. """
    chain = settings_chain(module)
    found = [assignment for path in chain for assigned, assignment in assignments(path) if assigned == name]
    if not found:
        for path in sorted((SETTINGS_DIR / "components").glob("*.py")):
            if path not in chain:
                found += [assignment for assigned, assignment in assignments(path) if assigned == name]
    return found


def names_assigned_in(path: Path) -> list[str]:
    """ Distinct setting-shaped names a file assigns or mutates, in file order """
    seen: list[str] = []
    for name, _ in assignments(path):
        if name.isupper() and name not in seen:
            seen.append(name)
    return seen
