"""
Repository facts that need neither Django nor a database: where the repo root is, which top-level
directories are first-party Python packages, and which files git says have changed.
"""
import os
import subprocess
from functools import cache
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

# Directories under the repo root that are never source, even if they grow a .py file
_NOT_PACKAGES = {"claude", "config", "data", "media_root", "paper", "scripts", "test_data", "node_modules"}
_PACKAGE_MARKERS = ("__init__.py", "apps.py", "urls.py", "models.py", "models")


@cache
def first_party_packages() -> tuple[str, ...]:
    """ Top-level importable packages in this repo (Django apps plus library/ and variantgrid/) """
    packages = []
    for entry in sorted(REPO_ROOT.iterdir()):
        if not entry.is_dir() or entry.name.startswith(".") or entry.name in _NOT_PACKAGES:
            continue
        if any((entry / marker).exists() for marker in _PACKAGE_MARKERS):
            packages.append(entry.name)
    return tuple(packages)


def iter_python_files(packages=None):
    """ Yields (module_name, path) for every .py under the given (default: all first-party) packages """
    for package in packages or first_party_packages():
        for dirpath, dirnames, filenames in os.walk(REPO_ROOT / package):
            dirnames[:] = sorted(d for d in dirnames if d != "__pycache__" and not d.startswith("."))
            for filename in sorted(filenames):
                if filename.endswith(".py"):
                    path = Path(dirpath) / filename
                    yield module_name_for(path), path


def module_name_for(path: Path) -> str | None:
    """ 'snpdb/models/models_variant.py' -> 'snpdb.models.models_variant'; packages drop '__init__' """
    try:
        relative = Path(path).resolve().relative_to(REPO_ROOT)
    except ValueError:
        return None
    if relative.suffix != ".py":
        return None
    parts = list(relative.with_suffix("").parts)
    if parts[-1] == "__init__":
        parts.pop()
    if not parts or parts[0] not in first_party_packages():
        return None
    return ".".join(parts)


def app_for_path(path) -> str | None:
    """ Top-level package a repo-relative or absolute path belongs to, or None """
    path = Path(path)
    if path.is_absolute():
        try:
            path = path.relative_to(REPO_ROOT)
        except ValueError:
            return None
    top = path.parts[0] if path.parts else None
    return top if top in first_party_packages() else None


def git(*args) -> str:
    return subprocess.run(["git", *args], cwd=REPO_ROOT, check=True, capture_output=True, text=True).stdout


def git_short_sha() -> str:
    return git("rev-parse", "--short", "HEAD").strip()


def default_base_ref() -> str:
    """ The merge base with origin/master when it exists (so a branch compares against what it forked
        from), otherwise HEAD (so only uncommitted work counts) """
    try:
        return git("merge-base", "HEAD", "origin/master").strip()
    except subprocess.CalledProcessError:
        return "HEAD"


def changed_files(base: str | None = None) -> list[str]:
    """ Repo-relative paths changed since `base`: committed, staged, unstaged and untracked """
    base = base or default_base_ref()
    changed = set(git("diff", "--name-only", base).splitlines())
    changed.update(git("ls-files", "--others", "--exclude-standard").splitlines())
    return sorted(path for path in changed if path and (REPO_ROOT / path).exists())
