"""
`vg map`: generated facts about the codebase, committed under claude/maps/ and checked in CI.

Each generator returns a list of MapTable; `write_map` renders it to claude/maps/<name>.md and `check`
reports the maps whose committed copy differs from a fresh render (the `makemigrations --check` idea).

Maps are canonical under the CI settings module (variantgrid.settings.env.github_actions) because
`models` and `urls` depend on INSTALLED_APPS and settings-gated URL includes. scripts/vg sets that
module for you; `manage.py vg map` under another settings module warns.
"""
from importlib import import_module

from library.vg.markdown import render_json, render_markdown
from library.vg.repo import REPO_ROOT

MAPS_DIR = REPO_ROOT / "claude" / "maps"
CANONICAL_SETTINGS_MODULE = "variantgrid.settings.env.github_actions"

# name -> needs Django. Generator modules are imported on demand (library.vg.maps.<name>) so the
# static ones can run from scripts/vg without Django on the path.
MAP_GENERATORS: dict[str, bool] = {
    "models": True,
    "urls": True,
    "commands": True,
    "tasks": False,
    "signals": False,
    "settings": False,
}

STATIC_MAPS = tuple(name for name, needs_django in MAP_GENERATORS.items() if not needs_django)


def render(name: str, as_json: bool = False, **kwargs) -> str:
    generator = import_module(f"library.vg.maps.{name}")
    tables = generator.generate(**kwargs)
    if as_json:
        return render_json(tables)
    return render_markdown(f"{name} map", generator.INTRO, tables)


def map_path(name: str):
    return MAPS_DIR / f"{name}.md"


def write_map(name: str, **kwargs) -> bool:
    """ Writes claude/maps/<name>.md; returns True if the content changed """
    content = render(name, **kwargs)
    path = map_path(name)
    changed = not path.exists() or path.read_text() != content
    if changed:
        MAPS_DIR.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
    return changed


def check(names) -> list[str]:
    """ Names of maps whose committed copy is missing or differs from a fresh render """
    stale = []
    for name in names:
        path = map_path(name)
        if not path.exists() or path.read_text() != render(name):
            stale.append(name)
    return stale
