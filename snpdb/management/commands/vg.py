"""
manage.py vg — introspection for driving VariantGrid from an agent or a terminal.

    vg status [--json]
    vg map [models|urls|commands|tasks|signals|settings|all] [--json] [--app X] [--counts]
    vg tests [--changed] [--base REF] [--run] [--parallel N] [--explain]
    vg page <url-or-url-name> [--as USER] [--kwargs k=v ...] [--text|--html|--links|--forms|--json] [--queries]
    vg page --create-user
    vg outline <file.py> [--min-lines N] | vg outline --coverage [package ...]
    vg settings [NAME] [--diff] [--json]
    vg docs check [doc.md|dir ...] [--all-plans]
    vg inspect <kind> <id> [--depth N] [--json]

All logic lives in library/vg/; this file only parses arguments. See claude/plans/agent_system.md §4.2.
"""
import json
import os
import re
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path

from django.conf import settings
from django.core.management.base import BaseCommand, CommandError

from library.vg import maps
from library.vg.docs import check_docs, render_report
from library.vg.inspect import KINDS, inspect, render_inspection
from library.vg.outline import outline, render_coverage, render_outline
from library.vg.page import AgentUserMissing, create_agent_user, render_page
from library.vg.repo import REPO_ROOT
from library.vg.settings_chain import (
    names_assigned_in,
    resolved_settings_module,
    settings_chain,
    trail,
)
from library.vg.status import gather_status, render_status
from library.vg.test_selection import select_tests

MAP_CHOICES = [*maps.MAP_GENERATORS, "all"]


class Command(BaseCommand):
    help = "Agent introspection: box status, generated maps, changed-file test selection, page rendering, outlines, settings"
    category = "dev"

    def add_arguments(self, parser):
        subparsers = parser.add_subparsers(dest="subcommand", required=True)

        status_parser = subparsers.add_parser("status", help="What is running on this box: db, builds, services, queues, errors")
        status_parser.add_argument("--json", action="store_true")

        map_parser = subparsers.add_parser("map", help="Generate claude/maps/*.md")
        map_parser.add_argument("name", nargs="?", default="all", choices=MAP_CHOICES)
        map_parser.add_argument("--json", action="store_true", help="Print JSON to stdout instead of writing")
        map_parser.add_argument("--print", action="store_true", help="Print Markdown to stdout instead of writing")
        map_parser.add_argument("--app", help="Restrict to one app (with --json/--print only)")
        map_parser.add_argument("--counts", action="store_true", help="models: add pg_class row estimates (queries)")

        tests_parser = subparsers.add_parser("tests", help="Select test modules for changed files")
        tests_parser.add_argument("--changed", action="store_true", default=True,
                                  help="Select from git changes (default and only mode for now)")
        tests_parser.add_argument("--base", help="Git ref to diff against (default: merge-base with origin/master)")
        tests_parser.add_argument("--run", action="store_true", help="Run the selected tests with --keepdb")
        tests_parser.add_argument("--parallel", type=int, help="Pass --parallel N to the test runner")
        tests_parser.add_argument("--explain", action="store_true", help="Show which changed file selected each label")

        page_parser = subparsers.add_parser("page", help="Render a page through the test client as claude_agent")
        page_parser.add_argument("url", nargs="?", help="A path (/variants/…) or URL name (view_variant)")
        page_parser.add_argument("--as", dest="username", help="Render as this user instead of claude_agent")
        page_parser.add_argument("--kwargs", nargs="*", default=[], metavar="k=v", help="URL name kwargs")
        page_parser.add_argument("--queries", action="store_true", help="Report production query count and repeats")
        page_parser.add_argument("--create-user", action="store_true", help="Create the claude_agent user and exit")
        output = page_parser.add_mutually_exclusive_group()
        output.add_argument("--text", action="store_true", help="Page as plain text")
        output.add_argument("--html", action="store_true", help="Raw HTML")
        output.add_argument("--links", action="store_true", help="Links as 'text → href'")
        output.add_argument("--forms", action="store_true", help="Forms and their fields")
        output.add_argument("--json", action="store_true", help="Everything as JSON")

        outline_parser = subparsers.add_parser("outline", help="Classes/functions of a module with line numbers and doc lines")
        outline_parser.add_argument("file", nargs="*", help="Module to outline (or packages, with --coverage)")
        outline_parser.add_argument("--min-lines", type=int, default=0, help="Hide methods shorter than this")
        outline_parser.add_argument("--coverage", action="store_true", help="Module docstring coverage per package")

        docs_parser = subparsers.add_parser("docs", help="Check that every path / path:Symbol citation in the docs resolves")
        docs_parser.add_argument("action", choices=["check"])
        docs_parser.add_argument("paths", nargs="*", help="Docs or directories to check (default: all agent docs)")
        docs_parser.add_argument("--all-plans", action="store_true", help="Also check landed / superseded / unstatused plans")

        inspect_parser = subparsers.add_parser("inspect", help="An object's whole graph by domain kind: " + ", ".join(KINDS))
        inspect_parser.add_argument("kind", choices=list(KINDS))
        inspect_parser.add_argument("key", help="pk, or a natural key (CA id, gene symbol, transcript accession, username, lab group name)")
        inspect_parser.add_argument("--depth", type=int, default=2, help="How far to follow relations (default 2)")
        inspect_parser.add_argument("--json", action="store_true")

        settings_parser = subparsers.add_parser("settings", help="Resolved value of a setting and every file that assigned it")
        settings_parser.add_argument("name", nargs="?", help="Setting name; omit to list the settings chain for this box")
        settings_parser.add_argument("--diff", action="store_true", help="List the settings this box's env file overrides")
        settings_parser.add_argument("--json", action="store_true")

    def handle(self, *args, **options):
        getattr(self, f"handle_{options['subcommand']}")(**options)

    # --- status ---

    def handle_status(self, json: bool, **_):  # pylint: disable=redefined-outer-name
        status = gather_status()
        self.stdout.write(_json_dumps(asdict(status)) if json else render_status(status))

    # --- map ---

    def handle_map(self, name, json: bool, print: bool, app, counts, **_):  # pylint: disable=redefined-outer-name
        names = list(maps.MAP_GENERATORS) if name == "all" else [name]
        settings_module = os.environ.get("DJANGO_SETTINGS_MODULE")
        if settings_module != maps.CANONICAL_SETTINGS_MODULE and any(maps.MAP_GENERATORS[n] for n in names):
            self.stderr.write(f"Note: maps are canonical under {maps.CANONICAL_SETTINGS_MODULE} "
                              f"(settings-gated apps/URLs differ); use scripts/vg map to generate them that way.")
        kwargs = {"app": app} if app else {}
        for map_name in names:
            map_kwargs = {**kwargs, "counts": counts} if map_name == "models" else kwargs
            if json or print:
                self.stdout.write(maps.render(map_name, as_json=json, **map_kwargs))
            else:
                changed = maps.write_map(map_name, **map_kwargs)
                self.stdout.write(f"{maps.map_path(map_name).relative_to(REPO_ROOT)}: {'updated' if changed else 'unchanged'}")

    # --- tests ---

    def handle_tests(self, base, run, parallel, explain, **_):
        selection = select_tests(base=base)
        if not selection.changed:
            self.stdout.write("No changed files.")
            return
        if explain:
            for label in selection.labels:
                self.stdout.write(f"{label}  <- {', '.join(selection.reasons[label])}")
            self.stdout.write("")
        if not selection.labels:
            self.stdout.write(f"{len(selection.changed)} changed file(s) select no tests.")
            return
        test_args = ["test", "--keepdb", *selection.labels]
        if parallel:
            test_args += ["--parallel", str(parallel)]
        self.stdout.write(f"{len(selection.labels)} test label(s) for {len(selection.changed)} changed file(s):")
        self.stdout.write("python3 manage.py " + " ".join(test_args))
        if run:
            self.stdout.flush()
            result = subprocess.run([sys.executable, str(REPO_ROOT / "manage.py"), *test_args], cwd=REPO_ROOT, check=False)
            if result.returncode:
                raise CommandError(f"Tests failed (exit {result.returncode})")

    # --- page ---

    def handle_page(self, url, username, kwargs, queries, create_user, text, html, links, forms, json: bool, **_):  # pylint: disable=redefined-outer-name
        if create_user:
            user, created = create_agent_user()
            self.stdout.write(f"User {user.username} {'created' if created else 'already existed'}; in all_users.")
            return
        if not url:
            raise CommandError("Give a path or URL name, e.g. `vg page /variantopedia/dashboard` or `vg page view_variant --kwargs variant_id=1`")
        url_kwargs = dict(pair.split("=", 1) for pair in kwargs)
        try:
            result = render_page(url, username=username, kwargs=url_kwargs, queries=queries)
        except (AgentUserMissing, ValueError) as e:
            raise CommandError(str(e)) from e

        if json:
            self.stdout.write(_json_dumps({k: v for k, v in asdict(result).items() if k != "html"}))
        elif html:
            self.stdout.write(result.html)
        elif text:
            self.stdout.write(result.text)
        elif links:
            self.stdout.write("\n".join(f"{t or '(no text)'} -> {h}" for t, h in result.links))
        elif forms:
            self.stdout.write("\n".join(result.forms) or "(no forms)")
        else:
            self.stdout.write(f"GET {result.url} -> {result.status}" + (f" -> {result.redirect_to}" if result.redirect_to else ""))
            if result.templates:
                self.stdout.write("templates: " + ", ".join(result.templates))
            for line in result.outline:
                self.stdout.write("  " + line)
        if queries and not json:
            self.stdout.write(f"queries: {result.production_queries} production ({result.total_queries} total)")
            for sql, n in result.repeated[:10]:
                self.stdout.write(f"  x{n}: {sql[:200]}")

    # --- outline ---

    def handle_outline(self, file, min_lines, coverage, **_):
        if coverage:
            self.stdout.write(render_coverage(file or None))
            return
        if not file:
            raise CommandError("Give a module to outline, or --coverage")
        path = Path(file[0])
        if not path.exists():
            raise CommandError(f"No such file: {file[0]}")
        self.stdout.write(render_outline(path.relative_to(REPO_ROOT) if path.is_absolute() else path, outline(path), min_lines))

    # --- docs ---

    def handle_docs(self, paths, all_plans, **_):
        report = check_docs(paths or None, all_plans=all_plans)
        self.stdout.write(render_report(report))
        if not report.ok:
            raise CommandError(f"{len(report.dead)} dead citation(s)")

    # --- inspect ---

    def handle_inspect(self, kind, key, depth, json: bool, **_):  # pylint: disable=redefined-outer-name
        try:
            inspection = inspect(kind, key, depth=depth)
        except LookupError as e:
            raise CommandError(str(e)) from e
        self.stdout.write(_json_dumps(inspection) if json else render_inspection(inspection))

    # --- settings ---

    def handle_settings(self, name, diff, json: bool, **_):  # pylint: disable=redefined-outer-name
        module = resolved_settings_module()
        chain = settings_chain(module)
        if not chain:
            raise CommandError(f"{module} is not a file under variantgrid/settings/, so there is no chain to walk")
        if diff:
            env_path = chain[-1]
            overridden = names_assigned_in(env_path)
            if json:
                self.stdout.write(_json_dumps({"module": module, "file": str(env_path.relative_to(REPO_ROOT)),
                                               "overrides": {n: _setting_repr(n) for n in overridden}}))
            else:
                self.stdout.write(f"{env_path.relative_to(REPO_ROOT)} sets or mutates {len(overridden)} settings:")
                for setting_name in overridden:
                    self.stdout.write(f"  {setting_name} = {_setting_repr(setting_name)}")
            return
        if not name:
            self.stdout.write(f"DJANGO_SETTINGS_MODULE resolves to {module}; assignments take effect in this order:")
            for path in chain:
                self.stdout.write(f"  {path.relative_to(REPO_ROOT)}")
            return
        name = name.upper()
        found = trail(name, module)
        if json:
            self.stdout.write(_json_dumps({"name": name, "value": _setting_repr(name),
                                           "trail": [{"file": a.label, "mutated": a.mutated, "source": a.source} for a in found]}))
            return
        self.stdout.write(f"{name} = {_setting_repr(name)}")
        if not found:
            self.stdout.write("  (no assignment found in variantgrid/settings/ - a Django default, or set outside the package)")
        for assignment in found:
            self.stdout.write(f"  {assignment.label}  {'~ ' if assignment.mutated else ''}{assignment.source}")


_SECRET_MARKERS = ("SECRET", "PASSWORD", "TOKEN", "KEY", "CREDENTIAL")
_URL_CREDENTIALS_RE = re.compile(r"://[^/@\s]+:[^/@\s]+@")


def _is_secret_name(name: str) -> bool:
    return any(marker in name.upper() for marker in _SECRET_MARKERS)


def _mask(value):
    """ Hides secret-looking dict keys at any depth and credentials embedded in URLs """
    if isinstance(value, dict):
        return {k: "(hidden)" if isinstance(k, str) and _is_secret_name(k) else _mask(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return type(value)(_mask(v) for v in value)
    if isinstance(value, str):
        return _URL_CREDENTIALS_RE.sub("://***@", value)
    return value


def _setting_repr(name: str, limit: int = 400) -> str:
    if not hasattr(settings, name):
        return "(not set)"
    if _is_secret_name(name):
        return "(secret, hidden)"
    value = repr(_mask(getattr(settings, name)))
    return value if len(value) <= limit else value[:limit] + f"… ({len(value)} chars)"


def _json_dumps(data) -> str:
    return json.dumps(data, indent=1, default=str)
