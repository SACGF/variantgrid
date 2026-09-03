"""
manage.py vg — introspection for driving VariantGrid from an agent or a terminal.

    vg map [models|urls|commands|tasks|signals|settings|all] [--check] [--json] [--app X] [--counts]
    vg tests [--changed] [--base REF] [--run] [--parallel N] [--explain]
    vg page <url-or-url-name> [--as USER] [--kwargs k=v ...] [--text|--html|--links|--forms|--json] [--queries]
    vg page --create-user

All logic lives in library/vg/; this file only parses arguments. See claude/plans/agent_system.md §4.2.
"""
import json
import os
import subprocess
import sys
from dataclasses import asdict

from django.core.management.base import BaseCommand, CommandError

from library.vg import maps
from library.vg.page import AgentUserMissing, create_agent_user, render_page
from library.vg.repo import REPO_ROOT
from library.vg.test_selection import select_tests

MAP_CHOICES = [*maps.MAP_GENERATORS, "all"]


class Command(BaseCommand):
    help = "Agent introspection: generated maps, changed-file test selection, page rendering"
    category = "dev"

    def add_arguments(self, parser):
        subparsers = parser.add_subparsers(dest="subcommand", required=True)

        map_parser = subparsers.add_parser("map", help="Generate claude/maps/*.md (or --check them)")
        map_parser.add_argument("name", nargs="?", default="all", choices=MAP_CHOICES)
        map_parser.add_argument("--check", action="store_true", help="Exit 1 if any committed map is stale")
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

    def handle(self, *args, **options):
        getattr(self, f"handle_{options['subcommand']}")(**options)

    # --- map ---

    def handle_map(self, name, check, json: bool, print: bool, app, counts, **_):  # pylint: disable=redefined-outer-name
        names = list(maps.MAP_GENERATORS) if name == "all" else [name]
        settings_module = os.environ.get("DJANGO_SETTINGS_MODULE")
        if settings_module != maps.CANONICAL_SETTINGS_MODULE and any(maps.MAP_GENERATORS[n] for n in names):
            self.stderr.write(f"Note: maps are canonical under {maps.CANONICAL_SETTINGS_MODULE} "
                              f"(settings-gated apps/URLs differ); use scripts/vg map to generate them that way.")
        kwargs = {"app": app} if app else {}
        if check:
            stale = maps.check(names)
            if stale:
                raise CommandError(f"Stale maps: {', '.join(stale)} — run `scripts/vg map` and commit claude/maps/")
            self.stdout.write(f"Maps up to date: {', '.join(names)}")
            return
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


def _json_dumps(data) -> str:
    return json.dumps(data, indent=1, default=str)
