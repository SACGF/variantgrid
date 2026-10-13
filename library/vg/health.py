"""
`vg health`: the ratchet numbers in one place, so a monthly look (or a before/after on a change) is one call.

Each number should only move one way: module docstring coverage up; dead doc citations, one-off management
commands, pylint messages (the errors listed, heading for zero) and deprecated calls down; the test suite's wall
time down. Entry points are `gather_health` and `render_health`. No Django: coverage, citations and command categories
are read from the tree by AST; the lint counts come from the last `scripts/linting/run_pylint.sh` run (`lint.txt`,
gitignored, so it is as old as that run - CI reports the errors fresh on every push, in `agent-maps.yml`);
deprecations are `library/vg/deprecations.py:find_deprecations` (basedpyright plus `manage.py check`, ~25s,
skipped with --no-deprecations); the suite time is the "Run unit tests" step of the last green Django Tests run on master, through `gh`.
A number that can't be had (no lint.txt, no `gh` or network) is reported as unavailable rather than failing.
"""
import ast
import json
import re
import subprocess
import warnings
from collections import Counter
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from library.vg.deprecations import find_deprecations
from library.vg.docs import check_docs
from library.vg.outline import docstring_coverage
from library.vg.repo import REPO_ROOT, first_party_packages, git_short_sha

LINT_REPORT = REPO_ROOT / "lint.txt"
CI_WORKFLOW = "django-tests.yml"
CI_TEST_STEP = "Run unit tests"
UNDECLARED = "(undeclared)"
LISTED_LINT_ERRORS = 20

_PYLINT_MESSAGE_RE = re.compile(r"^\S+:\d+:\d+: ([CRWEFI])\d{4}")
_PYLINT_SCORE_RE = re.compile(r"rated at (-?[\d.]+)/10")


@dataclass
class LintHealth:
    messages: int
    by_type: dict[str, int]
    score: str | None
    age_days: int
    errors: list[str] = field(default_factory=list)


@dataclass
class SuiteHealth:
    seconds: int
    run_id: int
    head_sha: str
    started_at: str


@dataclass
class DeprecationHealth:
    external: int
    ours: int
    import_in_repo: int
    import_in_packages: int


@dataclass
class Health:
    sha: str
    documented_modules: int
    counted_modules: int
    dead_citations: int
    docs_checked: int
    dead_in_plans: int
    commands_by_category: dict[str, int]
    lint: LintHealth | None
    deprecations: DeprecationHealth | None
    suite: SuiteHealth | None
    unavailable: dict[str, str] = field(default_factory=dict)

    @property
    def one_off_commands(self) -> int:
        return self.commands_by_category.get("one-off", 0)


def _command_category(path: Path) -> str | None:
    """ The `category` a module's `Command` class declares; None when the module defines no Command """
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", SyntaxWarning)  # an invalid escape in a command's string is not ours to report
            tree = ast.parse(path.read_bytes())
    except (OSError, SyntaxError):
        return None
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and node.name == "Command":
            for statement in node.body:
                if (isinstance(statement, ast.Assign) and isinstance(statement.value, ast.Constant)
                        and any(isinstance(t, ast.Name) and t.id == "category" for t in statement.targets)):
                    return statement.value.value
            return UNDECLARED
    return None


def command_categories() -> dict[str, int]:
    counts = Counter()
    for package in first_party_packages():
        for path in (REPO_ROOT / package / "management" / "commands").glob("*.py"):
            if path.name.startswith("_"):
                continue
            if category := _command_category(path):
                counts[category] += 1
    return dict(counts.most_common())


def lint_health(report: Path = LINT_REPORT) -> LintHealth | None:
    if not report.exists():
        return None
    text = report.read_text(errors="replace")
    by_type = Counter()
    errors = []
    for line in text.splitlines():
        if m := _PYLINT_MESSAGE_RE.match(line):
            by_type[m.group(1)] += 1
            if m.group(1) == "E":
                errors.append(line)
    score = _PYLINT_SCORE_RE.search(text)
    modified = datetime.fromtimestamp(report.stat().st_mtime, tz=UTC)
    return LintHealth(messages=sum(by_type.values()), by_type=dict(by_type.most_common()),
                      score=score.group(1) if score else None,
                      age_days=(datetime.now(tz=UTC) - modified).days, errors=errors)


def _gh_api(endpoint: str):
    result = subprocess.run(["gh", "api", endpoint], cwd=REPO_ROOT, check=True, capture_output=True, text=True,
                            timeout=20)
    return json.loads(result.stdout)


def _parse_github_time(value: str) -> datetime:
    return datetime.fromisoformat(value.replace("Z", "+00:00"))


def suite_health() -> SuiteHealth:
    """ Wall time of the test step in the last successful CI run on master (raises if gh can't say) """
    runs = _gh_api(f"repos/{{owner}}/{{repo}}/actions/workflows/{CI_WORKFLOW}/runs?branch=master&status=success&per_page=1")
    run = runs["workflow_runs"][0]
    jobs = _gh_api(f"repos/{{owner}}/{{repo}}/actions/runs/{run['id']}/jobs")
    for job in jobs["jobs"]:
        for step in job["steps"]:
            if step["name"] == CI_TEST_STEP and step.get("started_at") and step.get("completed_at"):
                elapsed = _parse_github_time(step["completed_at"]) - _parse_github_time(step["started_at"])
                return SuiteHealth(seconds=int(elapsed.total_seconds()), run_id=run["id"],
                                   head_sha=run["head_sha"][:9], started_at=step["started_at"])
    raise LookupError(f"run {run['id']} has no timed '{CI_TEST_STEP}' step")


def deprecation_health() -> tuple[DeprecationHealth, dict[str, str]]:
    report = find_deprecations()
    in_repo = len(report.runtime_in_repo)
    return DeprecationHealth(external=len(report.external), ours=len(report.ours), import_in_repo=in_repo,
                             import_in_packages=len(report.runtime) - in_repo), report.unavailable


def gather_health(ci: bool = True, deprecations: bool = True) -> Health:
    unavailable = {}
    coverage = docstring_coverage()
    docs = check_docs()
    lint = lint_health()
    if lint is None:
        unavailable["lint"] = "no lint.txt - run scripts/linting/run_pylint.sh"
    deprecated = None
    if deprecations:
        deprecated, missing = deprecation_health()
        if "static" in missing:
            deprecated = None
            unavailable["deprecations"] = missing["static"]
    else:
        unavailable["deprecations"] = "skipped (--no-deprecations)"
    suite = None
    if ci:
        try:
            suite = suite_health()
        except (OSError, subprocess.SubprocessError, LookupError, KeyError, IndexError, ValueError) as e:
            unavailable["suite"] = f"gh could not report the last CI run ({type(e).__name__})"
    else:
        unavailable["suite"] = "skipped (--no-ci)"
    return Health(sha=git_short_sha(),
                  documented_modules=sum(d for d, _ in coverage.values()),
                  counted_modules=sum(t for _, t in coverage.values()),
                  dead_citations=len(docs.dead), docs_checked=docs.docs_checked, dead_in_plans=len(docs.plan_dead),
                  commands_by_category=command_categories(), lint=lint, deprecations=deprecated, suite=suite, unavailable=unavailable)


def _minutes(seconds: int) -> str:
    return f"{seconds // 60}m{seconds % 60:02d}s"


def render_health(health: Health) -> str:
    total_commands = sum(health.commands_by_category.values())
    coverage = health.documented_modules / health.counted_modules if health.counted_modules else 0
    lines = [f"vg health at {health.sha}",
             f"  module docstrings  {coverage:.0%} ({health.documented_modules}/{health.counted_modules} modules)",
             f"  dead citations     {health.dead_citations} in {health.docs_checked} docs"
             + (f" (+{health.dead_in_plans} in plans)" if health.dead_in_plans else ""),
             f"  one-off commands   {health.one_off_commands} of {total_commands} ("
             + ", ".join(f"{c} {n}" for c, n in health.commands_by_category.items()) + ")"]
    if lint := health.lint:
        detail = ", ".join(f"{t} {n}" for t, n in lint.by_type.items())
        score = f", {lint.score}/10" if lint.score else ""
        lines.append(f"  pylint messages    {lint.messages} ({detail}{score}; lint.txt {lint.age_days} day(s) old)")
        lines.append(f"  pylint errors      {len(lint.errors)}")
        lines.extend(f"    {error}" for error in lint.errors[:LISTED_LINT_ERRORS])
        if len(lint.errors) > LISTED_LINT_ERRORS:
            lines.append(f"    ... and {len(lint.errors) - LISTED_LINT_ERRORS} more in lint.txt")
    else:
        lines.append(f"  pylint messages    unavailable: {health.unavailable['lint']}")
    if deprecated := health.deprecations:
        lines.append(f"  deprecated calls   {deprecated.external} stdlib/package + {deprecated.ours} of our own @deprecated"
                     f"; {deprecated.import_in_repo} warned at import (+{deprecated.import_in_packages} in packages)")
    else:
        lines.append(f"  deprecated calls   unavailable: {health.unavailable['deprecations']}")
    if suite := health.suite:
        lines.append(f"  test suite (CI)    {_minutes(suite.seconds)} (run {suite.run_id} on {suite.head_sha}, {suite.started_at[:10]})")
    else:
        lines.append(f"  test suite (CI)    unavailable: {health.unavailable['suite']}")
    return "\n".join(lines)
