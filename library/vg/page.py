"""
`vg page`: render one URL through the Django test client against the live database and describe the
result compactly — status, templates, an accessible-text outline, and (with queries=True) the number
of queries production would run plus the repeated ones (N+1 candidates).

The request runs as the `claude_agent` user (a plain member of all_users) inside a transaction that
is always rolled back, so a view with a write side-effect cannot dirty the database. Templates are
tracked the way Django's test runner does it: Template._render is swapped for the instrumented one
for the duration of the request.
"""
import re
from collections import Counter
from contextlib import contextmanager
from dataclasses import dataclass, field

from bs4 import BeautifulSoup
from django.contrib.auth.models import User
from django.db import connection, transaction
from django.template import Template
from django.test import Client
from django.test.signals import template_rendered
from django.test.utils import CaptureQueriesContext, instrumented_test_render
from django.urls import NoReverseMatch, reverse

from library.django_utils.unittest_utils import _normalize_sql, production_query_count
from library.guardian_utils import all_users_group

AGENT_USERNAME = "claude_agent"
_WHITESPACE_RE = re.compile(r"\s+")


@dataclass
class PageResult:
    url: str
    status: int
    templates: list[str]
    outline: list[str]
    text: str
    html: str
    links: list[tuple[str, str]]
    forms: list[str]
    redirect_to: str | None = None
    production_queries: int | None = None
    total_queries: int | None = None
    repeated: list[tuple[str, int]] = field(default_factory=list)


class AgentUserMissing(Exception):
    pass


def resolve_url(url_or_name: str, kwargs: dict | None = None) -> str:
    """ A path is used as given; anything else is treated as a URL name (namespace:name allowed) """
    if url_or_name.startswith("/"):
        return url_or_name
    try:
        return reverse(url_or_name, kwargs=kwargs or None)
    except NoReverseMatch as e:
        raise ValueError(f"'{url_or_name}' is neither a path nor a resolvable URL name "
                         f"(pass --kwargs key=value for its arguments): {e}") from e


def get_agent_user(username: str = AGENT_USERNAME) -> User:
    try:
        return User.objects.get(username=username)
    except User.DoesNotExist as e:
        raise AgentUserMissing(f"User '{username}' does not exist. Create the agent user with "
                               f"`manage.py vg page --create-user`, or pass --as <existing user>") from e


def create_agent_user() -> tuple[User, bool]:
    """ The agent user: member of all_users, no lab, cannot log in with a password """
    user, created = User.objects.get_or_create(username=AGENT_USERNAME,
                                               defaults={"first_name": "Claude", "last_name": "Agent"})
    if created:
        user.set_unusable_password()
        user.save()
    user.groups.add(all_users_group())
    return user, created


@contextmanager
def _track_templates(templates: list[str]):
    def on_render(sender, template, context, **kwargs):  # pylint: disable=unused-argument
        if template.name and template.name not in templates:
            templates.append(template.name)

    original_render = Template._render  # pylint: disable=protected-access
    Template._render = instrumented_test_render  # pylint: disable=protected-access
    template_rendered.connect(on_render)
    try:
        yield
    finally:
        template_rendered.disconnect(on_render)
        Template._render = original_render  # pylint: disable=protected-access


def _outline(soup: BeautifulSoup) -> tuple[list[str], list[tuple[str, str]], list[str]]:
    outline, links, forms = [], [], []
    title = soup.title.get_text(strip=True) if soup.title else ""
    if title:
        outline.append(f"title {title}")
    body = soup.body or soup
    for element in body.find_all(["h1", "h2", "h3", "h4", "table", "form", "a"]):
        name = element.name
        if name.startswith("h"):
            text = _WHITESPACE_RE.sub(" ", element.get_text(" ", strip=True))
            if text:
                outline.append(f"{name} {text}")
        elif name == "table":
            rows = element.find_all("tr")
            headers = [_WHITESPACE_RE.sub(" ", th.get_text(" ", strip=True)) for th in element.find_all("th")]
            width = max((len(r.find_all(["td", "th"])) for r in rows), default=len(headers))
            label = element.get("id") or element.get("data-datatable-url") or ""
            summary = f"table {len(rows)}x{width}"
            if label:
                summary += f" {label}"
            if headers:
                summary += " (" + ", ".join(h for h in headers if h)[:200] + ")"
            outline.append(summary)
        elif name == "form":
            fields = [f.get("name") for f in element.find_all(["input", "select", "textarea"])
                      if f.get("name") and f.get("type") != "hidden"]
            form_name = element.get("id") or element.get("name") or element.get("action") or "(unnamed)"
            method = (element.get("method") or "get").upper()
            line = f"form {form_name} [{method}]: {', '.join(dict.fromkeys(fields)) or '(no fields)'}"
            outline.append(line)
            forms.append(line)
        elif name == "a":
            href = element.get("href")
            text = _WHITESPACE_RE.sub(" ", element.get_text(" ", strip=True))
            if href and not href.startswith(("#", "javascript:")):
                links.append((text, href))
    return outline, links, forms


def render_page(url_or_name: str, username: str | None = None, kwargs: dict | None = None,
                queries: bool = False) -> PageResult:
    url = resolve_url(url_or_name, kwargs)
    user = get_agent_user(username or AGENT_USERNAME)
    client = Client(SERVER_NAME="localhost")
    templates: list[str] = []

    with transaction.atomic(), _track_templates(templates), CaptureQueriesContext(connection) as ctx:
        client.force_login(user)  # inside the transaction: login() writes last_login
        login_queries = len(ctx.captured_queries)
        response = client.get(url, follow=False)
        transaction.set_rollback(True)

    content = response.content.decode(response.charset or "utf-8", errors="replace") if hasattr(response, "content") else ""
    is_html = "html" in (response.get("Content-Type") or "")
    soup = BeautifulSoup(content, "lxml") if is_html and content else BeautifulSoup("", "lxml")
    outline, links, forms = _outline(soup) if is_html else ([], [], [])
    text = _WHITESPACE_RE.sub(" ", soup.get_text(" ", strip=True)) if is_html else content

    result = PageResult(url=url, status=response.status_code, templates=templates, outline=outline,
                        text=text, html=content, links=links, forms=forms,
                        redirect_to=response.get("Location") if 300 <= response.status_code < 400 else None)
    if queries:
        captured = ctx.captured_queries[login_queries:]
        sql = [q["sql"] for q in captured if not q["sql"].startswith(("SAVEPOINT", "RELEASE SAVEPOINT", "ROLLBACK TO"))]
        result.total_queries = len(sql)
        result.production_queries = production_query_count(captured)
        result.repeated = [(normalized, n) for normalized, n in Counter(_normalize_sql(s) for s in sql).most_common() if n > 1]
    return result

